"""Generate validated conformers for explicitly specified chemical states (Torsional Diffusion + ETKDG baseline).

Single state:
  python -m jmcg.condgen.run --model_dir <TD model dir> --out <dir> --state "CC(=O)[O-]" --parent "CC(=O)O" \
      --n_confs 8 --seed 0
Batch (csv with columns id,parent_smiles,state_smiles[,note]):
  python -m jmcg.condgen.run --model_dir <dir> --out <dir> --states_csv states.csv --n_confs 8 --seed 0

Outputs in --out:
  <id>__<method>.sdf     one record per conformer; SD properties carry state, method, seed, checkpoint, validity
  <id>.json              state description, torsion bookkeeping (parent heavy-atom keyed), generation settings,
                         per-conformer validation, diversity, failures
  summary.csv            one row per (state, method)
  run_info.json          checkpoint (path + SHA-256), TD repo, versions, settings, UTC/IST timestamps
A state that cannot be represented (parse error, not a protomer/tautomer of the parent, atom type outside the
checkpoint vocabulary, embedding failure) is recorded as a failure with its reason; nothing is silently substituted.
"""
import argparse
import datetime as dt
import json
import os
import platform
import time

import numpy as np
import pandas as pd
import rdkit
import torch
from rdkit import Chem

from jmcg.condgen import states as st
from jmcg.condgen import validate as va

IST = dt.timezone(dt.timedelta(hours=5, minutes=30))


def stamp():
    now = dt.datetime.now(dt.timezone.utc)
    return f"{now:%Y-%m-%d %H:%M:%S} UTC / {now.astimezone(IST):%Y-%m-%d %H:%M:%S} IST"


def write_sdf(path, mols, props):
    w = Chem.SDWriter(path)
    w.SetKekulize(True)
    for k, (m, p) in enumerate(zip(mols, props)):
        m = Chem.Mol(m)
        m.SetProp('_Name', f"{p['state_id']}|{p['method']}|conf{k}")
        for key, val in p.items():
            m.SetProp(key, json.dumps(val) if isinstance(val, (dict, list)) else str(val))
        w.write(m)
    w.close()


def run_state(gen, row, a):
    sid = row['id']
    rec = dict(state_id=sid, note=row.get('note', ''), seed=a.seed, n_confs_requested=a.n_confs, failures=[])
    try:
        mol, info = st.build_state(row['state_smiles'], row.get('parent_smiles') or None)
    except st.StateError as e:
        rec['category'] = 'D'
        rec['failures'].append(dict(stage='state', reason=str(e)))
        return rec, []
    rec['state'] = info
    target = info['state_smiles_canonical']
    hard, ood = gen.compatibility(mol) if gen else (['no TD checkpoint loaded'], [])
    rec['checkpoint_compatibility'] = dict(hard=hard, out_of_distribution=ood)
    rows = []
    has_ez = any(x['kind'] == 'double_bond' for x in info['stereo_elements'])
    methods = ['etkdg_v3'] + ([] if hard else ['td'] + (['td_lockdb'] if a.lock_arm and has_ez else []))
    # A in-distribution neutral | B charged, inside vocabulary and size | C outside training size/charge buckets |
    # D unsupported atom types (TD) - D states still get the ETKDG baseline
    rec['category'] = 'D' if hard else 'C' if gen.size_or_bucket_ood(mol) else 'B' if info['charged_atoms'] else 'A'
    for h in hard:
        rec['failures'].append(dict(stage='checkpoint', method='td', reason=h))
    for method in methods:
        t0 = time.time()
        try:
            if method == 'etkdg_v3':
                from jmcg.condgen.td_backend import etkdg_baseline
                mols, ginfo = etkdg_baseline(mol, a.n_confs, a.seed), {}
            else:
                mols, ginfo = gen.generate(mol, a.n_confs, a.seed, lock_stereo_double_bonds=method == 'td_lockdb')
                ginfo['torsions_parent_keyed'] = st.torsions_parent_keyed(
                    ginfo.pop('edge_index'), ginfo.pop('edge_mask'), info['heavy_map_state_to_parent'])
        except Exception as e:  # noqa: BLE001 - record and continue with the next method/state
            rec['failures'].append(dict(stage='generate', method=method, reason=f'{type(e).__name__}: {e}'[:300]))
            continue
        secs = time.time() - t0
        bounds = va.rdDistGeom.GetMoleculeBoundsMatrix(mols[0])
        checks = [va.check_conformer(m, target, bounds) for m in mols]
        div = va.diversity([m for m, c in zip(mols, checks) if c['finite']])
        summ = va.summarise(checks)
        rec[method] = dict(generation=ginfo, seconds=round(secs, 2), n_generated=len(mols), validation=checks,
                           summary=summ, diversity=div)
        props = [dict(state_id=sid, method=method, label='BASELINE RDKit ETKDGv3' if method == 'etkdg_v3' else
                      'Torsional Diffusion (frozen checkpoint)', state_smiles=target, state_class=info['state_class'],
                      net_charge=info['net_charge'], seed=a.seed, conf_index=k,
                      checkpoint_sha256='' if method == 'etkdg_v3' else gen.ckpt_sha256,
                      valid=c['valid'], stereo_ok=c.get('stereo'), clash_free=c.get('clash'),
                      mmff94s_energy=c.get('mmff_energy')) for k, c in enumerate(checks)]
        write_sdf(os.path.join(a.out, f'{sid}__{method}.sdf'), mols, props)
        rows.append(dict(state_id=sid, category=rec['category'], method=method, state_class=info['state_class'],
                         net_charge=info['net_charge'],
                         n_heavy=info['n_heavy'], n_rot_td=ginfo.get('n_rotatable_td'),
                         locked=len(ginfo.get('locked_double_bonds', [])),
                         td_model_applied=ginfo.get('model_applied'), n=len(mols), seconds=round(secs, 2),
                         ood='; '.join(ood) if method != 'etkdg_v3' else '', **summ,
                         **{k: v for k, v in div.items() if k != 'n'}))
    return rec, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model_dir', default=None, help='TD model dir (model_parameters.yml + ckpt); omit = ETKDG only')
    ap.add_argument('--ckpt', default='best_model.pt')
    ap.add_argument('--out', required=True)
    ap.add_argument('--state', default=None)
    ap.add_argument('--parent', default=None)
    ap.add_argument('--id', default='state0')
    ap.add_argument('--states_csv', default=None)
    ap.add_argument('--n_confs', type=int, default=8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--steps', type=int, default=20)
    ap.add_argument('--stereo_policy', choices=['require', 'enumerate'], default='require',
                    help='require: reject states with unassigned stereo; enumerate: expand them into one explicit '
                         'state per stereoisomer (ids <id>__st<k>). No stereoisomer is ever picked silently.')
    ap.add_argument('--lock_arm', action='store_true', help='also run TD with stereo double bonds removed from torsions')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', '4')))

    if a.states_csv:
        table = pd.read_csv(a.states_csv, comment='#').fillna('').to_dict('records')
    else:
        table = [dict(id=a.id, parent_smiles=a.parent or '', state_smiles=a.state)]
    if a.stereo_policy == 'enumerate':
        expanded = []
        for row in table:
            try:
                isomers = st.expand_stereo(row['state_smiles'])
            except st.StateError:
                isomers = [row['state_smiles']]  # build_state reports the parse error
            for k, smi in enumerate(isomers):
                expanded.append(dict(row, id=row['id'] if len(isomers) == 1 else f"{row['id']}__st{k}", state_smiles=smi,
                                     note=f"{row.get('note', '')} [stereo-enumerated {k + 1}/{len(isomers)}]"
                                     if len(isomers) > 1 else row.get('note', '')))
        table = expanded
    gen = None
    if a.model_dir:
        from jmcg.condgen.td_backend import TDGenerator
        gen = TDGenerator(a.model_dir, a.ckpt, steps=a.steps, device='cpu')
    info = dict(started=stamp(), args=vars(a), host=platform.node(), python=platform.python_version(),
                rdkit=rdkit.__version__, torch=torch.__version__, generator=gen.describe() if gen else None)
    print(json.dumps(info, indent=1), flush=True)

    summary = []
    for row in table:
        rec, rows = run_state(gen, row, a)
        with open(os.path.join(a.out, f"{row['id']}.json"), 'w') as f:
            json.dump(rec, f, indent=1, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
        summary += rows
        fails = '; '.join(f"{x.get('method', x['stage'])}: {x['reason']}" for x in rec['failures'])
        print(f"{row['id']:<24s} {rec.get('state', {}).get('state_class', '-'):<20s} "
              + ' '.join(f"{r['method']}: valid {r['n_valid']}/{r['n']} clusters {r.get('n_clusters_0p5A')}" for r in rows)
              + (f'  FAIL[{fails}]' if fails else ''), flush=True)
    pd.DataFrame(summary).to_csv(os.path.join(a.out, 'summary.csv'), index=False)
    info['finished'] = stamp()
    with open(os.path.join(a.out, 'run_info.json'), 'w') as f:
        json.dump(info, f, indent=1)
    print('done', info['finished'])


if __name__ == '__main__':
    main()
