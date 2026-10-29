"""Sequential baseline, stage 1 (states): parent molecule + pH -> enumerated microstates -> Uni-pKa free energies ->
model populations p(z | G, pH) and macro-pKa -> selected states for the conformer stage (jmcg.condgen.run, TD venv).

Run (jmcg venv):
  python -m jmcg.seq.state_stage --parents parents.csv --ph 7.4 --model_dir <unipka models> --out <dir>
     [--ph_grid 2,4,6,7.4,8,10,12] [--p_min 0.05] [--top_k 3] [--threads 8]
parents.csv: id,smiles (any protonation form; it is the enumeration seed and the parent for state validation).
Outputs: ensembles.csv, enumeration_failures.csv, energies.csv (state x fold), populations.csv (parent x pH x state,
  fold mean and sd), macro_pka.csv, selected_states.csv (condgen format: id,parent_smiles,state_smiles,note), run_info.json.
Selection at the primary --ph: states with fold-mean population >= p_min, at most top_k, always at least the most
populated one; the covered population mass is recorded. The selection is a sampling budget, not a claim that the
remaining states are absent.
"""
import argparse
import datetime as dt
import json
import os
import platform

import numpy as np
import pandas as pd
import rdkit

from jmcg.seq import enumerate as en
from jmcg.seq import score as sc

IST = dt.timezone(dt.timedelta(hours=5, minutes=30))


def stamp():
    now = dt.datetime.now(dt.timezone.utc)
    return f"{now:%Y-%m-%d %H:%M:%S} UTC / {now.astimezone(IST):%Y-%m-%d %H:%M:%S} IST"


def enumerate_parents(parents, lower=-2, upper=2):
    ens, fails = [], []
    for pid, smi in parents:
        try:
            states, f = en.enumerate_ensemble(smi, lower=lower, upper=upper)
        except Exception as e:  # noqa: BLE001 - recorded per parent
            fails.append(dict(parent=pid, smiles=smi, level='', reason=f'{type(e).__name__}: {e}'[:300]))
            continue
        ens += [dict(parent=pid, parent_smiles=smi, **s) for s in states]
        fails += [dict(parent=pid, **x) for x in f]
    return pd.DataFrame(ens), pd.DataFrame(fails, columns=['parent', 'smiles', 'level', 'reason'])


def per_parent(ens, energies, ph_values):
    """populations (fold mean/sd) and macro-pKa per parent."""
    E = energies.pivot_table(index='smiles', columns='fold', values='G')
    folds = list(E.columns)
    pops, pkas = [], []
    for pid, g in ens.groupby('parent'):
        G = E.loc[g.smiles].to_numpy()  # states x folds
        q = g.charge.to_numpy()
        for ph in ph_values:
            P = np.stack([sc.populations(G[:, k], q, ph) for k in range(len(folds))], 1)
            for (_, r), pm, ps in zip(g.iterrows(), P.mean(1), P.std(1)):
                pops.append(dict(parent=pid, pH=ph, state_id=r.state_id, smiles=r.smiles, charge=r.charge,
                                 p_mean=pm, p_sd_folds=ps))
        per_fold = [sc.macro_pkas(G[:, k], q) for k in range(len(folds))]
        for key in per_fold[0]:
            v = [d[key] for d in per_fold]
            pkas.append(dict(parent=pid, q_acid=key[0], q_base=key[1], pKa_mean=float(np.mean(v)), pKa_sd_folds=float(np.std(v))))
    return pd.DataFrame(pops), pd.DataFrame(pkas)


def select(pops, ph, p_min, top_k):
    rows = []
    for pid, g in pops[pops.pH == ph].groupby('parent'):
        g = g.sort_values('p_mean', ascending=False)
        keep = g[g.p_mean >= p_min].head(top_k)
        if keep.empty:
            keep = g.head(1)
        cov = float(keep.p_mean.sum())
        for _, r in keep.iterrows():
            rows.append(dict(parent=pid, state_id=r.state_id, smiles=r.smiles, charge=r.charge, p_mean=r.p_mean,
                             covered_mass=cov))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--parents', required=True)
    ap.add_argument('--ph', type=float, required=True)
    ap.add_argument('--ph_grid', default='')
    ap.add_argument('--model_dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--p_min', type=float, default=0.05)
    ap.add_argument('--top_k', type=int, default=3)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--lower', type=int, default=-2)
    ap.add_argument('--upper', type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    info = dict(started=stamp(), args=vars(a), host=platform.node(), python=platform.python_version(),
                rdkit=rdkit.__version__, enumerator=en.provenance(), PKA_MEAN=sc.PKA_MEAN,
                population_convention='p(z|G,pH) = softmax_z(-G_z - q_z ln10 (pH - PKA_MEAN)) over the enumerated '
                                      'ensemble; per fold, then averaged')
    par = pd.read_csv(a.parents)
    ens, fails = enumerate_parents(list(zip(par['id'], par['smiles'])), a.lower, a.upper)
    ens.to_csv(os.path.join(a.out, 'ensembles.csv'), index=False)
    fails.to_csv(os.path.join(a.out, 'enumeration_failures.csv'), index=False)
    energies = sc.free_energies(ens.smiles.tolist(), a.model_dir, threads=a.threads, task='seqstage')
    energies.to_csv(os.path.join(a.out, 'energies.csv'), index=False)
    grid = sorted({a.ph} | {float(x) for x in a.ph_grid.split(',') if x})
    pops, pkas = per_parent(ens, energies, grid)
    pops.to_csv(os.path.join(a.out, 'populations.csv'), index=False)
    pkas.to_csv(os.path.join(a.out, 'macro_pka.csv'), index=False)
    sel = select(pops, a.ph, a.p_min, a.top_k)
    psmi = dict(zip(par['id'], par['smiles']))
    pd.DataFrame(dict(id=[f'{p}__{s}' for p, s in zip(sel.parent, sel.state_id)],
                      parent_smiles=[psmi[p] for p in sel.parent], state_smiles=sel.smiles,
                      note=[f'pH {a.ph}: p={pm:.3f} (covered {c:.3f})' for pm, c in zip(sel.p_mean, sel.covered_mass)])
                 ).to_csv(os.path.join(a.out, 'selected_states.csv'), index=False)
    sel.to_csv(os.path.join(a.out, 'selection_detail.csv'), index=False)
    info.update(finished=stamp(), n_parents=len(par), n_states=len(ens), n_enumeration_failures=len(fails),
                n_selected=len(sel))
    json.dump(info, open(os.path.join(a.out, 'run_info.json'), 'w'), indent=1, default=str)
    print(json.dumps({k: info[k] for k in ('n_parents', 'n_states', 'n_enumeration_failures', 'n_selected')}))


if __name__ == '__main__':
    main()
