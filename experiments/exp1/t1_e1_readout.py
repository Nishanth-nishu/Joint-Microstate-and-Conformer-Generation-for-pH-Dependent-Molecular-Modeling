"""EXP1-V2-T1 E1 arms: post-process existing frozen-TD conformers (design: research/experiments/EXP1_V2_train/DESIGN.md).

  --method l     τ-frozen read-out of a trained V2 geometry head (``flexitors.sampling.l_readout``; ``--project off``
                 for the projector ablation), 20 steps, torsion restoration and guard, seed = the TD conformer
  --method mmff  non-learned control: MMFF94s minimisation (≤ 200 iterations), then the W dihedrals restored to TD's
  --method xtb   non-learned control: GFN2-xTB optimisation (``xtb`` on PATH, one thread), every W dihedral constrained
                 to TD's value (harmonic, force constant ``--xtb_fc`` Eh), then the W dihedrals restored exactly

Writes ``--out`` (dict smiles → list of RDKit mols, the evaluate_confs.py input format, same molecules and conformer
count as the TD input) and ``--out``.log.json (per conformer: guard counters, W-dihedral drift vs TD, seconds, errors).
A conformer whose processing raises is kept unchanged and logged as an error (no silent drop).

Run (TD venv): python experiments/exp1/t1_e1_readout.py --confs td0.pkl --method l --run_dir <T1 dir> --project on \
                   --out td0_Lproj.pkl --workers 12
"""
from __future__ import annotations

import json
import os
import pickle
import subprocess
import tempfile
import time
import traceback
from argparse import ArgumentParser, Namespace
from multiprocessing import Pool

import numpy as np

_G = {}


def _init(a):
    _G['a'] = a
    try:
        _load(a)
    except Exception as e:  # noqa: BLE001 - a raising Pool initializer respawns workers forever; fail per conformer
        _G['init_error'] = f'init {type(e).__name__}: {e}'


def _load(a):
    import torch
    import yaml
    torch.set_num_threads(1)
    if a.method == 'l':
        from flexitors.models.flexitors_score import build_model
        args = Namespace(**yaml.load(open(f'{a.run_dir}/model_parameters.yml'), Loader=yaml.Loader))
        m = build_model(args).eval()
        m.load_state_dict(torch.load(f'{a.run_dir}/{a.checkpoint}', map_location='cpu'), strict=True)
        _G['model'], _G['sigma_min'] = m, args.sigma_min


def _parts(mol):
    import torch
    from flexitors.geometry.l_step import Guard, w_quadruples
    from flexitors.td import import_baseline
    feat = import_baseline('utils.featurization')
    tors = import_baseline('utils.torsion')
    data = feat.featurize_mol(mol, 'qm9')
    em, mr = tors.get_transformation_mask(data)
    data.edge_mask, data.mask_rotate = torch.tensor(em), mr
    z = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    ei = data.edge_index.numpy()
    nb = [[] for _ in z]
    for u, v in ei[:, ::2].T:
        nb[u].append(int(v)); nb[v].append(int(u))
    x = mol.GetConformer().GetPositions().astype(np.float64)
    quads, kept, _, _ = w_quadruples(z, ei.T[np.asarray(em, dtype=bool)], nb, x)
    sides = np.asarray(mr, dtype=bool)[kept] if len(quads) else np.zeros((0, len(z)), dtype=bool)
    return data, x, quads, sides, Guard(mol, ei[:, ::2].T, x)


def _mmff(mol, x, quads, sides):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from flexitors.geometry.internal import dihedral
    from flexitors.geometry.l_step import restore_dihedrals
    m = Chem.Mol(mol)
    props = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
    ff = AllChem.MMFFGetMoleculeForceField(m, props)
    ff.Minimize(maxIts=200)
    y = np.array(ff.Positions()).reshape(-1, 3)
    _G['raw'] = y
    return restore_dihedrals(y, quads, sides, dihedral(x, quads)) if len(quads) else y


def _xtb(mol, x, quads, sides, fc, timeout):
    """Returns (positions, converged). Atoms in the RDKit order; xtb's constraint indices are 1-based."""
    from rdkit import Chem
    from flexitors.geometry.internal import dihedral
    from flexitors.geometry.l_step import restore_dihedrals
    chrg = Chem.GetFormalCharge(mol)
    uhf = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    with tempfile.TemporaryDirectory(prefix='t1xtb_') as d:
        with open(f'{d}/m.xyz', 'w') as f:
            f.write(f'{len(x)}\n\n' + ''.join(f'{a.GetSymbol()} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n'
                                              for a, p in zip(mol.GetAtoms(), x)))
        with open(f'{d}/c.inp', 'w') as f:
            f.write(f'$constrain\n force constant={fc}\n')
            for q, v in zip(quads, np.degrees(dihedral(x, quads)) if len(quads) else []):
                f.write(' dihedral: ' + ','.join(str(int(i) + 1) for i in q) + f',{v:.6f}\n')
            f.write('$end\n')
        r = subprocess.run(['xtb', 'm.xyz', '--gfn', '2', '--opt', 'normal', '--chrg', str(chrg), '--uhf', str(uhf),
                            '--input', 'c.inp'], cwd=d, capture_output=True, text=True, timeout=timeout,
                           env=dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_STACKSIZE='1G'))
        out = f'{d}/xtbopt.xyz' if os.path.exists(f'{d}/xtbopt.xyz') else f'{d}/xtblast.xyz'
        if not os.path.exists(out):
            raise RuntimeError(f'xtb exit {r.returncode}: {r.stdout[-300:]}')
        y = np.array([list(map(float, line.split()[1:4])) for line in open(out).read().splitlines()[2:2 + len(x)]])
        converged = r.returncode == 0 and 'GEOMETRY OPTIMIZATION CONVERGED' in r.stdout
    _G['raw'] = y
    return (restore_dihedrals(y, quads, sides, dihedral(x, quads)) if len(quads) else y), converged


def _one(job):
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from flexitors.geometry.internal import dihedral, wrap_angle
    smi, j, mol = job
    a = _G['a']
    t0 = time.time()
    rec = dict(smiles=smi, conf=j)
    out = Chem.Mol(mol)
    try:
        if 'init_error' in _G:
            raise RuntimeError(_G['init_error'])
        _G.pop('raw', None)
        data, x, quads, sides, guard = _parts(mol)
        rec['n_quads'] = int(len(quads))
        if a.method == 'l':
            from flexitors.sampling.l_readout import l_readout
            y, c = l_readout(_G['model'], data, x, quads, sides, guard, _G['sigma_min'], project=a.project == 'on')
            rec.update(c)
        elif a.method == 'xtb':
            y, rec['xtb_converged'] = _xtb(mol, x, quads, sides, a.xtb_fc, a.xtb_timeout)
        else:
            y = _mmff(mol, x, quads, sides)
        if 'raw' in _G and len(quads):  # before the exact restore: shows whether the method kept the torsions itself
            rec['w_drift_raw_deg'] = float(np.degrees(np.abs(wrap_angle(dihedral(_G['raw'], quads) - dihedral(x, quads)))).max())
        geo, chi = guard.check(y)
        rec.update(final_geometry_ok=bool(geo), final_chirality_ok=bool(chi),
                   w_drift_deg=float(np.degrees(np.abs(wrap_angle(dihedral(y, quads) - dihedral(x, quads))).max()))
                   if len(quads) else 0.0,
                   heavy_displacement_rms=float(np.sqrt(((y - x) ** 2).sum(1)[[at.GetAtomicNum() > 1 for at in
                                                                               mol.GetAtoms()]].mean())))
        conf = out.GetConformer()
        for i, p in enumerate(y):
            conf.SetAtomPosition(i, Point3D(*map(float, p)))
    except Exception as e:  # noqa: BLE001 - kept unchanged and logged
        rec['error'] = f'{type(e).__name__}: {e}'[:300]
        rec['traceback'] = traceback.format_exc()[-800:]
    rec['seconds'] = time.time() - t0
    return smi, j, out, rec


def main():
    ap = ArgumentParser()
    ap.add_argument('--confs', required=True)
    ap.add_argument('--method', choices=('l', 'mmff', 'xtb'), required=True)
    ap.add_argument('--xtb_fc', type=float, default=1.0, help='xtb dihedral-constraint force constant (Eh)')
    ap.add_argument('--xtb_timeout', type=float, default=600.0, help='seconds per conformer')
    ap.add_argument('--run_dir')
    ap.add_argument('--checkpoint', default='best_model.pt')
    ap.add_argument('--project', choices=('on', 'off'), default='on')
    ap.add_argument('--limit_mols', type=int, default=0, help='dry run on the first N molecules (pipeline check only)')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    confs = pickle.load(open(a.confs, 'rb'))
    smis = list(confs)[:a.limit_mols] if a.limit_mols else list(confs)
    jobs = [(s, j, m) for s in smis for j, m in enumerate(confs[s])]
    t0 = time.time()
    with Pool(a.workers, initializer=_init, initargs=(a,)) as pool:
        res = pool.map(_one, jobs, chunksize=4)
    out = {s: [None] * len(confs[s]) for s in smis}
    log = []
    for s, j, m, rec in res:
        out[s][j] = m
        log.append(rec)
    pickle.dump(out, open(a.out, 'wb'))
    ok = [r for r in log if 'error' not in r]
    summ = dict(method=a.method, project=a.project if a.method == 'l' else None, run_dir=a.run_dir,
                checkpoint=a.checkpoint, molecules=len(smis), conformers=len(log), errors=len(log) - len(ok),
                wall_seconds=time.time() - t0, workers=a.workers,
                cpu_seconds_per_conformer=float(np.mean([r['seconds'] for r in log])),
                max_w_drift_deg=max([r['w_drift_deg'] for r in ok], default=None),
                final_geometry_fail=sum(not r['final_geometry_ok'] for r in ok),
                final_chirality_fail=sum(not r['final_chirality_ok'] for r in ok),
                median_heavy_displacement_rms=float(np.median([r['heavy_displacement_rms'] for r in ok])) if ok else None)
    if a.method == 'l':
        summ.update(share_conformers_with_skip=float(np.mean([r['skips'] > 0 for r in ok])) if ok else None,
                    mean_skips=float(np.mean([r['skips'] for r in ok])) if ok else None,
                    mean_halvings=float(np.mean([r['halvings'] for r in ok])) if ok else None,
                    chirality_rejections=int(sum(r['chirality_rejections'] for r in ok)))
    if a.method == 'xtb':
        summ.update(xtb_not_converged=sum(not r['xtb_converged'] for r in ok))
    if a.method in ('xtb', 'mmff'):
        raw = [r['w_drift_raw_deg'] for r in ok if 'w_drift_raw_deg' in r]
        summ.update(w_drift_raw_deg_median=float(np.median(raw)) if raw else None,
                    w_drift_raw_deg_p95=float(np.percentile(raw, 95)) if raw else None,
                    w_drift_raw_deg_max=float(max(raw)) if raw else None)
    json.dump(dict(summary=summ, conformers=log), open(a.out + '.log.json', 'w'), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float))


if __name__ == '__main__':
    main()
