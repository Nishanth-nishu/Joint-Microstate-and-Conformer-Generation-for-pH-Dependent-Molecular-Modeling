"""EXP-1 M0 (revised, §10.3): oracle-velocity integration test of the V2 geometry sampler. CPU, no learning.

Question: can the sampler's projected, torsion-restored, constraint-guarded L update follow an ORACLE field to the
pure-L target over large, nonlinear ring moves? No model is involved; this gates any FlexiTors training.

Per (molecule, conformer k) of the VALIDATION paired cache with pair_ok:
  x_R = pos[k] (matched RDKit seed), x_G = gt_pos[k] (its own GT, Kabsch-aligned offline);
  x_G* = x_G with its W dihedrals set to x_R's by rigid side rotations, then Kabsch-aligned to x_R on heavy atoms
         (§10.1 item 3: x_G* − x_R carries no W-dihedral difference);
  N = 20 steps x ← x + (1/(N − n))·P_W(x)(x_G* − x), each with torsion restoration and the §2.6/§2.8 guard.
Metrics: recovery = 1 − median_pairs RMSD_heavy(x_N, x_G*) / median_pairs RMSD_heavy(x_R, x_G*); ring-dihedral RMSD to
x_G* (rings ≥ 4 atoms; share of pairs > 10°); W-dihedral drift; chirality rejections; halvings and skips.
Population (primary): molecules with ≥ 1 heavy rotatable bond AND a ring of ≥ 4 atoms. Sanity rows: rigid
(n_rot_heavy = 0) and m = 0. KILL if recovery < 50 % or > 5 % of primary pairs end with a skipped step or a flipped
centre. Test-set pairs are not used here.

Run (TD venv, needs torch_geometric to unpickle the cache):
  python experiments/exp1/m0_oracle.py --cache data/QM9/cache_paired.val --out <dir> [--procs 16 --limit 0]
"""
from __future__ import annotations

import json
import os
import pickle
import time
from argparse import ArgumentParser
from multiprocessing import Pool

import numpy as np
import pandas as pd

from flexitors.geometry.internal import dihedral, wrap_angle
from flexitors.geometry.l_step import Guard, kabsch_to, projected_oracle_flow, restore_dihedrals, w_quadruples
from flexitors.geometry.torsion_types import heavy_rotatable_bonds

DATA = None
N_STEPS = 20
KILL_RECOVERY, KILL_FAIL_SHARE = 0.50, 0.05


def rmsd_aligned(x: np.ndarray, y: np.ndarray, idx: np.ndarray) -> float:
    return float(np.sqrt(((kabsch_to(x, y, idx)[idx] - y[idx]) ** 2).sum(1).mean()))


def ring_quads(mol) -> np.ndarray:
    q = []
    for r in mol.GetRingInfo().AtomRings():
        if len(r) >= 4:
            q += [(r[i], r[(i + 1) % len(r)], r[(i + 2) % len(r)], r[(i + 3) % len(r)]) for i in range(len(r))]
    return np.array(q, dtype=np.int64).reshape(-1, 4)


def ring_rmsd_deg(x, y, rq) -> float | None:
    if not len(rq):
        return None
    return float(np.degrees(np.sqrt(np.mean(wrap_angle(dihedral(x, rq) - dihedral(y, rq)) ** 2))))


def run_molecule(i: int, restore=True, max_confs=None) -> list[dict]:
    d = DATA[i]
    mol = d.mol
    z = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    heavy = np.where(z > 1)[0]
    ei = d.edge_index.numpy()
    bonds = ei[:, ::2].T
    nbrs = [[] for _ in range(len(z))]
    for u, v in bonds:
        nbrs[u].append(int(v)); nbrs[v].append(int(u))
    em = np.asarray(d.edge_mask, dtype=bool)
    rot_edges = ei.T[em]
    sides = np.asarray(d.mask_rotate, dtype=bool)
    n_rot_heavy = len(heavy_rotatable_bonds(mol))
    ring4 = any(len(r) >= 4 for r in mol.GetRingInfo().AtomRings())
    stratum = ('m0' if len(rot_edges) == 0 else 'rigid') if n_rot_heavy == 0 else ('primary' if ring4 else 'flexible_no_ring4')
    rq = ring_quads(mol)
    rows = []
    ks = range(len(d.pos)) if max_confs is None else range(min(max_confs, len(d.pos)))
    for k in ks:
        base = dict(mol_idx=i, smiles=d.canonical_smi, conf=k, stratum=stratum, n_atoms=len(z), m_td=len(rot_edges),
                    n_rot_heavy=n_rot_heavy, pair_ok=bool(d.pair_ok[k]))
        if not d.pair_ok[k]:
            rows.append(base)
            continue
        xR = d.pos[k].numpy().astype(np.float64)
        xG = d.gt_pos[k].numpy().astype(np.float64)
        quads, kept, dropped, n_h = w_quadruples(z, rot_edges, nbrs, xR)
        sq = sides[kept]
        tau0 = dihedral(xR, quads) if len(quads) else np.zeros(0)
        xGs = restore_dihedrals(xG, quads, sq, tau0) if len(quads) else xG.copy()
        xGs = kabsch_to(xGs, xR, heavy)
        guard = Guard(mol, bonds, xR)
        t0 = time.time()
        xN, info = projected_oracle_flow(xR, xGs, quads, sq, guard, n_steps=N_STEPS, restore=restore)
        rows.append(dict(base, w_rows=len(quads), w_dropped=dropped, w_h_terminated=n_h,
                         n_stereo_centres=len(guard.centres), target_chirality_differs=not guard.chirality_ok(xGs),
                         rmsd_R_Gs=rmsd_aligned(xR, xGs, heavy), rmsd_N_Gs=rmsd_aligned(xN, xGs, heavy),
                         rmsd_R_G=rmsd_aligned(xR, xG, heavy), ring_R_Gs_deg=ring_rmsd_deg(xR, xGs, rq),
                         ring_N_Gs_deg=ring_rmsd_deg(xN, xGs, rq), seconds=time.time() - t0, **info))
    return rows


def _run(args):
    return run_molecule(*args)


def summarise(df: pd.DataFrame) -> dict:
    out = {}
    ok = df[df.pair_ok]
    for st, g in [('all', ok)] + list(ok.groupby('stratum')):
        if not len(g):
            continue
        fail = g.skips.gt(0) | g.flipped_centre_at_end.astype(bool)
        r = g.ring_N_Gs_deg.dropna()
        r0 = g.ring_R_Gs_deg.dropna()
        out[st] = dict(n_pairs=len(g), n_molecules=int(g.mol_idx.nunique()),
                       median_rmsd_R_Gs=float(g.rmsd_R_Gs.median()), median_rmsd_N_Gs=float(g.rmsd_N_Gs.median()),
                       recovery=float(1 - g.rmsd_N_Gs.median() / g.rmsd_R_Gs.median()) if g.rmsd_R_Gs.median() > 0 else None,
                       mean_rmsd_R_Gs=float(g.rmsd_R_Gs.mean()), mean_rmsd_N_Gs=float(g.rmsd_N_Gs.mean()),
                       share_pairs_with_skip=float(g.skips.gt(0).mean()),
                       share_pairs_flipped_at_end=float(g.flipped_centre_at_end.astype(bool).mean()),
                       share_fail_kill_rule=float(fail.mean()),
                       share_target_chirality_differs=float(g.target_chirality_differs.astype(bool).mean()),
                       share_fail_given_target_chirality_same=float(fail[~g.target_chirality_differs.astype(bool)].mean())
                       if (~g.target_chirality_differs.astype(bool)).any() else None,
                       mean_halvings_per_pair=float(g.halvings.mean()), chirality_rejections=int(g.chirality_rejections.sum()),
                       w_drift_final_max=float(g.final_w_drift.max()),
                       w_drift_before_restore_median=float(g.max_drift_before_restore.median()),
                       w_drift_before_restore_p95=float(g.max_drift_before_restore.quantile(0.95)),
                       w_rows_dropped=int(g.w_dropped.sum()), w_rows=int(g.w_rows.sum()),
                       w_rows_h_terminated=int(g.w_h_terminated.sum()),
                       ring_pairs=len(r), ring_share_gt10deg_R=float((r0 > 10).mean()) if len(r0) else None,
                       ring_share_gt10deg_N=float((r > 10).mean()) if len(r) else None,
                       ring_median_deg_R=float(r0.median()) if len(r0) else None,
                       ring_median_deg_N=float(r.median()) if len(r) else None,
                       seconds_per_pair_median=float(g.seconds.median()))
    return out


def main():
    global DATA
    ap = ArgumentParser()
    ap.add_argument('--cache', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--procs', type=int, default=16)
    ap.add_argument('--limit', type=int, default=0, help='first N molecules only (smoke)')
    ap.add_argument('--ablation_pairs', type=int, default=300, help='primary molecules re-run with restore off')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    DATA = pickle.load(open(a.cache, 'rb'))
    idx = list(range(len(DATA) if not a.limit else min(a.limit, len(DATA))))
    t0 = time.time()
    with Pool(a.procs) as p:
        rows = [r for rs in p.map(_run, [(i, True, None) for i in idx], chunksize=8) for r in rs]
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, 'm0_pairs.csv'), index=False)
    S = summarise(df)
    prim = S.get('primary', {})
    verdict = 'NO PRIMARY PAIRS'
    if prim:
        kill = []
        if prim['recovery'] < KILL_RECOVERY:
            kill.append(f"recovery {prim['recovery']:.3f} < {KILL_RECOVERY}")
        if prim['share_fail_kill_rule'] > KILL_FAIL_SHARE:
            kill.append(f"{prim['share_fail_kill_rule']:.3f} of pairs end with a skipped step or flipped centre > {KILL_FAIL_SHARE}")
        verdict = 'KILL: ' + '; '.join(kill) if kill else 'PASS'
    # restore_torsions = off, measured once (EXP-1 §10.2) on the first primary molecules, first conformer each
    pm = sorted(df[(df.stratum == 'primary') & df.pair_ok].mol_idx.unique())[:a.ablation_pairs]
    with Pool(a.procs) as p:
        abl = pd.DataFrame([r for rs in p.map(_run, [(i, False, 1) for i in pm], chunksize=4) for r in rs])
    on = df[df.mol_idx.isin(pm) & (df.conf == 0)]
    res = dict(cache=a.cache, n_molecules=len(idx), n_pairs_total=len(df), n_pairs_pair_ok=int(df.pair_ok.sum()),
               n_steps=N_STEPS, kill_rule=dict(recovery_min=KILL_RECOVERY, fail_share_max=KILL_FAIL_SHARE),
               verdict=verdict, strata=S, seconds=time.time() - t0,
               ablation_restore_off=dict(n_pairs=len(abl), summary=summarise(abl).get('primary'),
                                         same_pairs_restore_on=summarise(on).get('primary')))
    json.dump(res, open(os.path.join(a.out, 'm0_summary.json'), 'w'), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
