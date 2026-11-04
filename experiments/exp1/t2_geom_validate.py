"""EXP1-V2-T2 task 3.4-3.5 / task 5: geometry validation and error decomposition of every E1 arm (read only).

Per generated conformer of each arm (same TD seeds, same molecules):
  identity      connectivity re-perceived from the 3D coordinates (rdDetermineBonds.DetermineConnectivity) equals the
                molecular graph
  stereo        isomeric SMILES with stereo assigned from 3D equals the corrected SMILES; and (refined arms) equals that
                of the TD input conformer, i.e. refinement changed no stereo assignment
  torsions      not recomputed here: W drift is logged by t1_e1_readout.py with TD's own rotatable-bond mask
  clashes       minimum distance between atoms separated by >= 3 bonds (heavy atoms) and the count below 2.0 Å
  local error   vs the reference conformer nearest to it in heavy-atom RMSD (the RMSD matrix column minimum written
                by evaluate_confs.py), on heavy atoms mapped by substructure match:
                  bond-length MAE (Å), bond-angle MAE (deg), ring-dihedral MAE (deg, central bond in a ring),
                  acyclic-dihedral MAE (deg, central bond not in a ring; includes the W torsions TD sampled)
Molecule-level means are then averaged over the paired molecules (the 184 of the E1 compare), with a paired bootstrap of
arm − TD and arm − xTB over molecules.

Run (TD venv): python experiments/exp1/t2_geom_validate.py --in_dir <n200> --t1_dir <T1_e1> --true_mols <test_mols.pkl>
                   --seed 0 --workers 16 --out <json> [--arm KEY=PREFIX ...]
"""
from __future__ import annotations

import csv
import json
import pickle
from argparse import ArgumentParser
from multiprocessing import Pool

import numpy as np

ARMS = ('td', 'Lprojfix', 'Lnoprojfix', 'xtb', 'mmff', 'etkdg')
_G = {}


def _clean(smi, confs):
    """evaluate_confs.py clean_confs: reference conformers whose graph equals the corrected SMILES."""
    from rdkit import Chem
    s = Chem.MolToSmiles(Chem.MolFromSmiles(smi, sanitize=False), isomericSmiles=False)
    return [c for c in confs if Chem.MolToSmiles(Chem.RemoveHs(c, sanitize=False), isomericSmiles=False) == s]


def _dihedral(x, q):
    b0, b1, b2 = x[q[:, 0]] - x[q[:, 1]], x[q[:, 2]] - x[q[:, 1]], x[q[:, 3]] - x[q[:, 2]]
    b1n = b1 / np.linalg.norm(b1, axis=1, keepdims=True)
    v = b0 - (b0 * b1n).sum(1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(1, keepdims=True) * b1n
    return np.degrees(np.arctan2((np.cross(b1n, v) * w).sum(1), (v * w).sum(1)))


def _wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def _internals(mol):
    """Heavy-atom bonds, angles, ring dihedrals (central bond in a ring) and acyclic dihedrals, as index arrays."""
    bonds = np.array([[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mol.GetBonds()]).reshape(-1, 2)
    nb = [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]
    angles = np.array([[i, j, k] for j in range(mol.GetNumAtoms()) for ii, i in enumerate(nb[j]) for k in nb[j][ii + 1:]]
                      ).reshape(-1, 3)
    ring, acyc = [], []
    for b in mol.GetBonds():
        j, k = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        for i in nb[j]:
            for l in nb[k]:
                if i != k and l != j and i != l:
                    (ring if b.IsInRing() else acyc).append([i, j, k, l])
    return bonds, angles, np.array(ring).reshape(-1, 4), np.array(acyc).reshape(-1, 4)


def _graph_ok(mol):
    """Connectivity re-perceived from the 3D coordinates equals the molecular graph (all atoms)."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    ref = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    m2 = Chem.RWMol(mol)
    for b in list(m2.GetBonds()):
        m2.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    try:
        rdDetermineBonds.DetermineConnectivity(m2)
    except Exception:  # noqa: BLE001 - perception failure counts as not ok
        return False
    return ref == {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in m2.GetBonds()}


def _skeleton(mol):
    """Graph with every bond single and nothing aromatic or charged, so reference and generated heavy-atom graphs match
    by connectivity and element whatever their bond-order perception."""
    from rdkit import Chem
    m = Chem.RWMol(mol)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    for at in m.GetAtoms():
        at.SetIsAromatic(False); at.SetFormalCharge(0); at.SetNoImplicit(True); at.SetNumExplicitHs(0)
    return m.GetMol()


def _stereo_smiles(mol):
    """Canonical isomeric SMILES with tetrahedral and double-bond stereo assigned from the 3D coordinates."""
    from rdkit import Chem
    m = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(m)
    return Chem.MolToSmiles(Chem.RemoveHs(m))


def _mol(job):
    """All arms of one molecule. job = (csmi, {arm: [mols]}, [reference mols], {arm: rmsd matrix}, td_input [mols])."""
    from rdkit import Chem
    csmi, arms, refs, rm, td_in = job
    out = {}
    try:
        refH = [Chem.RemoveHs(r) for r in refs]
        target = Chem.MolToSmiles(Chem.MolFromSmiles(csmi))
        td_stereo = [_stereo_smiles(c) for c in td_in]
        for arm, confs in arms.items():
            probe = Chem.RemoveHs(confs[0])
            bonds, angles, rd, ad = _internals(probe)
            far = np.argwhere(np.triu(Chem.GetDistanceMatrix(probe) >= 3))
            sk = _skeleton(probe)
            maps = [_skeleton(r).GetSubstructMatch(sk) for r in refH]
            m_rmsd = rm[arm]
            rows = []
            for j, c in enumerate(confs):
                x = Chem.RemoveHs(c).GetConformer().GetPositions()
                st = _stereo_smiles(c)
                r = dict(identity=_graph_ok(c), stereo_vs_smiles=st == target,
                         stereo_vs_td_input=(st == td_stereo[j]) if arm != 'etkdg' else None)
                d = np.linalg.norm(x[far[:, 0]] - x[far[:, 1]], axis=1) if len(far) else np.array([np.inf])
                r['min_nonbonded'] = float(d.min()); r['clashes_lt2'] = int((d < 2.0).sum())
                col = m_rmsd[:, j]
                k = int(np.nanargmin(col)) if not np.isnan(col).all() else None
                if k is not None and len(maps[k]) == len(x):
                    y = refH[k].GetConformer().GetPositions()[list(maps[k])]

                    def bl(p):
                        return np.linalg.norm(p[bonds[:, 0]] - p[bonds[:, 1]], axis=1)

                    def ang(p):
                        u, w = p[angles[:, 0]] - p[angles[:, 1]], p[angles[:, 2]] - p[angles[:, 1]]
                        cos = (u * w).sum(1) / np.linalg.norm(u, axis=1) / np.linalg.norm(w, axis=1)
                        return np.degrees(np.arccos(np.clip(cos, -1, 1)))
                    r['bond_mae'] = float(np.abs(bl(x) - bl(y)).mean()) if len(bonds) else np.nan
                    r['angle_mae'] = float(np.abs(ang(x) - ang(y)).mean()) if len(angles) else np.nan
                    r['ring_dih_mae'] = float(np.abs(_wrap(_dihedral(x, rd) - _dihedral(y, rd))).mean()) if len(rd) else np.nan
                    r['acyclic_dih_mae'] = float(np.abs(_wrap(_dihedral(x, ad) - _dihedral(y, ad))).mean()) if len(ad) else np.nan
                    r['nearest_rmsd'] = float(col[k])
                rows.append(r)
            out[arm] = rows
    except Exception as e:  # noqa: BLE001 - reported per molecule
        return csmi, {'error': f'{type(e).__name__}: {e}'}
    return csmi, out


KEYS = ('bond_mae', 'angle_mae', 'ring_dih_mae', 'acyclic_dih_mae', 'nearest_rmsd', 'min_nonbonded')


def main():
    ap = ArgumentParser()
    ap.add_argument('--in_dir', required=True)
    ap.add_argument('--t1_dir', required=True)
    ap.add_argument('--true_mols', required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--arm', action='append', default=[], metavar='KEY=PREFIX',
                    help='extra arm read from PREFIX.pkl and PREFIX.eval.pkl (EXP1-V2-T2 final evaluation)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    s = a.seed
    extra = dict(x.split('=', 1) for x in a.arm)
    arms = ARMS + tuple(extra)
    subset = {r['corrected_smiles']: r['smiles'] for r in csv.DictReader(open(f'{a.in_dir}/subset.csv'))}
    true = pickle.load(open(a.true_mols, 'rb'))
    paths = {'td': (f'{a.in_dir}/td{s}.pkl', f'{a.in_dir}/td{s}.eval.pkl'),
             'etkdg': (f'{a.in_dir}/etkdg{s}.pkl', f'{a.in_dir}/etkdg{s}.eval.pkl'),
             **{k: (f'{a.t1_dir}/td{s}_{k}.pkl', f'{a.t1_dir}/td{s}_{k}.eval.pkl') for k in ARMS if k not in ('td', 'etkdg')},
             **{k: (f'{v}.pkl', f'{v}.eval.pkl') for k, v in extra.items()}}
    confs = {k: pickle.load(open(p, 'rb')) for k, (p, _) in paths.items()}
    ev = {k: {cs: r['rmsd'] for (_, cs), r in pickle.load(open(e, 'rb'))['results'].items()} for k, (_, e) in paths.items()}
    paired = sorted(set.intersection(*[{cs for cs, m in ev[k].items() if not np.isnan(m).any()} for k in arms]))
    jobs = [(cs, {k: confs[k][cs] for k in arms}, _clean(cs, true[subset[cs]]), {k: ev[k][cs] for k in arms}, confs['td'][cs])
            for cs in paired]
    with Pool(a.workers) as pool:
        res = dict(pool.map(_mol, jobs, chunksize=2))
    errors = {cs: v['error'] for cs, v in res.items() if 'error' in v}
    ok = [cs for cs in paired if cs not in errors]
    per = {k: {} for k in arms}
    conf_level = {k: dict(conformers=0, identity_fail=0, stereo_vs_smiles_fail=0, stereo_changed_vs_td_input=0,
                          clash_conformers=0, no_reference_match=0) for k in arms}
    for cs in ok:
        for k in arms:
            rows = res[cs][k]
            c = conf_level[k]
            c['conformers'] += len(rows)
            c['identity_fail'] += sum(not r['identity'] for r in rows)
            c['stereo_vs_smiles_fail'] += sum(not r['stereo_vs_smiles'] for r in rows)
            c['stereo_changed_vs_td_input'] += sum(r['stereo_vs_td_input'] is False for r in rows)
            c['clash_conformers'] += sum(r['clashes_lt2'] > 0 for r in rows)
            c['no_reference_match'] += sum('bond_mae' not in r for r in rows)
            per[k][cs] = {m: float(np.nanmean([r.get(m, np.nan) for r in rows])) for m in KEYS}
    rng = np.random.default_rng(0)
    summ = dict(seed=s, paired=len(paired), molecules_ok=len(ok), molecule_errors=errors, conformer_checks=conf_level,
                means={k: {m: float(np.nanmean([per[k][cs][m] for cs in ok])) for m in KEYS} for k in arms}, paired_diff={})
    for ref in ('td', 'xtb'):
        for k in arms:
            if k == ref:
                continue
            for m in KEYS:
                d = np.array([per[k][cs][m] - per[ref][cs][m] for cs in ok])
                d = d[~np.isnan(d)]
                bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(5000)]
                summ['paired_diff'][f'{k}-{ref}:{m}'] = dict(mean=float(d.mean()), n=int(len(d)),
                                                            ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])
    json.dump(dict(summary=summ, per_molecule=per), open(a.out, 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in summ.items() if k != 'paired_diff'}, indent=1, default=float))
    for key, v in summ['paired_diff'].items():
        print(key, round(v['mean'], 4), [round(x, 4) for x in v['ci95']])


if __name__ == '__main__':
    main()
