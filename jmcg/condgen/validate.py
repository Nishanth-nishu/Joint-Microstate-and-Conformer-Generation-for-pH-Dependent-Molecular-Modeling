"""Per-conformer validation and set diversity for generated conformers of one explicit state.

Checks (each conformer; thresholds are this project's choices, stated here, not taken from a benchmark):
  finite        all coordinates finite
  roundtrip     MolBlock (Kekule, as written to SDF) -> sanitised Mol keeps, atom by atom and in the same order:
                element, formal charge, H neighbours; bond by bond: bond type (aromaticity re-perceived); and the
                canonical isomeric SMILES with stereo read from the 3D coordinates equals the target state SMILES
  connectivity  bonds re-derived from 3D only (rdDetermineBonds.DetermineConnectivity, covalent radii) == state bonds
  bonds_angles  every 1-2 distance within the RDKit distance-geometry bounds +-0.10 A and every 1-3 distance within the
                bounds +-0.20 A (bounds from rdDistGeom.GetMoleculeBoundsMatrix: UFF-based ideal lengths/angles for
                the state's own bond orders/hybridisation). TD keeps the ETKDG seed's bonds and angles exactly, so for TD
                this verifies the seed, not the model.
  clash         no atom pair >= 4 bonds apart closer than 0.7 x sum of Bondi vdW radii (H-H: 0.6 x) - a lenient
                severe-overlap screen, not a quality metric. Exception: a polar H (bonded to N/O) and an N/O acceptor
                form a hydrogen bond, whose H...acceptor distance (~1.6-2.2 A) is below 0.7 x vdW sum (1.90 A for
                H...O); such pairs are counted as n_hbond_contacts and only flagged as a clash below 1.5 A
  stereo        stereo perceived from 3D (tetrahedral + double bond) == target, all specified elements
  mmff_energy   MMFF94s single point (kcal/mol) of the unrelaxed geometry: only comparable between conformers of the
                SAME state and method-to-method on that state; never across states (different atoms/charges)
Diversity (set): symmetry-aware heavy-atom pairwise RMSD (rdMolAlign.GetBestRMS), Butina clusters at 0.5 A and 0.1 A
(the 0.1 A count is reported as "unique heavy-atom conformers"; H-only differences such as OH/CH3 rotation are not
counted).
"""
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdMolAlign
from rdkit.ML.Cluster import Butina

try:
    from rdkit.Chem import rdDetermineBonds
except ImportError:  # RDKit < 2022.09
    rdDetermineBonds = None

PT = Chem.GetPeriodicTable()
TOL12, TOL13, CLASH, CLASH_HH, HBOND_MIN = 0.10, 0.20, 0.7, 0.6, 1.5


def _bonds(m):
    return {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in m.GetBonds()}


_RHP = Chem.RemoveHsParameters()
_RHP.removeDefiningBondStereo = True  # else an imine =N-H stays explicit and the SMILES never matches the target


def heavy(m):
    """Heavy-atom Mol. E/Z of a double bond defined only by an H (C=N-H imines) cannot be written without that H and
    is therefore NOT validated."""
    return Chem.RemoveHs(m, _RHP)


def canon_from_3d(m):
    c = Chem.Mol(m)
    Chem.AssignStereochemistryFrom3D(c)
    return Chem.MolToSmiles(heavy(c), isomericSmiles=True)


def roundtrip(m, target_smiles):
    rt = Chem.MolFromMolBlock(Chem.MolToMolBlock(m, kekulize=True), removeHs=False, sanitize=True)
    if rt is None or rt.GetNumAtoms() != m.GetNumAtoms():
        return dict(rt_atoms=False, rt_bonds=False, rt_isomeric=False)
    atoms = all(a.GetAtomicNum() == b.GetAtomicNum() and a.GetFormalCharge() == b.GetFormalCharge()
                and a.GetTotalNumHs(includeNeighbors=True) == b.GetTotalNumHs(includeNeighbors=True)
                for a, b in zip(m.GetAtoms(), rt.GetAtoms()))
    bt = lambda x: {(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))), b.GetBondType()) for b in x.GetBonds()}
    iso = Chem.MolToSmiles(heavy(rt), isomericSmiles=True)
    return dict(rt_atoms=atoms, rt_bonds=bt(m) == bt(rt), rt_isomeric=iso == target_smiles, rt_smiles=iso)


def check_conformer(m, target_smiles, bounds=None):
    pos = m.GetConformer().GetPositions()
    r = dict(finite=bool(np.isfinite(pos).all()))
    if not r['finite']:
        return dict(r, valid=False)
    r.update(roundtrip(m, target_smiles))
    r['connectivity'] = None
    if rdDetermineBonds is not None:
        g = Chem.RWMol(m)
        for b in list(g.GetBonds()):
            g.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        try:
            rdDetermineBonds.DetermineConnectivity(g, charge=Chem.GetFormalCharge(m))
            r['connectivity'] = _bonds(g) == _bonds(m)
        except Exception as e:  # noqa: BLE001
            r['connectivity_error'] = str(e)[:80]
    topo = Chem.GetDistanceMatrix(m)
    d = np.linalg.norm(pos[:, None] - pos[None], axis=-1)
    bm = rdDistGeom.GetMoleculeBoundsMatrix(m) if bounds is None else bounds
    iu = np.triu_indices(len(pos), 1)
    lb, ub = bm.T[iu], bm[iu]  # RDKit: upper bounds in the upper triangle, lower bounds in the lower triangle
    dist, tp = d[iu], topo[iu]
    dev = np.maximum(lb - dist, dist - ub)
    r['max_12_dev'] = float(dev[tp == 1].max()) if (tp == 1).any() else 0.0
    r['max_13_dev'] = float(dev[tp == 2].max()) if (tp == 2).any() else 0.0
    r['bonds_angles'] = r['max_12_dev'] <= TOL12 and r['max_13_dev'] <= TOL13
    z = np.array([a.GetAtomicNum() for a in m.GetAtoms()])
    rv = np.array([PT.GetRvdw(int(x)) for x in z])
    f = np.where((z[:, None] == 1) & (z[None] == 1), CLASH_HH, CLASH)
    bad = (topo >= 4) & (d < f * (rv[:, None] + rv[None]))
    # polar H (on N/O) ... N/O acceptor: a hydrogen bond, not an overlap -> only flagged below HBOND_MIN
    polar_h = np.array([a.GetAtomicNum() == 1 and any(n.GetAtomicNum() in (7, 8) for n in a.GetNeighbors())
                        for a in m.GetAtoms()])
    acc = np.isin(z, (7, 8))
    hb = (polar_h[:, None] & acc[None]) | (acc[:, None] & polar_h[None])
    r['n_hbond_contacts'] = int(np.triu(bad & hb & (d >= HBOND_MIN), 1).sum())
    bad &= ~hb | (d < HBOND_MIN)
    r['n_clash_pairs'] = int(np.triu(bad, 1).sum())
    r['clash_pairs'] = [(int(i), int(j), m.GetAtomWithIdx(int(i)).GetSymbol() + m.GetAtomWithIdx(int(j)).GetSymbol(),
                         round(float(d[i, j]), 3)) for i, j in zip(*np.nonzero(np.triu(bad, 1)))][:5]
    r['min_nonlocal_dist'] = float(d[topo >= 4].min()) if (topo >= 4).any() else None
    r['clash'] = r['n_clash_pairs'] == 0
    r['stereo_smiles_3d'] = canon_from_3d(m)
    r['stereo'] = r['stereo_smiles_3d'] == target_smiles
    try:
        props = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
        r['mmff_energy'] = float(AllChem.MMFFGetMoleculeForceField(m, props).CalcEnergy()) if props else None
    except Exception:  # noqa: BLE001
        r['mmff_energy'] = None
    r['valid'] = bool(r['finite'] and r['rt_atoms'] and r['rt_bonds'] and r['rt_isomeric'] and r['connectivity']
                      and r['bonds_angles'] and r['clash'] and r['stereo'])
    return r


def diversity(mols):
    n = len(mols)
    if n < 2:
        return dict(n=n)
    heavy = [Chem.RemoveHs(m) for m in mols]
    dm, full = [], np.zeros((n, n))
    for i in range(1, n):
        for j in range(i):
            v = rdMolAlign.GetBestRMS(Chem.Mol(heavy[i]), Chem.Mol(heavy[j]), maxMatches=1000)
            dm.append(v); full[i, j] = full[j, i] = v
    iu = np.triu_indices(n, 1)
    return dict(n=n, pairwise_rmsd_min=float(full[iu].min()), pairwise_rmsd_median=float(np.median(full[iu])),
                pairwise_rmsd_max=float(full[iu].max()),
                n_clusters_0p5A=len(Butina.ClusterData(dm, n, 0.5, isDistData=True, reordering=True)),
                n_unique_0p1A=len(Butina.ClusterData(dm, n, 0.1, isDistData=True, reordering=True)))


CHECKS = ['valid', 'finite', 'rt_atoms', 'rt_bonds', 'rt_isomeric', 'connectivity', 'bonds_angles', 'clash', 'stereo']


def summarise(rows):
    out = {f'n_{k}': int(sum(r.get(k) is True for r in rows)) for k in CHECKS}
    out['n_connectivity_unchecked'] = int(sum(r.get('connectivity') is None for r in rows))
    e = [r['mmff_energy'] for r in rows if r.get('mmff_energy') is not None]
    out['mmff_median'] = float(np.median(e)) if e else None
    out['mmff_min'] = float(np.min(e)) if e else None
    return out
