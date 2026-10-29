"""State-conditioned conformer-ensemble descriptors phi(X_z) for the joint prototype, and the 2D control psi(z).
Definitions are fixed in research/experiments/E2_joint/PREREG.md.

phi (from X_z ~ p(x|z): ETKDGv3, n conformers with per-conformer seeds, MMFF94s minimisation, every conformer
validated with jmcg.condgen.validate.check_conformer; Boltzmann weights at 298.15 K on MMFF94s energies relative to the
state's own minimum, so no cross-state energy comparison is made):
  phi1 = -ln N_eff, N_eff = exp(-sum w ln w)         (dimensionless)
  phi2 = <n_IMHB>  polar H (on N/O) .. N/O acceptor in [1.5, 2.5] A, acceptor >= 3 bonds from the donor heavy atom
  phi3 = <polar SASA> / 100 A^2  (N, O and H on N/O; FreeSASA, probe 1.4 A)
psi (2D only): psi1 = rotatable bonds, psi2 = topological donor-acceptor pairs 4-7 bonds apart, psi3 = TPSA / 100 A^2.
"""
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdDistGeom, rdFreeSASA, rdMolDescriptors

from jmcg.condgen import validate as va

KT_298 = 0.0019872043 * 298.15  # kcal/mol (MMFF energies are kcal/mol)
DONOR_ACCEPTOR = ('N', 'O')


def conformer_seeds(seed, n):
    """Same per-conformer RDKit seeds as jmcg.condgen.td_backend.conformer_seeds (that module needs the TD venv)."""
    return [int(x) for x in np.random.default_rng(seed).integers(1, 2 ** 31 - 1, n)]


def _imhb_count(m, conf_id, dist_path):
    pos = m.GetConformer(conf_id).GetPositions()
    n = 0
    for h in m.GetAtoms():
        if h.GetSymbol() != 'H':
            continue
        d = h.GetNeighbors()[0]
        if d.GetSymbol() not in DONOR_ACCEPTOR:
            continue
        for acc in m.GetAtoms():
            if acc.GetSymbol() in DONOR_ACCEPTOR and acc.GetIdx() != d.GetIdx() and dist_path[d.GetIdx(), acc.GetIdx()] >= 3:
                r = float(np.linalg.norm(pos[h.GetIdx()] - pos[acc.GetIdx()]))
                if 1.5 <= r <= 2.5:
                    n += 1
    return n


def _polar_sasa(m, conf_id):
    radii = rdFreeSASA.classifyAtoms(m)
    rdFreeSASA.CalcSASA(m, radii, confIdx=conf_id)
    tot = 0.0
    for a in m.GetAtoms():
        polar = a.GetSymbol() in DONOR_ACCEPTOR or (a.GetSymbol() == 'H' and a.GetNeighbors()[0].GetSymbol() in DONOR_ACCEPTOR)
        if polar:
            tot += float(a.GetProp('SASA'))
    return tot


def ensemble_features(smiles, n_confs=30, seed=0, max_iters=2000, return_ensemble=False):
    """dict(phi1, phi2, phi3, n_valid, n_generated, failures) for one microstate SMILES. return_ensemble=True adds
    'ensemble': (Mol holding the minimised conformers, valid conformer ids, Boltzmann weights) for the sampler."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    target = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    out = dict(n_generated=0, n_valid=0, failures=[])
    if not AllChem.MMFFHasAllMoleculeParams(m):
        out['failures'].append('MMFF94s parameters missing')
        return out
    cids = []
    for s in conformer_seeds(seed, n_confs):
        p = AllChem.ETKDGv3()
        p.randomSeed = s
        t = Chem.Mol(m)
        if AllChem.EmbedMolecule(t, p) == 0:
            cids.append(m.AddConformer(t.GetConformer(), assignId=True))
    out['n_generated'] = len(cids)
    if not cids:
        out['failures'].append('ETKDGv3 embedding failed for every seed')
        return out
    bounds = rdDistGeom.GetMoleculeBoundsMatrix(m)
    dist_path = Chem.GetDistanceMatrix(m)
    # MMFF typing rewrites bond types to MMFF aromaticity in place (e.g. 2-pyridone ring -> Kekule single/double), so
    # it gets its own copy; validation and phi use the RDKit-perceived state with the minimised coordinates.
    mf = Chem.Mol(m)
    props = AllChem.MMFFGetMoleculeProperties(mf, mmffVariant='MMFF94s')
    E, keep = [], []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mf, props, confId=cid)
        ff.Minimize(maxIts=max_iters)
        conf = m.GetConformer(cid)
        for i, xyz in enumerate(mf.GetConformer(cid).GetPositions()):
            conf.SetAtomPosition(i, xyz.tolist())
        single = Chem.Mol(m, confId=cid)
        chk = va.check_conformer(single, target, bounds)
        if chk['valid']:
            E.append(ff.CalcEnergy())
            keep.append(cid)
        else:
            out['failures'].append('invalid after MMFF: ' + ','.join(k for k in va.CHECKS[1:] if chk.get(k) is False))
    out['n_valid'] = len(keep)
    if not keep:
        return out
    E = np.array(E)
    w = np.exp(-(E - E.min()) / KT_298)
    w /= w.sum()
    H = float(-(w[w > 0] * np.log(w[w > 0])).sum())
    out['phi1'] = -H  # -ln N_eff
    out['phi2'] = float(sum(wi * _imhb_count(m, c, dist_path) for wi, c in zip(w, keep)))
    out['phi3'] = float(sum(wi * _polar_sasa(m, c) for wi, c in zip(w, keep))) / 100.0
    out['E_range_kcal'] = float(E.max() - E.min())
    if return_ensemble:
        out['ensemble'] = (m, keep, w)
    return out


def graph_features(smiles):
    m = Chem.MolFromSmiles(smiles)
    dm = Chem.GetDistanceMatrix(m)
    don = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() in DONOR_ACCEPTOR and a.GetTotalNumHs() > 0]
    acc = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() in DONOR_ACCEPTOR]
    pairs = sum(1 for d in don for a in acc if a != d and 4 <= dm[d, a] <= 7)
    return dict(psi1=float(rdMolDescriptors.CalcNumRotatableBonds(m)), psi2=float(pairs),
                psi3=float(Descriptors.TPSA(m)) / 100.0)
