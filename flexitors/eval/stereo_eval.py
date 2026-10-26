"""Exp 3 endpoints per molecule (EXP-3 §4, §9.3; R-CHEM conformance note, round5/IMPLEMENTATION.md, items 1-2).

* **E/Z fidelity** (DRUGS primary): share of generated conformers whose geometric label equals the SMILES target on
  every specified structural bond (cos(τ_b − φ_b) > 0 on TD's quadruple; ``round5/data/ez_targets.csv``).
* **GT-ensemble consistency** (QM9 primary, the 37 structural molecules): share of generated conformers whose joint
  cis/trans state over the structural bonds occurs in the GT ensemble. GT atom order is not the generated one, so each
  GT conformer is mapped onto TD's molecule by graph matching (atomic numbers and connectivity only, so GEOM bond
  orders and charges do not matter); every match counts, so symmetric reference atoms cannot fake a mismatch
  (plan_P2 V-E3-5).
* **Planarity**: share of (conformer, bond) pairs with |δ_b| ≤ 20° to the nearest of {0, π}, over EVERY TD-rotatable
  double bond (structural, non-stereo and imine N–H), i.e. the "planarity on the 102" (§9.3), and the mean deviation.
Molecules whose GT ensemble mixes states are flagged (EXP-3 §4: "3 of 38").
"""
from __future__ import annotations

import pickle
from argparse import ArgumentParser

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from flexitors.data.ez_targets import load_targets
from flexitors.data.qm9_val import clean_confs
from flexitors.eval.validity import PLANARITY_DEG, ez_state_array
from flexitors.geometry.internal import dihedral
from flexitors.geometry.torsion_types import is_structural_stereo, td_reference_quadruple, td_rotatable_bonds

RDLogger.DisableLog('rdApp.*')
MAX_MATCHES = 2000


def double_bonds(mh: Chem.Mol) -> list[tuple[int, tuple[int, int, int, int], bool]]:
    """(bond_idx, TD quadruple, structural?) for every TD-rotatable double bond of TD's molecule."""
    out = []
    for k in td_rotatable_bonds(mh):
        b = mh.GetBondWithIdx(k)
        if b.GetBondType() == Chem.BondType.DOUBLE:
            out.append((k, td_reference_quadruple(mh, k), is_structural_stereo(b)))
    return out


def _skeleton(mol: Chem.Mol) -> Chem.Mol:
    """Graph-only copy (all bonds single, no aromaticity or charges) for order-independent matching."""
    m = Chem.RWMol(mol)
    for a in m.GetAtoms():
        a.SetFormalCharge(0)
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    m.UpdatePropertyCache(strict=False)
    return m.GetMol()


def gt_states(mh: Chem.Mol, gts: list[Chem.Mol], quads: np.ndarray) -> set[tuple[bool, ...]]:
    """Set of joint cis states (cos τ > 0) of the GT conformers on ``quads`` (TD indices), over every graph match."""
    query = _skeleton(mh)
    states = set()
    for g in gts:
        matches = _skeleton(g).GetSubstructMatches(query, uniquify=False, useChirality=False, maxMatches=MAX_MATCHES)
        pos = g.GetConformer().GetPositions()
        for mt in matches:
            mapped = np.array(mt)[quads]
            states.add(tuple(bool(x) for x in np.cos(dihedral(pos, mapped)) > 0))
    return states


def molecule_endpoints(mh: Chem.Mol, gen: list[Chem.Mol], gts: list[Chem.Mol], targets: dict[int, float]) -> dict:
    """Per-molecule E/Z fidelity, GT-ensemble consistency and planarity of the generated conformers ``gen``."""
    dbs = double_bonds(mh)
    out = dict(n_double=len(dbs), n_structural=sum(s for _, _, s in dbs),
               n_specified=sum(k in targets for k, _, s in dbs if s), n_conf=len(gen))
    if not dbs:
        return out
    pos = np.stack([m.GetConformer().GetPositions() for m in gen])
    quads = np.array([q for _, q, _ in dbs])
    tau = dihedral(pos, quads)
    phi = np.array([targets.get(k, np.nan) if s else np.nan for k, _, s in dbs])
    match, dev = ez_state_array(tau, phi)
    out['planarity_share'] = float(np.mean(dev <= PLANARITY_DEG))
    out['planarity_mean_deg'] = float(dev.mean())
    spec = ~np.isnan(phi)
    if spec.any():
        out['ez_fidelity'] = float(np.mean(match[:, spec].all(axis=1)))
    struct = np.array([s for _, _, s in dbs])
    if struct.any() and gts:
        allowed = gt_states(mh, gts, quads[struct])
        joint = [tuple(bool(x) for x in row) for row in np.cos(tau[:, struct]) > 0]
        out['ez_gt_consistency'] = float(np.mean([j in allowed for j in joint])) if allowed else np.nan
        out['gt_mixed'] = len(allowed) > 1
    return out


def evaluate(confs: dict, csv: pd.DataFrame, true_mols: dict, targets: dict[str, dict[int, float]]) -> pd.DataFrame:
    """Per-molecule endpoints for every csv row that has generated conformers and a TD-rotatable double bond."""
    rows = []
    for r in csv.itertuples(index=False):
        gen = confs.get(r.corrected_smiles)
        m = Chem.MolFromSmiles(r.corrected_smiles)
        if not gen or m is None:
            continue
        mh = Chem.AddHs(m)
        if not double_bonds(mh):
            continue
        gts = clean_confs(r.corrected_smiles, true_mols.get(r.smiles, []))
        rows.append(dict(smiles=r.smiles, corrected_smiles=r.corrected_smiles,
                         **molecule_endpoints(mh, gen, gts, targets.get(r.corrected_smiles, {}))))
    return pd.DataFrame(rows)


def main() -> None:
    ap = ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--confs', required=True)
    ap.add_argument('--csv', required=True)
    ap.add_argument('--true_mols', required=True)
    ap.add_argument('--ez_table', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    with open(a.confs, 'rb') as f:
        confs = pickle.load(f)
    with open(a.true_mols, 'rb') as f:
        true_mols = pickle.load(f)
    D = evaluate(confs, pd.read_csv(a.csv), true_mols, load_targets(a.ez_table))
    D.to_csv(a.out, index=False)
    mean = {c: float(D[c].mean()) if c in D else float('nan')
            for c in ('planarity_share', 'ez_gt_consistency', 'ez_fidelity')}
    print(f'SUMMARY stereo_eval molecules {len(D)}; planarity share {mean["planarity_share"]:.4f}; '
          f'GT consistency {mean["ez_gt_consistency"]:.4f}; fidelity {mean["ez_fidelity"]:.4f}')


if __name__ == '__main__':
    main()
