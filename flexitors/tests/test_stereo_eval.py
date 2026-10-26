"""flexitors.eval.stereo_eval and the F-flip arm (EXP-3 §3.4, §4, §9.3; R-CHEM conformance items 1-2)."""
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from flexitors.eval.stereo_eval import double_bonds, gt_states, molecule_endpoints
from flexitors.geometry.internal import dihedral
from flexitors.sampling.stereo_fix import f_flip

SMI = 'C/C=C/C(C)C=CCC=N'   # specified C=C, unspecified C=C, imine C=N-H


def confs_of(smiles, n, seed=1):
    mh = Chem.AddHs(Chem.MolFromSmiles(smiles))
    cids = list(AllChem.EmbedMultipleConfs(mh, n, randomSeed=seed))
    AllChem.MMFFOptimizeMoleculeConfs(mh)
    return mh, [Chem.Mol(mh, confId=c) for c in cids]


def set_db(m, quad, value):
    m = Chem.Mol(m)
    rdMolTransforms.SetDihedralRad(m.GetConformer(), *map(int, quad), value)
    return m


def test_double_bond_inventory():
    mh = Chem.AddHs(Chem.MolFromSmiles(SMI))
    dbs = double_bonds(mh)
    assert len(dbs) == 3 and sum(s for _, _, s in dbs) == 2  # imine N-H is non-stereo


def test_gt_consistency_survives_atom_reordering():
    mh, gen = confs_of(SMI, 4)
    order = list(np.random.default_rng(0).permutation(mh.GetNumAtoms()))
    gts = [Chem.RenumberAtoms(g, [int(i) for i in order]) for g in gen[:2]]  # GT in another atom order
    struct = np.array([q for _, q, s in double_bonds(mh) if s])
    allowed = gt_states(mh, gts, struct)
    own = {tuple(bool(x) for x in np.cos(dihedral(g.GetConformer().GetPositions(), struct)) > 0) for g in gen[:2]}
    assert own <= allowed
    r = molecule_endpoints(mh, gen[:2], gts, {})
    assert r['ez_gt_consistency'] == 1.0
    k, q, _ = [d for d in double_bonds(mh) if d[2]][1]
    tau = dihedral(gen[0].GetConformer().GetPositions(), np.array([q]))[0]
    wrong = set_db(gen[0], q, tau + math.pi)
    if tuple(bool(x) for x in np.cos(dihedral(wrong.GetConformer().GetPositions(), struct)) > 0) not in allowed:
        assert molecule_endpoints(mh, [wrong], gts, {})['ez_gt_consistency'] == 0.0


def test_fidelity_planarity_and_f_flip():
    mh, gen = confs_of(SMI, 3)
    (k, q, _) = [d for d in double_bonds(mh) if d[2]][0]
    targets = {k: math.pi}
    r = molecule_endpoints(mh, gen, [], targets)
    assert r['ez_fidelity'] == 1.0 and r['planarity_share'] == 1.0 and r['n_specified'] == 1
    bad = [set_db(gen[0], q, 0.0), set_db(gen[1], q, math.pi / 2)]
    r = molecule_endpoints(mh, bad, [], targets)
    assert r['ez_fidelity'] == 0.0 and r['planarity_share'] < 1.0
    fixed, n = f_flip(bad, targets)
    assert n == 2
    for m in fixed:
        assert math.cos(dihedral(m.GetConformer().GetPositions(), np.array([q]))[0] - math.pi) > 0.999
    b0 = [rdMolTransforms.GetBondLength(bad[0].GetConformer(), b.GetBeginAtomIdx(), b.GetEndAtomIdx())
          for b in mh.GetBonds()]
    b1 = [rdMolTransforms.GetBondLength(fixed[0].GetConformer(), b.GetBeginAtomIdx(), b.GetEndAtomIdx())
          for b in mh.GetBonds()]
    assert np.allclose(b0, b1, atol=1e-6)
    assert molecule_endpoints(mh, fixed, [], targets)['ez_fidelity'] == 1.0
