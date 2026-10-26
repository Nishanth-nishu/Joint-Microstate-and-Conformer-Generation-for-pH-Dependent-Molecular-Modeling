"""flexitors.sampling.stereo_guard (EXP-3 §3.2-§3.3, §6 tests 4-6, §9.2; DECISION fix 3). Positions are moved by a
test-only replica of TD's ``modify_conformer`` rotation (utils/torsion.py:57-75); the same operators are checked against
TD's own code on gnode118 (test_reverse.py)."""
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

from flexitors.geometry.internal import dihedral, wrap_angle
from flexitors.geometry.torsion_types import td_rotatable_bonds, td_reference_quadruple, td_rotated_side
from flexitors.sampling.stereo_guard import (LATE_FOLD_SIGMA, UPDATE_SIGN, StereoHooks, constrainable_bonds,
                                             fold_half_circle)

SMI = 'C/C=C/C(C)C=CC(C)C=N'  # specified C=C, unspecified C=C, imine C=N-H (non-stereo)


def make_case(smiles=SMI, seed=4):
    mh = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mh, randomSeed=seed) == 0
    rot = td_rotatable_bonds(mh)
    sides = [td_rotated_side(mh, k) for k in rot]
    quads = [td_reference_quadruple(mh, k) for k in rot]
    targets = {}
    for k in rot:
        b = mh.GetBondWithIdx(k)
        if b.GetStereo() != Chem.BondStereo.STEREONONE:
            targets[k] = np.pi  # C/C=C/C: trans on TD's methyl-methyl quadruple
    return mh, mh.GetConformer().GetPositions(), rot, sides, quads, targets


def apply(pos, quads, sides, upd):
    """modify_conformer semantics: rotate the v side about (pos[u] - pos[v]) by the update, about pos[v]"""
    p = pos.copy()
    for q, side, u_ in zip(quads, sides, upd):
        if u_ == 0:
            continue
        u, v = q[1], q[2]
        axis = (p[u] - p[v]) / np.linalg.norm(p[u] - p[v])
        rot = Rotation.from_rotvec(axis * u_).as_matrix()
        idx = sorted(side)
        p[idx] = (p[idx] - p[v]) @ rot.T + p[v]
    return p


def test_update_sign_matches_modify_conformer():
    mh, pos, rot, sides, quads, _ = make_case()
    for s in range(len(rot)):
        upd = np.zeros(len(rot))
        upd[s] = 0.1
        d = wrap_angle(dihedral(apply(pos, quads, sides, upd), np.array([quads[s]])) -
                       dihedral(pos, np.array([quads[s]])))[0]
        assert d == pytest.approx(UPDATE_SIGN * 0.1, abs=1e-9)


def test_fold_range_and_involution():
    x = np.random.default_rng(0).uniform(-10, 10, 200000)
    f = fold_half_circle(x)
    assert np.all(np.abs(f) <= np.pi / 2 + 1e-12)
    assert np.allclose(fold_half_circle(f), f)
    assert np.allclose(fold_half_circle(np.array([np.pi / 2, -np.pi / 2, 0.3])), [np.pi / 2, -np.pi / 2, 0.3])


def test_bond_classes():
    mh, _, _, _, _, targets = make_case()
    bonds, m = constrainable_bonds(mh, targets)
    assert sorted(b.kind for b in bonds) == ['nonstereo', 'specified', 'unspecified']
    assert m == len(td_rotatable_bonds(mh))
    assert {b.kind for b in constrainable_bonds(Chem.AddHs(Chem.MolFromSmiles('CC#CC')), {})[0]} == {'triple'}


def walk(arm, steps=20, n_conf=4, seed=0, start_wrong=False):
    """a random torsion walk standing in for the reverse SDE: seed stage, then per step update -> fold -> verify"""
    mh, pos, rot, sides, quads, targets = make_case()
    bonds, m = constrainable_bonds(mh, targets)
    hooks = StereoHooks(arm, bonds, m, steps)
    spec = [b for b in bonds if b.kind == 'specified'][0]
    rng = np.random.default_rng(seed)
    P = np.stack([pos] * n_conf)
    if start_wrong:  # seeds in the wrong isomer
        flip = np.zeros(m)
        flip[spec.slot] = np.pi
        P = np.stack([apply(p, quads, sides, flip) for p in P])
    tau0 = dihedral(P, np.array([spec.quad]))[:, 0]
    seed_upd = hooks.seed_updates(P, rng.uniform(-np.pi, np.pi, (n_conf, m)))
    P = np.stack([apply(p, quads, sides, u) for p, u in zip(P, seed_upd)])
    P = np.stack([apply(p, quads, sides, u) for p, u in zip(P, hooks.fold_updates(list(range(n_conf)), P, None, None))])
    sig = 10 ** np.linspace(np.log10(np.pi), np.log10(0.01 * np.pi), steps + 1)[:-1]
    for k, s in enumerate(sig):
        upd = hooks.mask_update(rng.normal(0, 0.6, (n_conf, m)))
        P = np.stack([apply(p, quads, sides, u) for p, u in zip(P, upd)])
        P = np.stack([apply(p, quads, sides, u) for p, u in zip(P, hooks.fold_updates(list(range(n_conf)), P, s, k))])
        hooks.verify(list(range(n_conf)), P, k)
    return hooks, bonds, P, spec, tau0


def test_k1_freezes_and_seed_checks():
    hooks, bonds, P, spec, tau0 = walk('K1', start_wrong=True)
    tau = dihedral(P, np.array([spec.quad]))[:, 0]
    assert np.allclose(np.cos(tau - np.pi), 1, atol=1e-6)  # flipped to the target once, then never moved
    assert hooks.counts['ez_violations'].sum() == 0
    hooks, _, P, spec, tau0 = walk('K1')
    assert np.allclose(wrap_angle(dihedral(P, np.array([spec.quad]))[:, 0] - tau0), 0, atol=1e-9)


def test_k1_masks_both_seed_randomisation_and_updates():
    mh, pos, rot, sides, quads, targets = make_case()
    bonds, m = constrainable_bonds(mh, targets)
    hooks = StereoHooks('K1', bonds, m, 20)
    spec = [b for b in bonds if b.kind == 'specified'][0]
    u = hooks.seed_updates(np.stack([pos]), np.ones((1, m)))
    assert u[0, spec.slot] == 0 and (np.delete(u[0], spec.slot) == 1).all()
    assert hooks.mask_update(np.ones((2, m)))[:, spec.slot].tolist() == [0, 0]


def test_k1_wrong_freezes_the_opposite_isomer():
    hooks, _, P, spec, _ = walk('K1_wrong')
    assert np.allclose(np.cos(dihedral(P, np.array([spec.quad]))[:, 0]), 1, atol=1e-6)  # cis on the methyls


def test_k1_canonical_freezes_all_double_bonds():
    hooks, bonds, P, _, _ = walk('K1_canonical')
    assert len(hooks.frozen) == len(bonds)


def test_k2_inf_folds_specified_and_late_folds_the_rest():
    hooks, bonds, P, spec, _ = walk('K2_inf', start_wrong=True)
    assert np.all(np.cos(dihedral(P, np.array([spec.quad]))[:, 0] - np.pi) >= -1e-9)  # in I around the target
    late = [b for b in bonds if b.kind in ('unspecified', 'nonstereo')]
    assert hooks.late == late and len(hooks.late_phi) == 4
    tau = dihedral(P, np.array([b.quad for b in late]))
    for i in range(4):
        assert np.all(np.cos(tau[i] - hooks.late_phi[i]) >= -1e-9)
    assert hooks.counts['ez_violations'].sum() == 0 and hooks.counts['fold_corrections'].sum() > 0
    assert LATE_FOLD_SIGMA == pytest.approx(0.1 * np.pi)


def test_o_counts_violations_and_moves_freely():
    hooks, _, _, _, _ = walk('O', start_wrong=True)
    assert hooks.counts['ez_violations'].sum() > 0 and hooks.counts['fold_corrections'].sum() == 0
