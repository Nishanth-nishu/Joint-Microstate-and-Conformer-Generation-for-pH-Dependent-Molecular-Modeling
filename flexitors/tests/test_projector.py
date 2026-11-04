"""flexitors.geometry.projector / l_step: EXP-1 §5.3 test_projector, test_invariance and test_constraints items that
M0 needs (no model, no TD environment)."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from flexitors.geometry.internal import dihedral, signed_volume, wrap_angle
from flexitors.geometry.l_step import Guard, guarded_step, projected_oracle_flow, restore_dihedrals, w_quadruples
from flexitors.geometry.projector import dihedral_gradients, dihedral_preserving_projector, rigid_jacobian
from flexitors.geometry.torsion_types import td_rotated_side, td_rotatable_bonds


def _build(smiles, seed=11):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(m, randomSeed=seed) == 0
    pos = m.GetConformer().GetPositions()
    z = np.array([a.GetAtomicNum() for a in m.GetAtoms()])
    bonds = np.array([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()])
    nbrs = [[n.GetIdx() for n in a.GetNeighbors()] for a in m.GetAtoms()]
    edges, sides = [], []
    for k in td_rotatable_bonds(m):
        side = td_rotated_side(m, k)
        b = m.GetBondWithIdx(k)
        s, e = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        edges.append((e, s) if s in side else (s, e))  # second atom on the rotated side, as TD's masked edge
        row = np.zeros(len(z), dtype=bool); row[list(side)] = True
        sides.append(row)
    quads, kept, _, _ = w_quadruples(z, np.array(edges).reshape(-1, 2), nbrs, pos)
    return m, pos, z, bonds, quads, np.array(sides).reshape(-1, len(z))[kept]


def test_gradients_match_finite_differences():
    _, pos, _, _, quads, _ = _build('CC(O)CC1CCNC1C(=O)N')
    W = dihedral_gradients(pos, quads)
    h, fd = 1e-6, np.zeros_like(W)
    for j in range(pos.size):
        p, q = pos.copy().reshape(-1), pos.copy().reshape(-1)
        p[j] += h; q[j] -= h
        fd[:, j] = wrap_angle(dihedral(p.reshape(-1, 3), quads) - dihedral(q.reshape(-1, 3), quads)) / (2 * h)
    assert np.allclose(W, fd, atol=1e-6)


def test_projector_properties():
    _, pos, _, _, quads, _ = _build('CC(O)CC1CCNC1C(=O)N')
    P, rank = dihedral_preserving_projector(pos, quads)
    assert np.allclose(P, P.T, atol=1e-10) and np.allclose(P @ P, P, atol=1e-10)
    assert np.abs(dihedral_gradients(pos, quads) @ P).max() < 1e-6
    assert np.abs(rigid_jacobian(pos) @ P).max() < 1e-6
    assert rank == pos.size - 6 - len(quads)


def test_m0_molecule_rigid_complement():
    _, pos, _, _, quads, _ = _build('c1ccccc1')
    assert len(quads) == 0
    P, rank = dihedral_preserving_projector(pos, quads)
    assert rank == pos.size - 6 and np.abs(rigid_jacobian(pos) @ P).max() < 1e-8


def test_drift_is_second_order_and_restoration_is_exact():
    _, pos, _, _, quads, sides = _build('CCCC(O)CC(=O)NC')
    P, _ = dihedral_preserving_projector(pos, quads)
    rng = np.random.default_rng(0)
    v = P @ rng.normal(size=pos.size)
    v /= np.sqrt((v ** 2).reshape(-1, 3).sum(1).mean())  # unit RMS per atom
    tau0 = dihedral(pos, quads)
    drift = [np.abs(wrap_angle(dihedral(pos + s * v.reshape(-1, 3), quads) - tau0)).max() for s in (0.01, 0.02, 0.05)]
    assert drift[2] < 1e-2 and 3.0 < drift[1] / drift[0] < 5.0, drift  # ~ step^2
    x = restore_dihedrals(pos + 0.05 * v.reshape(-1, 3), quads, sides, tau0)
    assert np.abs(wrap_angle(dihedral(x, quads) - tau0)).max() < 1e-9


def test_restoration_keeps_bonds_and_angles():
    m, pos, _, bonds, quads, sides = _build('CCCC(O)CC(=O)NC')
    x = restore_dihedrals(pos, quads, sides, dihedral(pos, quads) + 0.7)
    L0 = np.linalg.norm(pos[bonds[:, 0]] - pos[bonds[:, 1]], axis=1)
    L1 = np.linalg.norm(x[bonds[:, 0]] - x[bonds[:, 1]], axis=1)
    assert np.abs(L0 - L1).max() < 1e-9
    assert np.abs(wrap_angle(dihedral(x, quads) - dihedral(pos, quads) - 0.7)).max() < 1e-9


def test_guard_flags_enantiomer_and_backtracks_stretch():
    m, pos, _, bonds, quads, sides = _build('C[C@H](N)C(=O)O')
    g = Guard(m, bonds, pos)
    assert len(g.centres) == 1
    assert g.check(pos) == (True, True)
    assert g.check(-pos)[1] is False  # inversion flips the signed volume
    delta = np.zeros_like(pos); delta[bonds[0, 1]] = (pos[bonds[0, 1]] - pos[bonds[0, 0]]) * 0.5  # +50 % bond stretch
    x, info = guarded_step(pos, delta, quads, sides, dihedral(pos, quads), g)
    assert info['halvings'] >= 2 and not info['skipped']
    L = np.linalg.norm(x[bonds[0, 0]] - x[bonds[0, 1]]) / np.linalg.norm(pos[bonds[0, 0]] - pos[bonds[0, 1]])
    assert abs(L - 1) <= 0.15


def test_oracle_flow_reaches_a_pure_L_target():
    m, pos, z, bonds, quads, sides = _build('CC1CCC(O)CC1C')
    m2 = Chem.Mol(m)
    assert AllChem.MMFFOptimizeMolecule(m2, maxIters=500) in (0, 1)
    tgt = restore_dihedrals(m2.GetConformer().GetPositions(), quads, sides, dihedral(pos, quads))
    a, b = tgt - tgt.mean(0), pos - pos.mean(0)  # Kabsch: P_W removes rigid motion, so the target must be aligned
    U, _, Vt = np.linalg.svd(a.T @ b)
    tgt = a @ (U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt) + pos.mean(0)
    g = Guard(m, bonds, pos)
    x, info = projected_oracle_flow(pos, tgt, quads, sides, g)
    heavy = z > 1
    r0 = np.sqrt(((pos - tgt)[heavy] ** 2).sum(1).mean())
    r1 = np.sqrt(((x - tgt)[heavy] ** 2).sum(1).mean())
    assert r1 < 0.5 * r0 and info['final_w_drift'] < 1e-9 and not info['flipped_centre_at_end'], (r0, r1, info)
    assert np.all(np.sign(signed_volume(x, g.centres)) == g.sign0)
