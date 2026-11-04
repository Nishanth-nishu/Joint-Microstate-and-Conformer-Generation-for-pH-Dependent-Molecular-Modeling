"""Guarded local-geometry (L) step for EXP-1 V2: W-quadruple rule, exact torsion restoration, constraint layer.

* W rule, heavy-first (EXP-1 §10.1 item 2): for each TD-rotatable directed edge (b, c) take the lowest-index heavy
  neighbour of b (≠ c) and of c (≠ b), falling back to H only if a side has no heavy neighbour; a neighbour whose bond
  angle with (b, c) has sin < sin 10° is skipped for the next one; if a side has no acceptable neighbour the row is
  dropped and counted.
* Exact torsion restoration (EXP-1 §10.2): after a geometry step the W dihedrals are measured with atan2 and each
  moving side (TD's ``mask_rotate`` row, the side of c) is rotated back about the bond. Side rotations change no bond
  length or angle and cannot flip a tetrahedral centre (E-JOINT-026).
* Constraint layer (EXP-1 §2.6 / §2.8): every bonded pair within ±15 % of its seed length; no non-bonded heavy-atom
  pair closer than 0.75 × the sum of covalent radii; every tetrahedral stereocentre keeps the sign of its signed
  volume at the seed. A failing step is halved up to 3 times, then skipped; backtracks, skips and chirality
  rejections are counted.
"""
from __future__ import annotations

import math

import numpy as np
from rdkit import Chem

from flexitors.geometry.internal import dihedral, signed_volume, wrap_angle
from flexitors.geometry.projector import dihedral_preserving_projector

SIN_MIN = math.sin(math.radians(10.0))
BOND_WINDOW, CLASH_FACTOR, MAX_HALVINGS = 0.15, 0.75, 3
PT = Chem.GetPeriodicTable()


def _angle_sin(pos, a, b, c):
    u, v = pos[a] - pos[b], pos[c] - pos[b]
    return np.linalg.norm(np.cross(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))


def w_quadruples(z: np.ndarray, rot_edges: np.ndarray, neighbours: list[list[int]], pos: np.ndarray):
    """Heavy-first W quadruples (a, b, c, d) for the directed rotatable edges (b, c); returns (quads, kept row index
    into rot_edges, number of dropped rows, number of rows using an H end)."""
    quads, kept, n_h = [], [], 0
    for r, (b, c) in enumerate(np.asarray(rot_edges)):
        ends = []
        for x, y in ((b, c), (c, b)):
            cand = sorted((v for v in neighbours[x] if v != y), key=lambda v: (z[v] == 1, v))
            ok = [v for v in cand if _angle_sin(pos, v, x, y) >= SIN_MIN]
            ends.append(ok[0] if ok else None)
        if None in ends:
            continue
        n_h += int(z[ends[0]] == 1 or z[ends[1]] == 1)
        quads.append((ends[0], b, c, ends[1]))
        kept.append(r)
    return np.array(quads, dtype=np.int64).reshape(-1, 4), np.array(kept, dtype=np.int64), len(rot_edges) - len(kept), n_h


def kabsch_to(mobile: np.ndarray, ref: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """``mobile`` rigidly moved onto ``ref`` by the Kabsch fit of the atoms ``idx`` (all atoms are moved)."""
    a, b = mobile[idx], ref[idx]
    ca, cb = a.mean(0), b.mean(0)
    U, _, Vt = np.linalg.svd((a - ca).T @ (b - cb))
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1.0, 1.0, d]) @ Vt
    return (mobile - ca) @ R + cb


def _rotate(pos, side, origin, axis, angle):
    k = axis / np.linalg.norm(axis)
    v = pos[side] - origin
    c, s = math.cos(angle), math.sin(angle)
    pos[side] = origin + v * c + np.cross(k, v) * s + np.outer(v @ k, k) * (1 - c)


def restore_dihedrals(pos: np.ndarray, quads: np.ndarray, sides: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rotate the moving side (boolean row of ``sides``, which contains c and d of quad (a, b, c, d)) of each W bond so
    that its dihedral equals ``target`` (right-handed rotation of d about b→c raises the dihedral). Returns a copy."""
    p = np.array(pos, dtype=np.float64)
    for (a, b, c, d), side, t in zip(quads, sides, target):
        delta = wrap_angle(t - dihedral(p, np.array([[a, b, c, d]]))[0])
        if delta != 0.0:
            _rotate(p, side, p[c], p[c] - p[b], delta)
    return p


def stereo_centres(mol: Chem.Mol) -> np.ndarray:
    """Rows (centre, n1, n2, n3) for every potential tetrahedral stereocentre (RDKit FindPotentialStereo, which also
    flags ring-junction centres), neighbours in index order; the sign of the signed volume is what is guarded."""
    rows = []
    for el in Chem.FindPotentialStereo(Chem.Mol(mol)):
        if el.type == Chem.StereoType.Atom_Tetrahedral:
            c = el.centeredOn
            nb = sorted(n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors())
            if len(nb) >= 3:
                rows.append((c, nb[0], nb[1], nb[2]))
    return np.array(rows, dtype=np.int64).reshape(-1, 4)


class Guard:
    """Constraint layer of one molecule, fixed at its seed conformer."""

    def __init__(self, mol: Chem.Mol, bonds: np.ndarray, seed_pos: np.ndarray):
        self.bonds = np.asarray(bonds, dtype=np.int64)
        self.seed_len = np.linalg.norm(seed_pos[self.bonds[:, 0]] - seed_pos[self.bonds[:, 1]], axis=1)
        z = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
        heavy = np.where(z > 1)[0]
        bonded = {tuple(sorted(b)) for b in self.bonds.tolist()}
        pairs = [(i, j) for ii, i in enumerate(heavy) for j in heavy[ii + 1:] if (i, j) not in bonded]
        self.pairs = np.array(pairs, dtype=np.int64).reshape(-1, 2)
        rc = np.array([PT.GetRcovalent(int(x)) for x in z])
        self.min_d = CLASH_FACTOR * (rc[self.pairs[:, 0]] + rc[self.pairs[:, 1]]) if len(self.pairs) else np.zeros(0)
        self.centres = stereo_centres(mol)
        self.sign0 = np.sign(signed_volume(seed_pos, self.centres)) if len(self.centres) else np.zeros(0)

    def chirality_ok(self, pos):
        return not len(self.centres) or bool(np.all(np.sign(signed_volume(pos, self.centres)) == self.sign0))

    def check(self, pos) -> tuple[bool, bool]:
        """(geometry ok: bond window and clash, chirality ok)."""
        L = np.linalg.norm(pos[self.bonds[:, 0]] - pos[self.bonds[:, 1]], axis=1)
        ok = bool(np.all(np.abs(L / self.seed_len - 1.0) <= BOND_WINDOW))
        if ok and len(self.pairs):
            ok = bool(np.all(np.linalg.norm(pos[self.pairs[:, 0]] - pos[self.pairs[:, 1]], axis=1) >= self.min_d))
        return ok, self.chirality_ok(pos)


def guarded_step(pos, delta, quads, sides, tau0, guard: Guard, restore=True):
    """One L step x + s·delta (s = 1, 1/2, 1/4, 1/8) followed by torsion restoration; the first s that passes the
    guard is taken. Returns (new pos, dict(halvings, skipped, chirality_rejections, drift_before_restore))."""
    info = dict(halvings=0, skipped=False, chirality_rejections=0, drift_before_restore=0.0)
    s = 1.0
    for h in range(MAX_HALVINGS + 1):
        x = pos + s * delta
        if len(quads):
            info['drift_before_restore'] = max(info['drift_before_restore'],
                                               float(np.abs(wrap_angle(dihedral(x, quads) - tau0)).max()))
            if restore:
                x = restore_dihedrals(x, quads, sides, tau0)
        geo, chi = guard.check(x)
        if geo and chi:
            info['halvings'] = h
            return x, info
        info['chirality_rejections'] += int(not chi)
        s *= 0.5
    info['halvings'] = MAX_HALVINGS
    info['skipped'] = True
    return pos, info


def projected_oracle_flow(pos0, target, quads, sides, guard: Guard, n_steps=20, restore=True):
    """M0 integration x ← x + (1/(N − n))·P_W(x)(target − x), n = 0 … N − 1, each step guarded and torsion-restored
    (EXP-1 §10.3). Returns (final pos, per-run counters)."""
    x = np.array(pos0, dtype=np.float64)
    tau0 = dihedral(x, quads) if len(quads) else np.zeros(0)
    tot = dict(halvings=0, skips=0, chirality_rejections=0, max_drift_before_restore=0.0)
    for n in range(n_steps):
        P, _ = dihedral_preserving_projector(x, quads)
        delta = (P @ (target - x).reshape(-1)).reshape(x.shape) / (n_steps - n)
        x, info = guarded_step(x, delta, quads, sides, tau0, guard, restore=restore)
        tot['halvings'] += info['halvings']
        tot['skips'] += int(info['skipped'])
        tot['chirality_rejections'] += info['chirality_rejections']
        tot['max_drift_before_restore'] = max(tot['max_drift_before_restore'], info['drift_before_restore'])
    tot['final_w_drift'] = float(np.abs(wrap_angle(dihedral(x, quads) - tau0)).max()) if len(quads) else 0.0
    tot['flipped_centre_at_end'] = not guard.chirality_ok(x)
    return x, tot
