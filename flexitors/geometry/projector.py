"""Dihedral-preserving projector P_W for local-geometry (L) steps (EXP-1 §2.2 revised after D-004; §10.1 item 2).

At a conformer x (n atoms) with one chosen dihedral per rotatable bond (the W rows):
  W(x)     (m, 3n)  analytic dihedral gradients (Wilson B-matrix rows) of the m chosen quadruples;
  J_rig(x) (6, 3n)  three translations and three rotations about the centroid (TD Eq. 10, PDF p. 7);
  C = [W; J_rig],  P_W = I − Cᵀ(CCᵀ)⁺C,
the orthogonal projector onto ker W ∩ rigid^⊥ (rank 3n − 6 − m when C has full rank). A step δx = P_W v changes the
chosen dihedrals by zero to first order and adds no rigid motion. Dense float64 numpy, one molecule at a time: this is
the M0 (CPU, no learning) implementation; a training run needs a batched torch version (EXP-1 §5.2, D-006).
"""
from __future__ import annotations

import numpy as np


def dihedral_gradients(pos: np.ndarray, quads: np.ndarray) -> np.ndarray:
    """Analytic gradient of each dihedral (i, j, k, l) w.r.t. all coordinates, shape (m, 3n).

    Blondel & Karplus, J. Comput. Chem. 17, 1132 (1996), eq. 27, with F = x_i − x_j, G = x_j − x_k, H = x_l − x_k,
    A = F × G, B = H × G; sign convention of ``flexitors.geometry.internal.dihedral`` (RDKit), checked against central
    finite differences in ``tests/test_projector.py``.
    """
    p = np.asarray(pos, dtype=np.float64)
    q = np.asarray(quads, dtype=np.int64).reshape(-1, 4)
    W = np.zeros((len(q), p.size))
    for r, (i, j, k, l) in enumerate(q):
        F, G, H = p[i] - p[j], p[j] - p[k], p[l] - p[k]
        A, B = np.cross(F, G), np.cross(H, G)
        a2, b2, g = A @ A, B @ B, np.linalg.norm(G)
        gi = -g / a2 * A
        gl = g / b2 * B
        gj = g / a2 * A + (F @ G) / (a2 * g) * A - (H @ G) / (b2 * g) * B
        gk = -g / b2 * B - (F @ G) / (a2 * g) * A + (H @ G) / (b2 * g) * B
        for atom, gv in ((i, gi), (j, gj), (k, gk), (l, gl)):
            W[r, 3 * atom:3 * atom + 3] += gv
    return W


def rigid_jacobian(pos: np.ndarray) -> np.ndarray:
    """Rows spanning rigid motions at ``pos``: 3 translations, 3 rotations e × (x − centroid); shape (6, 3n)."""
    p = np.asarray(pos, dtype=np.float64)
    n = len(p)
    J = np.zeros((6, 3 * n))
    c = p - p.mean(0)
    for ax in range(3):
        J[ax, ax::3] = 1.0
        e = np.zeros(3)
        e[ax] = 1.0
        J[3 + ax] = np.cross(e, c).reshape(-1)
    return J


def dihedral_preserving_projector(pos: np.ndarray, quads: np.ndarray, rcond: float = 1e-10) -> tuple[np.ndarray, int]:
    """P_W = I − Cᵀ(CCᵀ)⁺C with C = [W; J_rig] (m = 0 gives the rigid-motion complement); returns (P_W, rank P_W)."""
    p = np.asarray(pos, dtype=np.float64)
    C = rigid_jacobian(p)
    if len(quads):
        C = np.vstack([dihedral_gradients(p, quads), C])
    M = C @ C.T
    P = np.eye(p.size) - C.T @ np.linalg.pinv(M, rcond=rcond, hermitian=True) @ C
    return P, int(round(np.trace(P)))
