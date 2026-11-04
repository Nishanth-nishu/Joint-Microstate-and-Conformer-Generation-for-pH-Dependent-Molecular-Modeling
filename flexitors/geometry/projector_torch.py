"""Batched torch version of the dihedral-preserving projector P_W (EXP-1 §2.2, §5.2 D-006: no per-molecule Python loop).

For a PyG batch (positions ``pos`` (N_tot, 3), graph index ``batch``, W quadruples ``quads`` (M_tot, 4) in batch atom
numbering) every molecule b gets C_b = [W_b; J_rig,b] with W_b the analytic dihedral gradients of its quadruples and
J_rig,b its 6 rigid-motion rows. All C_b are padded into one tensor (B, K, 3·N_max) with K = max_b m_b + 6, so

    P_b v = v − C_bᵀ (C_b C_bᵀ)⁺ C_b v,   rank P_b = 3 n_b − rank(C_b C_bᵀ)

is computed for the whole batch with one batched pseudo-inverse. Zero padding rows and columns do not change the
projector of the real block. Same formulas as ``flexitors.geometry.projector`` (numpy, checked against it in
``tests/test_projector_torch.py``); float64 throughout because CCᵀ can be ill-conditioned near-linear quadruples.
"""
from __future__ import annotations

import torch


def dihedral_gradients_torch(pos: torch.Tensor, quads: torch.Tensor) -> torch.Tensor:
    """Gradient of each dihedral (i, j, k, l) w.r.t. its four atoms, shape (M, 4, 3) (Blondel & Karplus 1996, eq. 27;
    sign convention of ``flexitors.geometry.internal.dihedral``)."""
    p = pos.double()
    i, j, k, l = quads.long().unbind(1)
    F, G, H = p[i] - p[j], p[j] - p[k], p[l] - p[k]
    A, B = torch.cross(F, G, dim=1), torch.cross(H, G, dim=1)
    a2, b2 = (A * A).sum(1, keepdim=True), (B * B).sum(1, keepdim=True)
    g = G.norm(dim=1, keepdim=True)
    fg, hg = (F * G).sum(1, keepdim=True), (H * G).sum(1, keepdim=True)
    gi = -g / a2 * A
    gl = g / b2 * B
    gj = g / a2 * A + fg / (a2 * g) * A - hg / (b2 * g) * B
    gk = -g / b2 * B - fg / (a2 * g) * A + hg / (b2 * g) * B
    return torch.stack([gi, gj, gk, gl], 1)


class BatchProjector:
    """P_W for every molecule of a batch at fixed positions; ``project`` maps node vectors (N_tot, 3) to node vectors."""

    def __init__(self, pos: torch.Tensor, batch: torch.Tensor, quads: torch.Tensor, rcond: float = 1e-10):
        dev = pos.device
        p = pos.double()
        n_graphs = int(batch.max()) + 1 if len(batch) else 0
        counts = torch.bincount(batch, minlength=n_graphs)
        ptr = torch.cat([counts.new_zeros(1), counts.cumsum(0)[:-1]])
        self.batch, self.local = batch, torch.arange(len(batch), device=dev) - ptr[batch]
        self.n_graphs, self.n_max, self.counts = n_graphs, int(counts.max()) if n_graphs else 0, counts
        qb = batch[quads[:, 1].long()] if len(quads) else batch.new_zeros(0)
        m = torch.bincount(qb, minlength=n_graphs)
        m_max = int(m.max()) if len(quads) else 0
        K, D = m_max + 6, 3 * self.n_max
        C = torch.zeros(n_graphs, K, D, dtype=torch.float64, device=dev)
        if len(quads):
            qptr = torch.cat([m.new_zeros(1), m.cumsum(0)[:-1]])
            order = torch.argsort(qb, stable=True)
            quads, qb = quads[order], qb[order]
            row = torch.arange(len(quads), device=dev) - qptr[qb]
            grad = dihedral_gradients_torch(p, quads)  # (M, 4, 3)
            for t in range(4):
                atom = quads[:, t].long()
                for ax in range(3):
                    C.index_put_((qb, row, 3 * self.local[atom] + ax), grad[:, t, ax], accumulate=True)
        # rigid rows: translations, then rotations e × (x − centroid)
        cent = torch.zeros(n_graphs, 3, dtype=torch.float64, device=dev).index_add_(0, batch, p) / counts[:, None]
        c = p - cent[batch]
        for ax in range(3):
            C[batch, m_max + ax, 3 * self.local + ax] = 1.0
            e = torch.zeros(3, dtype=torch.float64, device=dev)
            e[ax] = 1.0
            rot = torch.cross(e.expand_as(c), c, dim=1)
            for d in range(3):
                C[batch, m_max + 3 + ax, 3 * self.local + d] = rot[:, d]
        M = C @ C.transpose(1, 2)
        Mp = torch.linalg.pinv(M, rcond=rcond, hermitian=True)
        self.C, self.Mp = C, Mp
        self.rank = 3 * counts - torch.diagonal(Mp @ M, dim1=1, dim2=2).sum(1).round().long()

    def _dense(self, v: torch.Tensor) -> torch.Tensor:
        out = torch.zeros(self.n_graphs, self.n_max, 3, dtype=torch.float64, device=v.device)
        out[self.batch, self.local] = v.double()
        return out.reshape(self.n_graphs, -1)

    def project(self, v: torch.Tensor) -> torch.Tensor:
        """P_b v_b for every molecule; v (N_tot, 3) in batch order; returns float64 (N_tot, 3), differentiable in v."""
        d = self._dense(v)
        cv = torch.einsum('bkd,bd->bk', self.C, d)
        d = d - torch.einsum('bkd,bk->bd', self.C, torch.einsum('bkl,bl->bk', self.Mp, cv))
        return d.reshape(self.n_graphs, self.n_max, 3)[self.batch, self.local]
