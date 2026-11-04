"""EXP-1 §10.1 V2 training transform (``PairedJointTransform``) on the round-2 paired cache (``cache_paired.*``).

One sample of a molecule (datapoint with ``pos``/``gt_pos``/``pair_ok`` lists, ``edge_mask``, ``mask_rotate``):
1. conformer k uniformly (TD's ``random.choice``, ``utils/dataset.py:48-55``); x_R = pos[k], x_G = gt_pos[k];
2. W quadruples (heavy-first rule, §10.1 item 2) on x_R; pure-L target x_G* = x_G with its W dihedrals set to x_R's by
   rigid side rotations, Kabsch-aligned to x_R on heavy atoms (§10.1 item 3);
3. torsion clock σ ~ logU[σ_min, σ_max], Δτ ~ N(0, σ²) per TD-rotatable bond, applied to BOTH endpoints with TD's
   ``modify_conformer`` (§2.2 step 3): x_R′, x_G*′;
4. geometry clock t ~ U[0, 1]; x_t = interp(x_R′, x_G*′, t) + γ(t)·P_W(x_R′) z, γ(t) = σ_b √(t(1−t)), z ~ N(0, I);
   ``interp`` is TD's ``interp_x_torch`` (terminal bond lengths kept) (§10.1 item 4);
5. chirality (§10.1 item 5): if x_t flips a potential tetrahedral centre relative to x_R′, z is redrawn (≤ 5 tries);
   if it still flips, the geometry loss of the sample is masked and counted;
6. ``pair_ok`` False (§10.1 item 1): kept with the torsion loss only (x_t = x_R′, geometry loss masked).
Outputs on the Data object: ``pos`` = x_t, ``geo_target`` = x_G*′, ``node_lambda`` = t (the λ-embedding path,
§10.1 item 6), ``node_sigma``, ``edge_rotate`` = Δτ (TD's torsion score target), ``w_quad_index`` (shape (4, m); PyG
increments it by the node offset and concatenates on the last dimension when batching), ``geo_mask`` and the counters ``n_chirality_redraws``, ``chirality_masked``,
``pair_masked``. RNG: numpy's global generator, as TD's transform; DataLoader workers seed it with
``seed_numpy_worker``.
"""
from __future__ import annotations

import math
import random

import numpy as np
import torch

from flexitors.geometry.internal import dihedral
from flexitors.geometry.l_step import Guard, kabsch_to, restore_dihedrals, w_quadruples
from flexitors.geometry.projector import dihedral_preserving_projector
from flexitors.td import import_baseline

MAX_Z_TRIES = 5


def seed_numpy_worker(worker_id: int) -> None:
    """DataLoader ``worker_init_fn``: seed numpy (and ``random``) from the per-worker torch seed, which torch derives
    from the main-process RNG every epoch; without it all workers share numpy's state."""
    torch.set_num_threads(1)  # workers do numpy/RDKit work; inherited OMP threads oversubscribe the CPUs
    s = torch.initial_seed() % 2 ** 32
    np.random.seed(s)
    random.seed(s)


class PairedJointTransform:
    """§10.1 V2 transform; ``sigma_b`` in Å (γ amplitude), σ range as TD's training (log-uniform)."""

    def __init__(self, sigma_b: float, sigma_min: float = 0.01 * math.pi, sigma_max: float = math.pi):
        self.sigma_b, self.sigma_min, self.sigma_max = sigma_b, sigma_min, sigma_max
        self._torsion = import_baseline('utils.torsion')
        self._dataset = import_baseline('utils.dataset')

    def _modify(self, pos, rot_edges, mask_rotate, dtau):
        return self._torsion.modify_conformer(pos.copy(), torch.as_tensor(rot_edges), mask_rotate, dtau, as_numpy=True)

    def __call__(self, data):
        n = len(data.pos)
        k = random.randrange(n)
        xR = data.pos[k].numpy().astype(np.float64)
        xG = data.gt_pos[k].numpy().astype(np.float64)
        pair_ok = bool(data.pair_ok[k])
        mol = data.mol
        z = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
        ei = data.edge_index.numpy()
        em = np.asarray(data.edge_mask, dtype=bool)
        rot_edges = ei.T[em]
        mask_rotate = np.asarray(data.mask_rotate, dtype=bool)
        nbrs = [[] for _ in range(len(z))]
        for u, v in ei[:, ::2].T:
            nbrs[u].append(int(v)); nbrs[v].append(int(u))
        quads, kept, _, _ = w_quadruples(z, rot_edges, nbrs, xR)
        heavy = np.where(z > 1)[0]

        sigma = float(np.exp(np.random.uniform(np.log(self.sigma_min), np.log(self.sigma_max))))
        dtau = np.random.normal(0.0, sigma, size=int(em.sum()))
        xRp = self._modify(xR, rot_edges, mask_rotate, dtau)
        t = float(np.random.uniform())
        redraws, chi_masked = 0, False
        if pair_ok:
            xGs = restore_dihedrals(xG, quads, mask_rotate[kept], dihedral(xR, quads)) if len(quads) else xG.copy()
            xGs = kabsch_to(xGs, xR, heavy)
            xGsp = self._modify(xGs, rot_edges, mask_rotate, dtau)
            base = self._dataset.interp_x_torch(torch.as_tensor(xRp), torch.as_tensor(xGsp), t, data.edge_index).numpy()
            P, _ = dihedral_preserving_projector(xRp, quads)
            guard = Guard(mol, ei[:, ::2].T, xRp)
            gamma = self.sigma_b * math.sqrt(t * (1.0 - t))
            for _ in range(MAX_Z_TRIES):
                xt = base + gamma * (P @ np.random.normal(size=xRp.size)).reshape(xRp.shape)
                if guard.chirality_ok(xt):
                    break
                redraws += 1
            else:
                chi_masked = True
            target = xGsp
        else:
            xt, target = xRp, xRp
        out = data
        out.pos = torch.as_tensor(xt, dtype=torch.float32)
        out.geo_target = torch.as_tensor(target, dtype=torch.float32)
        out.node_sigma = sigma * torch.ones(len(z))
        out.node_lambda = t * torch.ones(len(z))
        out.edge_rotate = torch.as_tensor(dtau)
        out.w_quad_index = torch.as_tensor(quads, dtype=torch.long).reshape(-1, 4).T.contiguous()  # (4, m): PyG cats *index* keys on dim -1
        out.geo_mask = torch.tensor([float(pair_ok and not chi_masked)])
        out.n_chirality_redraws = torch.tensor([redraws])
        out.chirality_masked = torch.tensor([float(chi_masked)])
        out.pair_masked = torch.tensor([float(not pair_ok)])
        for key in ('gt_pos', 'pair_ok', 'weights', 'mol'):
            if key in out:
                del out[key]
        return out
