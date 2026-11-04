"""τ-frozen geometry read-out of a V2 model (EXP-1 §10.2; the M1 / A6 / A6b sampler).

From a seed conformer x (no torsion randomisation, torsion clock σ = σ_min throughout), N = 20 uniform steps on
t ∈ [0, 1]: one network call gives Δ̂, x̂₁ = x + P_W(x) Δ̂, then
    x ← x + (Δt / (1 − t)) P_W(x)(x̂₁ − x)  for n < N − 1,   x ← x + P_W(x)(x̂₁ − x)  at the last step,
each step followed by exact torsion restoration and the §2.6/§2.8 guard (``flexitors.geometry.l_step.guarded_step``).
Deterministic given the seed conformer and the model (no noise is drawn). One molecule per call (CPU or GPU).
"""
from __future__ import annotations

import numpy as np
import torch
from torch_geometric.data import Batch, Data

from flexitors.geometry.internal import dihedral
from flexitors.geometry.l_step import Guard, guarded_step
from flexitors.geometry.projector import dihedral_preserving_projector


@torch.no_grad()
def l_readout(model, data, seed_pos: np.ndarray, quads: np.ndarray, sides: np.ndarray, guard: Guard,
              sigma_min: float, n_steps: int = 20, device='cpu', project: bool = True):
    """Returns (final positions, counters). ``data``: a featurised molecule (TD Data: x, edge_index, edge_attr,
    edge_mask); ``quads``/``sides``: W quadruples and their moving sides; ``guard`` fixed at the seed.
    ``project=False``: the projector-off ablation, step Δ̂ unprojected (torsion restoration and guard unchanged)."""
    model.eval()
    x = np.array(seed_pos, dtype=np.float64)
    tau0 = dihedral(x, quads) if len(quads) else np.zeros(0)
    tot = dict(halvings=0, skips=0, chirality_rejections=0)
    for n in range(n_steps):
        t = n / n_steps
        d = Data(x=data.x, edge_index=data.edge_index, edge_attr=data.edge_attr,
                 edge_mask=torch.as_tensor(data.edge_mask), pos=torch.as_tensor(x, dtype=torch.float32))
        d.node_sigma = sigma_min * torch.ones(len(x))
        d.node_lambda = t * torch.ones(len(x))
        b = Batch.from_data_list([d]).to(device)
        dhat = model(b).l_pred.double().cpu().numpy()
        if project:
            P, _ = dihedral_preserving_projector(x, quads)
            step = (P @ dhat.reshape(-1)).reshape(x.shape)
        else:
            step = dhat  # P_W(x)(x̂₁ − x) = P_W(x) P_W(x) Δ̂ = P_W(x) Δ̂
        if n < n_steps - 1:
            step = step * (1.0 / n_steps) / (1.0 - t)
        x, info = guarded_step(x, step, quads, sides, tau0, guard)
        tot['halvings'] += info['halvings']; tot['skips'] += int(info['skipped'])
        tot['chirality_rejections'] += info['chirality_rejections']
    return x, tot
