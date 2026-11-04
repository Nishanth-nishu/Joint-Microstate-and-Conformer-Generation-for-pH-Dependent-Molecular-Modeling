"""EXP-1 losses: TD's torsion DSM loss (verbatim formula of ``utils/training.py:25-31``) and the §10.1 item 4
projected endpoint loss of the V2 geometry head."""
from __future__ import annotations

import torch

from flexitors.geometry.projector_torch import BatchProjector
from flexitors.td import import_baseline


def td_loss(data) -> tuple[torch.Tensor, torch.Tensor]:
    """(loss, base loss) exactly as TD's ``train_epoch``: mean over rotatable bonds of (score − pred)² / score_norm."""
    torus = import_baseline('diffusion.torus')
    pred = data.edge_pred
    if pred.numel() == 0:  # every molecule in the batch is rigid (see FlexiTorsScoreModel.forward)
        z = pred.sum() * 0.0
        return z, z.detach()
    score = torch.tensor(torus.score(data.edge_rotate.cpu().numpy(), data.edge_sigma.cpu().numpy()), device=pred.device)
    norm = torch.tensor(torus.score_norm(data.edge_sigma.cpu().numpy()), device=pred.device)
    return ((score - pred) ** 2 / norm).mean(), (score ** 2 / norm).mean()


def projected_l_loss(data, project: bool = True) -> tuple[torch.Tensor, torch.Tensor, int]:
    """(loss, base loss, number of molecules with an active geometry loss).

    x̂₁ = x_t + P_W(x_t) Δ̂; per molecule ‖P_W(x_t)(x̂₁ − x_G*′)‖² / rank P_W, averaged over the molecules whose
    ``geo_mask`` is 1 (pair_ok and no unresolved chirality flip). The base loss is the same with Δ̂ = 0.
    ``project=False`` (EXP1-V2-T1 ablation): P_W = I, rank = 3N; x̂₁ = x_t + Δ̂, loss ‖x̂₁ − x_G*′‖² / 3N.
    """
    resid = (data.pos - data.geo_target).double()
    n_graphs = int(data.batch.max()) + 1
    if project:
        proj = BatchProjector(data.pos.detach(), data.batch, data.w_quad_index.T)
        r = proj.project(resid + proj.project(data.l_pred))
        r0 = proj.project(resid)
        rank = proj.rank.clamp(min=1).double()
    else:
        r, r0 = resid + data.l_pred.double(), resid
        rank = 3.0 * torch.zeros(n_graphs, dtype=torch.float64, device=r.device).index_add_(
            0, data.batch, torch.ones(len(r), dtype=torch.float64, device=r.device))
    per = torch.zeros(n_graphs, dtype=torch.float64, device=r.device).index_add_(0, data.batch, (r ** 2).sum(1)) / rank
    per0 = torch.zeros(n_graphs, dtype=torch.float64, device=r.device).index_add_(0, data.batch, (r0 ** 2).sum(1)) / rank
    mask = data.geo_mask.double()
    n = int(mask.sum().item())
    denom = max(n, 1)
    return (per * mask).sum() / denom, (per0 * mask).sum() / denom, n
