"""One TD reverse-SDE loop with hooks (round5/DECISION.md graft 1; plan_P2 §2.5 ``reverse_loop``).

Port of the SDE branch of ``torsional-diffusion/diffusion/sampling.py:116-246`` (no ODE, likelihood, pdb or Particle
Guidance; PG arms call TD's own ``sample()``). With ``noise=None`` and ``hooks=None`` it executes the same statements in
the same order as TD, so it is bit-identical to ``sample()`` on CPU (golden test ``test_reverse.py`` on gnode118).

Two additions, both default-off:
* ``noise``: a callable giving the step noise z per conformer from its own stream (``flexitors.sampling.seeding``),
  so a conformer's trajectory does not depend on its batch (EXP-3 §9.1; EXP-5 §12.1). z is still drawn for every
  rotatable bond before any hook acts, so constrained arms share the noise of all other bonds (common random numbers);
* ``hooks``: a ``StereoHooks``: the update mask after the score step, a fold correction applied with TD's own
  ``perturb_batch``, and the per-step verification (EXP-3 §3.2-§3.3).
NFE: one score-model call per conformer per step (EXP-5 §4.1), counted.
"""
from __future__ import annotations

import copy
from typing import Callable

import numpy as np

from flexitors.td import import_baseline


def reverse_diffusion(conformers: list, model, sigma_max: float, sigma_min: float, steps: int, batch_size: int,
                      device, noise: Callable[[list[int], int], np.ndarray] | None = None, hooks=None) -> int:
    """Run the reverse SDE in place on ``conformers`` (TD ``Data`` objects of one molecule); return the NFE spent.

    Args:
        noise: ``noise(conformer_indices, step) -> array (len(indices), m)``; None = TD's global ``torch.normal``.
        hooks: ``StereoHooks`` of the arm, or None.
    """
    import torch
    sampling = import_baseline('diffusion.sampling')
    conf_dataset = sampling.InferenceDataset(conformers)
    loader = sampling.DataLoader(conf_dataset, batch_size=batch_size, shuffle=False)  # as sampling.py:122
    sigma_schedule = 10 ** np.linspace(np.log10(sigma_max), np.log10(sigma_min), steps + 1)[:-1]  # sampling.py:125
    eps = 1 / steps
    nfe = 0
    for data in loader:
        data_gpu = copy.deepcopy(data).to(device)
        idx = [int(i) for i in data.idx]
        n_conf = len(idx)
        for sigma_idx, sigma in enumerate(sigma_schedule):
            data_gpu.node_sigma = sigma * torch.ones(data.num_nodes, device=device)
            with torch.no_grad():
                data_gpu = model(data_gpu)
            nfe += n_conf
            g = sigma * torch.sqrt(torch.tensor(2 * np.log(sigma_max / sigma_min)))
            if noise is None:
                z = torch.normal(mean=0, std=1, size=data_gpu.edge_pred.shape)
            else:
                z = torch.from_numpy(noise(idx, sigma_idx).reshape(data_gpu.edge_pred.shape)).float()
            score = data_gpu.edge_pred.cpu()
            perturb = g ** 2 * eps * score + g * np.sqrt(eps) * z  # sampling.py:204
            if hooks is not None:
                perturb = torch.from_numpy(hooks.mask_update(perturb.numpy().reshape(n_conf, -1)).reshape(-1)).float()
            conf_dataset.apply_torsion_and_update_pos(data, perturb.numpy())
            if hooks is not None:
                pos = data.pos.numpy().reshape(n_conf, -1, 3).astype(np.float64)
                corr = hooks.fold_updates(idx, pos, float(sigma), sigma_idx)
                if np.any(corr != 0):
                    conf_dataset.apply_torsion_and_update_pos(data, corr.reshape(-1))
                hooks.verify(idx, data.pos.numpy().reshape(n_conf, -1, 3).astype(np.float64), sigma_idx)
            data_gpu.pos = data.pos.to(device)
    return nfe
