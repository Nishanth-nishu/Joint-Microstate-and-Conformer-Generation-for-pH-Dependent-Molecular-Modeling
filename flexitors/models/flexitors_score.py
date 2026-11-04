"""FlexiTors V2 score model: TD's ``TensorProductScoreModel`` plus a per-atom ``1o`` local-geometry head (EXP-1 §2.4,
§5.1 ``flexitors_score.py`` / ``heads.py``, §10.1 items 4 and 6).

* Torsion head: TD's, unchanged (``diffusion/score_model.py:119-132``, ``0o`` pseudoscalar per rotatable bond).
* Geometry clock t: through TD's λ-embedding path (``lambda_embed_dim`` = 32, ``score_model.py:192-199``, the round-2 S4
  mechanism); the transform writes t to ``data.node_lambda``.
* Geometry head: an e3nn ``o3.Linear`` from the last atom layer's irreps to ``1x1o``, so it reads only the ``1o``
  channels (EXP-1 §2.4 table: a Cartesian L velocity/score transforms as ``1o``, v(−x) = −v(x)). It returns the
  per-atom displacement Δ̂; x̂₁ = x_t + P_W(x_t) Δ̂ is formed in the loss/sampler (§10.1 item 4).
TD's ``forward`` is not copied: the last conv layer's output is captured with a forward hook, so with the geometry
head ignored this model computes exactly TD's ``edge_pred`` (golden test ``tests/test_v2.py``).
"""
from __future__ import annotations

from argparse import ArgumentParser, Namespace

import torch
from e3nn import o3

from flexitors.td import import_baseline


def _td_score_model():
    return import_baseline('diffusion.score_model').TensorProductScoreModel


def add_model_args(ap: ArgumentParser) -> None:
    """TD qm9 hyper-parameters (``workdir/pretrained/qm9_default/model_parameters.yml``) + the geometry clock."""
    ap.add_argument('--in_node_features', type=int, default=44)
    ap.add_argument('--in_edge_features', type=int, default=4)
    ap.add_argument('--sigma_embed_dim', type=int, default=32)
    ap.add_argument('--lambda_embed_dim', type=int, default=32)
    ap.add_argument('--radius_embed_dim', type=int, default=50)
    ap.add_argument('--num_conv_layers', type=int, default=4)
    ap.add_argument('--max_radius', type=float, default=5.0)
    ap.add_argument('--ns', type=int, default=32)
    ap.add_argument('--nv', type=int, default=8)
    ap.add_argument('--sigma_min', type=float, default=0.01 * 3.14)
    ap.add_argument('--sigma_max', type=float, default=3.14)
    ap.add_argument('--scale_by_sigma', type=int, default=1)
    ap.add_argument('--use_second_order_repr', type=int, default=1)
    ap.add_argument('--no_residual', type=int, default=0)
    ap.add_argument('--no_batch_norm', type=int, default=0)
    ap.add_argument('--no_parity', type=int, default=0)
    ap.add_argument('--l_head_init', choices=('default', 'zero'), default='default',
                    help='zero: geometry head weights start at 0, so Δ̂ = 0 (the base loss) at initialisation')


def build_model(args: Namespace):
    """The V2 model from TD-style hyper-parameters (``utils/utils.py:get_model`` arguments) with
    ``args.lambda_embed_dim`` > 0 for the geometry clock."""
    base = _td_score_model()

    class FlexiTorsScoreModel(base):
        def __init__(self, l_head_init='default', **kw):
            super().__init__(**kw)
            self.l_head = o3.Linear(self.conv_layers[-1].out_irreps, o3.Irreps('1x1o'))
            if l_head_init == 'zero':
                with torch.no_grad():
                    self.l_head.weight.zero_()
            self._last_atom_layer = None
            self.conv_layers[-1].register_forward_hook(self._keep_last)

        def _keep_last(self, module, inputs, output):
            self._last_atom_layer = output

        def forward(self, data):
            mask = data.edge_mask
            rigid = not bool(mask.any())
            # rigid: TD's bond conv does `out_nodes = out_nodes or n_atoms`, so edge_mask.sum() == 0 silently becomes
            # one output row per atom and the σ-scaling then fails. The atom layers never read edge_mask, so marking
            # one bond only to run TD's forward leaves Δ̂ exact; the torsion output is the empty tensor it should be.
            if rigid:
                assert mask.numel() > 0, 'molecule without bonds'
                data.edge_mask = torch.zeros_like(mask)
                data.edge_mask[0] = True
            data = super().forward(data)
            if rigid:
                data.edge_mask, data.edge_pred = mask, data.edge_pred.reshape(-1)[:0]
            data.l_pred = self.l_head(self._last_atom_layer)
            self._last_atom_layer = None
            return data

    assert getattr(args, 'lambda_embed_dim', 0) > 0, 'V2 needs the geometry clock: lambda_embed_dim > 0'
    return FlexiTorsScoreModel(in_node_features=args.in_node_features, in_edge_features=args.in_edge_features,
                               ns=args.ns, nv=args.nv, sigma_embed_dim=args.sigma_embed_dim, sigma_min=args.sigma_min,
                               sigma_max=args.sigma_max, num_conv_layers=args.num_conv_layers,
                               max_radius=args.max_radius, radius_embed_dim=args.radius_embed_dim,
                               scale_by_sigma=args.scale_by_sigma, use_second_order_repr=args.use_second_order_repr,
                               residual=not args.no_residual, batch_norm=not args.no_batch_norm,
                               parity=not args.no_parity, lambda_embed_dim=args.lambda_embed_dim,
                               l_head_init=args.l_head_init)
