"""EXP-1 V2 training loop (thin, after TD's ``train.py`` / ``utils/training.py:8-49``): torsion DSM loss (weight w_τ)
plus the projected endpoint geometry loss (weight w_geo) on the paired cache, with its own configuration file and
checkpoint directory (``--log_dir``; nothing is written next to TD's checkpoints).

Outputs in ``--log_dir``: ``model_parameters.yml`` (all arguments), ``train_log.jsonl`` (every ``--log_every``
iterations: losses, per-module gradient norms, masked shares), ``epochs.jsonl`` (per epoch: train/val losses and base
losses, timing), ``best_model.pt`` (state dict, lowest validation objective w_τ·L_τ + w_geo·L_geo), ``last_model.pt``
(model, optimizer, scheduler, epoch, best value, python/numpy/torch/CUDA RNG states; written atomically).
Validation noise is reproducible: the validation loader is rebuilt every epoch with a fixed generator seed.

Run (TD venv, FLEXITORS_TD_ROOT = the frozen TD checkout with the torus tables):
  python -m flexitors.training.loop --cache_train <cache_paired.train> --cache_val <cache_paired.val> --log_dir <dir> \
      --sigma_b 0.15 --w_tau 0 --w_geo 1 --n_epochs 1 [--limit_train_iters N --limit_val_iters N --resume]
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import os
import pickle
import random
import time

import numpy as np
import torch
import yaml
from torch_geometric.data import Dataset
from torch_geometric.loader import DataLoader

from flexitors.models.flexitors_score import add_model_args, build_model
from flexitors.training.losses import projected_l_loss, td_loss
from flexitors.training.transforms import PairedJointTransform, seed_numpy_worker

VAL_SEED = 20261029
GRAD_GROUPS = {'l_head': ('l_head.',), 'torsion_head': ('bond_conv.', 'final_linear.', 'final_edge_embedding.'),
               'conv_layers': ('conv_layers.',), 'embeddings': ('node_embedding.', 'edge_embedding.')}


class PairedCache(Dataset):
    """The round-2 paired cache (a pickled list of PyG Data); ``get`` deep-copies like TD's ``ConformerDataset``."""

    def __init__(self, path: str, transform, limit: int = 0):
        super().__init__(transform=transform)
        with open(path, 'rb') as f:
            self.points = pickle.load(f)
        if limit:
            self.points = self.points[:limit]

    def len(self):
        return len(self.points)

    def get(self, idx):
        return copy.deepcopy(self.points[idx])


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache_train', required=True)
    ap.add_argument('--cache_val', required=True)
    ap.add_argument('--log_dir', required=True)
    ap.add_argument('--sigma_b', type=float, default=0.15, help='γ amplitude (Å), EXP-1 §2.5 sweep {0.05, 0.15, 0.3}')
    ap.add_argument('--w_tau', type=float, default=1.0)
    ap.add_argument('--w_geo', type=float, default=1.0)
    ap.add_argument('--geo_projector', choices=('on', 'off'), default='on',
                    help='off: ablation without the dihedral-preserving projector in the loss (P_W = I)')
    ap.add_argument('--geo_norm', choices=('none', 'base'), default='none',
                    help='base: optimise L_geo / (its Δ̂ = 0 base, detached, per batch), so it is O(1) like L_τ')
    ap.add_argument('--n_epochs', type=int, default=1)
    ap.add_argument('--batch_size', type=int, default=32)
    ap.add_argument('--lr', type=float, default=1e-3)
    ap.add_argument('--scheduler_patience', type=int, default=20)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--num_workers', type=int, default=6)
    ap.add_argument('--limit_train_iters', type=int, default=0)
    ap.add_argument('--limit_val_iters', type=int, default=0)
    ap.add_argument('--limit_train_mols', type=int, default=0)
    ap.add_argument('--log_every', type=int, default=50)
    ap.add_argument('--resume', action='store_true')
    add_model_args(ap)
    return ap.parse_args(argv)


def grad_norms(model) -> dict:
    out = {}
    for g, prefixes in GRAD_GROUPS.items():
        sq = [p.grad.detach().double().pow(2).sum() for n, p in model.named_parameters()
              if p.grad is not None and n.startswith(prefixes)]
        out[g] = float(torch.stack(sq).sum().sqrt()) if sq else 0.0
    return out


def geo_objective(lg, bg, geo_norm):
    """The optimised geometry term: raw L_geo, or L_geo / base (base detached; raw when no molecule is active)."""
    if geo_norm == 'base' and float(bg) > 0:
        return lg / bg.detach()
    return lg


def run_batch(model, data, args):
    data = model(data)
    lt, bt = td_loss(data)
    lg, bg, n_geo = projected_l_loss(data, project=args.geo_projector == 'on')
    go = geo_objective(lg, bg, args.geo_norm)
    loss = args.w_tau * lt + args.w_geo * go.float()
    return loss, dict(loss=float(loss), td=float(lt), td_base=float(bt), geo=float(lg), geo_base=float(bg),
                      geo_obj=float(go), n_geo=n_geo,
                      pair_masked=float(data.pair_masked.mean()), chirality_masked=float(data.chirality_masked.mean()),
                      redraws=float(data.n_chirality_redraws.float().mean()))


def _mean(rows, key):
    return float(np.mean([r[key] for r in rows])) if rows else float('nan')


def rng_state() -> dict:
    st = dict(py_rng=random.getstate(), np_rng=np.random.get_state(), torch_rng=torch.get_rng_state())
    if torch.cuda.is_available():
        st['cuda_rng'] = torch.cuda.get_rng_state_all()
    return st


def save_atomic(obj, path):
    torch.save(obj, path + '.tmp')
    os.replace(path + '.tmp', path)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.log_dir, exist_ok=True)
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = build_model(args).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, mode='min', factor=0.7, patience=args.scheduler_patience,
                                                       min_lr=args.lr / 100)
    tf = PairedJointTransform(args.sigma_b, args.sigma_min, args.sigma_max)
    train_ds = PairedCache(args.cache_train, tf, args.limit_train_mols)
    val_ds = PairedCache(args.cache_val, tf)
    loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=args.num_workers,
                        worker_init_fn=seed_numpy_worker, persistent_workers=False)
    start, best = 0, math.inf
    last = os.path.join(args.log_dir, 'last_model.pt')
    if args.resume and os.path.exists(last):
        s = torch.load(last, map_location='cpu')
        model.load_state_dict(s['model'], strict=True)
        opt.load_state_dict(s['optimizer'])
        sched.load_state_dict(s['scheduler'])
        random.setstate(s['py_rng']); np.random.set_state(s['np_rng']); torch.set_rng_state(s['torch_rng'])
        if 'cuda_rng' in s and torch.cuda.is_available():
            torch.cuda.set_rng_state_all(s['cuda_rng'])
        start, best = s['epoch'] + 1, s['best_val']
        print(f'RESUMED at epoch {start} (best val {best})', flush=True)
    with open(os.path.join(args.log_dir, 'model_parameters.yml'), 'w') as f:
        yaml.dump(dict(vars(args), n_parameters=sum(p.numel() for p in model.parameters()),
                       device=str(device), torch=str(torch.__version__)), f)
    print(f'model parameters {sum(p.numel() for p in model.parameters())}, train mols {len(train_ds)}, '
          f'val mols {len(val_ds)}, device {device}', flush=True)
    log = open(os.path.join(args.log_dir, 'train_log.jsonl'), 'a')
    ep_log = open(os.path.join(args.log_dir, 'epochs.jsonl'), 'a')
    for epoch in range(start, args.n_epochs):
        model.train()
        t0, rows, it = time.time(), [], 0
        for data in loader:
            data = data.to(device)
            opt.zero_grad()
            loss, rec = run_batch(model, data, args)
            if not math.isfinite(rec['loss']):
                raise FloatingPointError(f'non-finite loss at epoch {epoch} iteration {it}: {rec}')
            loss.backward()
            g = grad_norms(model)
            opt.step()
            rec.update(epoch=epoch, it=it, **{f'grad_{k}': v for k, v in g.items()})
            rows.append(rec)
            if it % args.log_every == 0:
                log.write(json.dumps(rec) + '\n'); log.flush()
                print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)
            it += 1
            if args.limit_train_iters and it >= args.limit_train_iters:
                break
        t_train = time.time() - t0
        model.eval()
        vrows = []
        vloader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers,
                             worker_init_fn=seed_numpy_worker, generator=torch.Generator().manual_seed(VAL_SEED))
        state = (random.getstate(), np.random.get_state())
        random.seed(VAL_SEED); np.random.seed(VAL_SEED)  # workers=0 path draws from the main-process RNGs
        with torch.no_grad():
            for j, data in enumerate(vloader):
                vrows.append(run_batch(model, data.to(device), args)[1])
                if args.limit_val_iters and j + 1 >= args.limit_val_iters:
                    break
        random.setstate(state[0]); np.random.set_state(state[1])
        val = args.w_tau * _mean(vrows, 'td') + args.w_geo * _mean(vrows, 'geo_obj')
        sched.step(val)
        rec = dict(epoch=epoch, train_seconds=t_train, iters=it, val_objective=val,
                   **{f'train_{k}': _mean(rows, k) for k in ('loss', 'td', 'td_base', 'geo', 'geo_base', 'pair_masked',
                                                            'chirality_masked', 'redraws')},
                   **{f'val_{k}': _mean(vrows, k) for k in ('td', 'td_base', 'geo', 'geo_base')},
                   val_batches=len(vrows), lr=opt.param_groups[0]['lr'])
        if val <= best:
            best = val
            save_atomic(model.state_dict(), os.path.join(args.log_dir, 'best_model.pt'))
            rec['best'] = True
        save_atomic(dict(epoch=epoch, model=model.state_dict(), optimizer=opt.state_dict(), scheduler=sched.state_dict(),
                         best_val=best, **rng_state()), last)
        ep_log.write(json.dumps(rec) + '\n'); ep_log.flush()
        print('EPOCH ' + json.dumps(rec), flush=True)


if __name__ == '__main__':
    main()
