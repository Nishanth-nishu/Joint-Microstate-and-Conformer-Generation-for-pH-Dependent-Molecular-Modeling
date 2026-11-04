"""EXP-1 V2 training smoke: checks of a 1-epoch run (validation data only; no test-set molecule is touched).

Criteria (declared here before the run, 2026-10-29; a smoke gate, not a model-selection rule):
  S1 finite       every logged train loss and every validation loss is finite
  S2 learning     validation geometry loss < 0.9 × its Δ̂ = 0 base, and the mean train geometry loss of the last 5
                  logged iterations < that of the first 5
  S3 gradients    l_head, conv_layers and embeddings have finite, non-zero gradient norms at every logged iteration;
                  torsion_head too when w_τ > 0 (with w_τ = 0 it must be exactly 0: checked by the CPU resume run)
  S4 checkpoint   best_model.pt and last_model.pt load strictly into a freshly built model; a save → load round trip
                  gives bit-identical outputs; last_model.pt carries optimizer state; the CPU resume run continues at
                  epoch 1 (checked from its epochs.jsonl in the sbatch script)
  S5 determinism  the τ-frozen geometry read-out (§10.2) on 10 validation molecules is bit-identical across two CPU
                  runs; the GPU repeat difference is reported, not gated (GPU kernels are not deterministic)
Also reported, not gated: read-out RMSD to the pure-L target x_G* on 200 validation pairs (seed = the cache's matched
RDKit conformer, NOT a fresh ETKDG embedding), GPU time per epoch.

Run (TD venv): python experiments/exp1/v2_smoke.py --log_dir <trained dir> --cache_val <cache_paired.val> --out <json>
"""
from __future__ import annotations

import io
import json
import math
import pickle
from argparse import ArgumentParser, Namespace

import numpy as np
import torch
import yaml
from torch_geometric.data import Batch

from flexitors.geometry.internal import dihedral, wrap_angle
from flexitors.geometry.l_step import Guard, kabsch_to, restore_dihedrals, w_quadruples
from flexitors.models.flexitors_score import build_model
from flexitors.sampling.l_readout import l_readout
from flexitors.training.transforms import PairedJointTransform


def rows(path):
    return [json.loads(x) for x in open(path)]


def molecule_parts(d, k):
    z = np.array([a.GetAtomicNum() for a in d.mol.GetAtoms()])
    ei = d.edge_index.numpy()
    nb = [[] for _ in z]
    for u, v in ei[:, ::2].T:
        nb[u].append(int(v)); nb[v].append(int(u))
    xR = d.pos[k].numpy().astype(np.float64)
    quads, kept, _, _ = w_quadruples(z, ei.T[np.asarray(d.edge_mask, dtype=bool)], nb, xR)
    sides = np.asarray(d.mask_rotate, dtype=bool)[kept]
    heavy = np.where(z > 1)[0]
    xGs = restore_dihedrals(d.gt_pos[k].numpy().astype(np.float64), quads, sides, dihedral(xR, quads)) if len(quads) \
        else d.gt_pos[k].numpy().astype(np.float64)
    return xR, kabsch_to(xGs, xR, heavy), quads, sides, heavy, Guard(d.mol, ei[:, ::2].T, xR)


def rmsd(x, y, idx):
    return float(np.sqrt(((kabsch_to(x, y, idx)[idx] - y[idx]) ** 2).sum(1).mean()))


def ring_rmsd(mol, x, y):
    q = [(r[i], r[(i + 1) % len(r)], r[(i + 2) % len(r)], r[(i + 3) % len(r)])
         for r in mol.GetRingInfo().AtomRings() if len(r) >= 4 for i in range(len(r))]
    if not q:
        return None
    q = np.array(q)
    return float(np.degrees(np.sqrt(np.mean(wrap_angle(dihedral(x, q) - dihedral(y, q)) ** 2))))


def main():
    ap = ArgumentParser()
    ap.add_argument('--log_dir', required=True)
    ap.add_argument('--cache_val', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--n_readout', type=int, default=200)
    a = ap.parse_args()
    R = {}
    tr, ep = rows(f'{a.log_dir}/train_log.jsonl'), rows(f'{a.log_dir}/epochs.jsonl')
    # yaml.Loader: runs before the loop.py fix dumped torch.__version__ as a python object (our own file)
    args = Namespace(**yaml.load(open(f'{a.log_dir}/model_parameters.yml'), Loader=yaml.Loader))
    vals = [e[k] for e in ep for k in ('val_td', 'val_geo', 'val_objective')]
    R['S1_finite'] = all(math.isfinite(r['loss']) for r in tr) and all(math.isfinite(v) for v in vals)
    last = ep[-1]
    first5, last5 = np.mean([r['geo'] for r in tr[:5]]), np.mean([r['geo'] for r in tr[-5:]])
    R['S2_detail'] = dict(val_geo=last['val_geo'], val_geo_base=last['val_geo_base'], val_td=last['val_td'],
                          val_td_base=last['val_td_base'], train_geo_first5=first5, train_geo_last5=last5)
    R['S2_learning'] = bool(last['val_geo'] < 0.9 * last['val_geo_base'] and last5 < first5)
    groups = ['l_head', 'conv_layers', 'embeddings'] + (['torsion_head'] if args.w_tau > 0 else [])
    R['S3_detail'] = {g: [min(r[f'grad_{g}'] for r in tr), max(r[f'grad_{g}'] for r in tr)] for g in
                      ('l_head', 'conv_layers', 'embeddings', 'torsion_head')}
    R['S3_gradients'] = all(math.isfinite(r[f'grad_{g}']) and r[f'grad_{g}'] > 0 for r in tr for g in groups)
    # S4: strict loads + round trip
    val = pickle.load(open(a.cache_val, 'rb'))
    np.random.seed(0)
    import random; random.seed(0)
    tf = PairedJointTransform(args.sigma_b, args.sigma_min, args.sigma_max)
    import copy
    b = Batch.from_data_list([tf(copy.deepcopy(val[i])) for i in range(16)])
    m1, m2 = build_model(args).eval(), build_model(args).eval()
    m1.load_state_dict(torch.load(f'{a.log_dir}/best_model.pt', map_location='cpu'), strict=True)
    s = torch.load(f'{a.log_dir}/last_model.pt', map_location='cpu')
    build_model(args).load_state_dict(s['model'], strict=True)
    buf = io.BytesIO(); torch.save(m1.state_dict(), buf); buf.seek(0)
    m2.load_state_dict(torch.load(buf), strict=True)
    with torch.no_grad():
        o1, o2 = m1(b.clone()), m2(b.clone())
    R['S4_detail'] = dict(optimizer_steps=int(next(iter(s['optimizer']['state'].values()))['step']),
                          epoch=s['epoch'], best_val=s['best_val'])
    R['S4_checkpoint'] = bool(torch.equal(o1.edge_pred, o2.edge_pred) and torch.equal(o1.l_pred, o2.l_pred)
                              and R['S4_detail']['optimizer_steps'] > 0)
    # S5: determinism of the read-out, and read-out quality (informative)
    picks = [(i, 0) for i in range(len(val)) if bool(val[i].pair_ok[0])][:a.n_readout]
    det_cpu, det_gpu, qual = [], [], []
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    mg = build_model(args).to(dev).eval(); mg.load_state_dict(m1.state_dict())
    for j, (i, k) in enumerate(picks):
        d = val[i]
        xR, xGs, quads, sides, heavy, guard = molecule_parts(d, k)
        x1, c1 = l_readout(m1, d, xR, quads, sides, guard, args.sigma_min)
        if j < 10:
            x2, _ = l_readout(m1, d, xR, quads, sides, guard, args.sigma_min)
            det_cpu.append(float(np.abs(x1 - x2).max()))
            g1, _ = l_readout(mg, d, xR, quads, sides, guard, args.sigma_min, device=dev)
            g2, _ = l_readout(mg, d, xR, quads, sides, guard, args.sigma_min, device=dev)
            det_gpu.append(float(np.abs(g1 - g2).max()))
        qual.append(dict(rmsd_seed=rmsd(xR, xGs, heavy), rmsd_readout=rmsd(x1, xGs, heavy),
                         ring_seed=ring_rmsd(d.mol, xR, xGs), ring_readout=ring_rmsd(d.mol, x1, xGs), skips=c1['skips']))
    R['S5_detail'] = dict(cpu_repeat_max_abs_diff=max(det_cpu), gpu_repeat_max_abs_diff=max(det_gpu), device=dev)
    R['S5_determinism'] = max(det_cpu) == 0.0
    q = [r for r in qual]
    rs = [r for r in q if r['ring_seed'] is not None]
    R['readout_informative'] = dict(
        n=len(q), median_rmsd_seed_to_target=float(np.median([r['rmsd_seed'] for r in q])),
        median_rmsd_readout_to_target=float(np.median([r['rmsd_readout'] for r in q])),
        share_readout_closer=float(np.mean([r['rmsd_readout'] < r['rmsd_seed'] for r in q])),
        ring_pairs=len(rs), ring_gt10_seed=float(np.mean([r['ring_seed'] > 10 for r in rs])) if rs else None,
        ring_gt10_readout=float(np.mean([r['ring_readout'] > 10 for r in rs])) if rs else None,
        share_with_skip=float(np.mean([r['skips'] > 0 for r in q])),
        note='seed = cache matched RDKit conformer (not fresh ETKDG); validation molecules; not a gate')
    R['epoch_seconds'] = last['train_seconds']
    R['PASS'] = all(R[k] for k in ('S1_finite', 'S2_learning', 'S3_gradients', 'S4_checkpoint', 'S5_determinism'))
    json.dump(R, open(a.out, 'w'), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float))


if __name__ == '__main__':
    main()
