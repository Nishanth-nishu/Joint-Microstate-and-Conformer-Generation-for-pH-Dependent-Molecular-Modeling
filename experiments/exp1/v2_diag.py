"""EXP-1 V2 focused S2 diagnosis (validation molecules only; no training, no test molecule).

Per sample of the §10.1 transform (σ_b as the runs, fixed seed), with r(Δ̂) = ‖P_W(x_t)(x_t + P_W Δ̂ − x_G*′)‖² / rank:
  base         r(0), the S2 reference
  signal       ‖P_W(x_t)(interp(x_R′, x_G*′, t) − x_G*′)‖² / rank   (bridge noise removed)
  noise        ‖P_W(x_t) γ P_W(x_R′) z‖² / rank                     (γ = σ_b √(t(1 − t)))
  mmff         r(x_MMFF − x_t): MMFF94s relaxation of x_t (≤ 200 iterations), a non-learned predictor of the target
               from x_t alone (GEOM targets are xTB, so this is a reference, not an optimum)
  model:<dir>  r(Δ̂) of each trained checkpoint (best_model.pt)
Reported overall and by t quintile, as ratios to base (mean of per-sample values / mean base).

Run (TD venv): python experiments/exp1/v2_diag.py --cache_val <cache_paired.val> --runs <log_dir> ... --out <json>
"""
from __future__ import annotations

import copy
import json
import pickle
import random
from argparse import ArgumentParser, Namespace

import numpy as np
import torch
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem
from torch_geometric.data import Batch

from flexitors.geometry.projector import dihedral_preserving_projector
from flexitors.models.flexitors_score import build_model
from flexitors.training.transforms import PairedJointTransform


def proj_loss(P, rank, v):
    return float(((P @ v.reshape(-1)) ** 2).sum() / max(rank, 1))


def mmff_relax(mol, x, iters=200):
    m = Chem.Mol(mol)
    conf = m.GetConformer()
    for i, p in enumerate(x):
        conf.SetAtomPosition(i, p.tolist())
    props = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
    if props is None:
        return None
    ff = AllChem.MMFFGetMoleculeForceField(m, props)
    ff.Minimize(maxIts=iters)
    return np.array(ff.Positions()).reshape(-1, 3)


def main():
    ap = ArgumentParser()
    ap.add_argument('--cache_val', required=True)
    ap.add_argument('--runs', nargs='*', default=[])
    ap.add_argument('--sigma_b', type=float, default=0.15)
    ap.add_argument('--n', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    val = pickle.load(open(a.cache_val, 'rb'))
    rng = np.random.default_rng(a.seed)
    idx = [int(i) for i in rng.permutation(len(val))]
    models = {}
    for d in a.runs:
        args = Namespace(**yaml.load(open(f'{d}/model_parameters.yml'), Loader=yaml.Loader))
        args.l_head_init = getattr(args, 'l_head_init', 'default')  # runs before the flag existed
        m = build_model(args).eval()
        m.load_state_dict(torch.load(f'{d}/best_model.pt', map_location='cpu'), strict=True)
        models[d] = m
    tf = PairedJointTransform(a.sigma_b)
    interp = tf._dataset.interp_x_torch
    rows = []
    for i in idx:
        if len(rows) >= a.n:
            break
        d0 = val[i]
        mol = d0.mol
        random.seed(a.seed + i); np.random.seed(a.seed + i)
        out = tf(copy.deepcopy(d0))
        if float(out.geo_mask) != 1.0 or out.w_quad_index.shape[1] == 0:
            continue
        # replay the transform's draws to separate the bridge mean from its noise
        random.seed(a.seed + i); np.random.seed(a.seed + i)
        k = random.randrange(len(d0.pos))
        sigma = float(np.exp(np.random.uniform(np.log(tf.sigma_min), np.log(tf.sigma_max))))
        dtau = np.random.normal(0.0, sigma, size=int(np.asarray(d0.edge_mask, dtype=bool).sum()))
        t = float(np.random.uniform())
        assert np.allclose(dtau, out.edge_rotate.numpy()) and abs(t - float(out.node_lambda[0])) < 1e-6
        xt = out.pos.double().numpy()
        xg = out.geo_target.double().numpy()
        quads = out.w_quad_index.T.numpy()
        em = np.asarray(out.edge_mask, dtype=bool)
        rot = out.edge_index.numpy().T[em]
        xR = d0.pos[k].numpy().astype(np.float64)
        xRp = tf._modify(xR, rot, np.asarray(out.mask_rotate, dtype=bool), dtau)
        mean = interp(torch.as_tensor(xRp), torch.as_tensor(xg), t, out.edge_index).numpy()
        P, rank = dihedral_preserving_projector(xt, quads)
        r = dict(i=i, t=t, n_atoms=len(xt), rank=int(rank), base=proj_loss(P, rank, xt - xg),
                 signal=proj_loss(P, rank, mean - xg), noise=proj_loss(P, rank, xt - mean))
        xm = mmff_relax(mol, xt)
        r['mmff'] = proj_loss(P, rank, xt + (P @ (xm - xt).reshape(-1)).reshape(xt.shape) - xg) if xm is not None else None
        if models:
            b = Batch.from_data_list([out])
            with torch.no_grad():
                for dname, m in models.items():
                    dh = m(b.clone()).l_pred.double().numpy()
                    r['model:' + dname.rstrip('/').split('/')[-1]] = proj_loss(P, rank, xt + (P @ dh.reshape(-1)).reshape(xt.shape) - xg)
        rows.append(r)
    keys = [k for k in rows[0] if k in ('signal', 'noise', 'mmff') or k.startswith('model:')]
    q = np.quantile([r['t'] for r in rows], [0.2, 0.4, 0.6, 0.8])

    def summary(rs):
        base = np.mean([r['base'] for r in rs])
        s = dict(n=len(rs), base=base)
        for k in keys:
            v = [r[k] for r in rs if r[k] is not None]
            s[k + '/base'] = float(np.mean(v) / np.mean([r['base'] for r in rs if r[k] is not None])) if v else None
        return s
    R = dict(sigma_b=a.sigma_b, n=len(rows), overall=summary(rows),
             by_t_quintile=[summary([r for r in rows if lo <= r['t'] < hi])
                            for lo, hi in zip([-1.0, *q], [*q, 2.0])],
             t_edges=[float(x) for x in q],
             mmff_failures=sum(r['mmff'] is None for r in rows))
    json.dump(dict(summary=R, rows=rows), open(a.out, 'w'), default=float)
    print(json.dumps(R, indent=1, default=float))


if __name__ == '__main__':
    main()
