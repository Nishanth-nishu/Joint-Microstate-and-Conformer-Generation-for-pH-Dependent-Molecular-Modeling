"""Tautomer pilot stage 2: refit the Uni-pKa free-energy head on cached features (research/experiments/T_pilot/PREREG.md).

G(s) = mean over the 11 coordinate sets of head(CLS feature)   (released inference protocol)
macro-pKa = (LSE_A(-G) - LSE_B(-G)) / ln10 + PKA_MEAN;  log10 K_taut = log10([T2]/[T1]) = (G_T1 - G_T2) / ln10
Arms: R0 released head; R-pKa (Dwar-iBonD only); R-joint (+ measured Tautobase train-fold pairs);
R-joint+E (+ 'E' pairs, secondary); R-taut (tautomers only, secondary). Tautomer metrics are out-of-fold over 5
scaffold-grouped folds of the measured pairs; macro-pKa metrics on the 5 external sets (fold models averaged).
"""
import argparse
import copy
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import torch

from jmcg.experiments import n1_pilot as n1
from jmcg.data.build_sensitive_set import MEASURED, build_tautobase

LN10 = math.log(10)
HP = {'lr': 1e-4, 'epochs': 20, 'pka_batch': 64, 'taut_batch': 16, 'lam': 1.0, 'seeds': (0, 1, 2), 'cv_folds': 5}


class Store:
    """Cached features, one (11, d) block per canonical microstate SMILES."""

    def __init__(self, path):
        z = np.load(path)
        order = np.lexsort((z['conformer_index'], z['smiles_c']))
        smi, feats = z['smiles_c'][order], z['feats'][order].astype(np.float32)
        uniq, start = np.unique(smi, return_index=True)
        self.idx = {s: i for i, s in enumerate(uniq)}
        self.x = torch.from_numpy(feats.reshape(len(uniq), -1, feats.shape[1]))

    def ids(self, smiles):
        return [self.idx[n1.canon(s)] for s in smiles]


def energies(head, x):
    n, c, d = x.shape
    return head(x.reshape(n * c, 1, d)).view(n, c).mean(1)


def pka(head, store, trans):
    """trans: list of (acid_ids, base_ids). Returns tensor of predicted macro-pKa."""
    ids = sorted({i for a, b in trans for i in a + b})
    pos = {i: k for k, i in enumerate(ids)}
    g = energies(head, store.x[ids])
    out = [(torch.logsumexp(-g[[pos[i] for i in a]], 0) - torch.logsumexp(-g[[pos[i] for i in b]], 0)) / LN10 + n1.PKA_MEAN
           for a, b in trans]
    return torch.stack(out)


def logk(head, store, pairs):
    ids = sorted({i for p in pairs for i in p})
    pos = {i: k for k, i in enumerate(ids)}
    g = energies(head, store.x[ids])
    return torch.stack([(g[pos[t1]] - g[pos[t2]]) / LN10 for t1, t2 in pairs])


def load_trans(path, store):
    df = pd.read_csv(path, sep='\t', index_col=0)
    trans, y, miss = [], [], 0
    for smi, t in zip(df['SMILES'], df['TARGET']):
        a, b = smi.split('>>')
        try:
            trans.append((store.ids(a.split(',')), store.ids(b.split(','))))
            y.append(float(t))
        except KeyError:
            miss += 1
    return trans, torch.tensor(y), miss


def fit(head0, store, pka_data, taut_data, seed, use_pka=True):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    head = copy.deepcopy(head0).train()
    opt = torch.optim.Adam(head.parameters(), lr=HP['lr'])
    trans, y = pka_data
    pairs, k = taut_data if taut_data else ([], None)
    # same step count for every arm, so tautomer-only and joint fits see the same number of tautomer batches
    n_steps = math.ceil(len(trans) / HP['pka_batch'])
    for _ in range(HP['epochs']):
        perm = rng.permutation(len(trans)) if use_pka else None
        for s in range(n_steps):
            loss = torch.zeros(())
            if use_pka:
                b = perm[s * HP['pka_batch']:(s + 1) * HP['pka_batch']]
                loss = loss + torch.mean((pka(head, store, [trans[i] for i in b]) - y[b]) ** 2)
            if pairs:
                b = rng.integers(0, len(pairs), HP['taut_batch'])
                loss = loss + HP['lam'] * torch.mean((logk(head, store, [pairs[i] for i in b]) - k[b]) ** 2)
            opt.zero_grad(); loss.backward(); opt.step()
    return head.eval()


def scaffold_folds(groups, n, seed=0):
    uniq = np.array(sorted(set(groups)))
    np.random.default_rng(seed).shuffle(uniq)
    fold_of = {g: i % n for i, g in enumerate(uniq)}
    return np.array([fold_of[g] for g in groups])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features', required=True)
    ap.add_argument('--head', required=True, help='head_released_fold0.pt')
    ap.add_argument('--data', required=True, help='ext/unipka_orig/dataset')
    ap.add_argument('--tautobase', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--smoke', action='store_true', help='1 epoch, 1 seed')
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    if a.smoke:
        HP.update(epochs=1, seeds=(0,))
    os.makedirs(a.out, exist_ok=True)

    store = Store(a.features)
    head0 = torch.load(a.head, weights_only=False).eval()
    z = np.load(a.features)
    with torch.no_grad():
        g = head0(torch.from_numpy(z['feats'].astype(np.float32))[:, None, :]).view(-1).numpy()
    print(f'released head on cached features: max |dG| = {np.abs(g - z["g_released"]).max():.2e}', flush=True)

    dwar, y_dwar, miss = load_trans(os.path.join(a.data, 'dwar-iBond.tsv'), store)
    ext = {s: load_trans(os.path.join(a.data, f'{s}.tsv'), store) for s in n1.SETS}
    rows, _ = build_tautobase(a.tautobase, kinds=MEASURED + ('e', 'c'))
    tb = pd.DataFrame(rows)
    tb['ok'] = [all(n1.canon(s) in store.idx for s in (t1, t2)) for t1, t2 in zip(tb.smiles_t1, tb.smiles_t2)]
    meas = tb[(tb.label_kind == 'M') & tb.ok].reset_index(drop=True)
    est = tb[(tb.label_kind == 'E') & tb.ok].reset_index(drop=True)
    # acyclic molecules have an empty Murcko scaffold: group them by their own tautomer pair instead
    groups = [sc if isinstance(sc, str) and sc else min(t1, t2) for sc, t1, t2 in zip(meas.scaffold, meas.smiles_t1, meas.smiles_t2)]
    meas['group'] = groups
    meas['cv'] = scaffold_folds(groups, HP['cv_folds'])
    pm = [tuple(store.ids([t1, t2])) for t1, t2 in zip(meas.smiles_t1, meas.smiles_t2)]
    pe = [tuple(store.ids([t1, t2])) for t1, t2 in zip(est.smiles_t1, est.smiles_t2)]
    km, ke = torch.tensor(meas.log_k.values, dtype=torch.float32), torch.tensor(est.log_k.values, dtype=torch.float32)
    info = {'hp': HP, 'dwar_transitions': len(dwar), 'dwar_missing': miss, 'tautobase_measured': len(meas),
            'tautobase_E': len(est), 'tautobase_missing': int((~tb.ok).sum())}
    print(json.dumps(info), flush=True)

    oof = {arm: np.zeros((len(HP['seeds']), len(meas))) for arm in ['R0', 'R-pKa', 'R-joint', 'R-joint+E', 'R-taut']}
    ext_pred = {arm: {s: [] for s in n1.SETS} for arm in oof}
    with torch.no_grad():
        r0 = logk(head0, store, pm).numpy()
        r0_ext = {s: pka(head0, store, ext[s][0]).numpy() for s in n1.SETS}
    for si, seed in enumerate(HP['seeds']):
        oof['R0'][si] = r0
        for s in n1.SETS:
            ext_pred['R0'][s].append(r0_ext[s])
        h_pka = fit(head0, store, (dwar, y_dwar), None, seed)
        with torch.no_grad():
            oof['R-pKa'][si] = logk(h_pka, store, pm).numpy()
            for s in n1.SETS:
                ext_pred['R-pKa'][s].append(pka(h_pka, store, ext[s][0]).numpy())
        for f in range(HP['cv_folds']):
            tr, te = np.where(meas.cv != f)[0], np.where(meas.cv == f)[0]
            train_m = ([pm[i] for i in tr], km[tr])
            train_me = ([pm[i] for i in tr] + pe, torch.cat([km[tr], ke]))
            for arm, td, use_pka in [('R-joint', train_m, True), ('R-joint+E', train_me, True), ('R-taut', train_m, False)]:
                h = fit(head0, store, (dwar, y_dwar), td, seed, use_pka=use_pka)
                with torch.no_grad():
                    oof[arm][si, te] = logk(h, store, [pm[i] for i in te]).numpy()
                    for s in n1.SETS:
                        ext_pred[arm][s].append(pka(h, store, ext[s][0]).numpy())
            print(f'seed {seed} cv fold {f} done', flush=True)

    out = meas[['row_nr', 'smiles_t1', 'smiles_t2', 'log_k', 'group', 'cv', 'stratum_v2']].copy()
    for arm, v in oof.items():
        out[f'logk_{arm}'] = v.mean(0)
        for si, seed in enumerate(HP['seeds']):
            out[f'logk_{arm}_seed{seed}'] = v[si]
    out.to_csv(os.path.join(a.out, 'tautomer_oof.csv'), index=False)
    ext_rows = []
    for s in n1.SETS:
        for i, t in enumerate(ext[s][1].numpy()):
            ext_rows.append({'set': s, 'i': i, 'target': float(t),
                             **{f'pka_{arm}': float(np.mean([p[i] for p in ext_pred[arm][s]])) for arm in ext_pred}})
    pd.DataFrame(ext_rows).to_csv(os.path.join(a.out, 'external_pka.csv'), index=False)
    with open(os.path.join(a.out, 'fit_info.json'), 'w') as fh:
        json.dump(info, fh, indent=1)


if __name__ == '__main__':
    sys.exit(main())
