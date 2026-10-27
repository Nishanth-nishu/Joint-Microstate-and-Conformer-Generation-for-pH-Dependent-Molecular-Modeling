"""Pre-registered analysis for the tautomer pilot (research/experiments/T_pilot/PREREG.md)."""
import argparse
import json
import os

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import spearmanr

from jmcg.experiments.n1_analyze import REPORTED_MAE

TAUT_ARMS = ['R0', 'R-pKa', 'R-joint', 'R-joint+E', 'R-taut']


def cluster_boot(diff, groups, n=10000, seed=0):
    """Mean of per-item diff with a bootstrap that resamples whole groups (scaffolds)."""
    rng = np.random.default_rng(seed)
    g = pd.Series(diff).groupby(np.asarray(groups))
    s, c = g.sum().to_numpy(), g.size().to_numpy()
    idx = rng.integers(0, len(s), size=(n, len(s)))
    bs = s[idx].sum(1) / c[idx].sum(1)
    return float(np.mean(diff)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def slope(pred, true, groups, n=2000, seed=0):
    """OLS slope of predicted on measured log K, with a group (scaffold) bootstrap CI."""
    pred, true, groups = np.asarray(pred), np.asarray(true), np.asarray(groups)
    fit = lambda p, t: float(np.polyfit(t, p, 1)[0])
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    members = {g: np.where(groups == g)[0] for g in uniq}
    bs = []
    for _ in range(n):
        idx = np.concatenate([members[g] for g in rng.choice(uniq, len(uniq))])
        bs.append(fit(pred[idx], true[idx]))
    return fit(pred, true), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def stage0(path, dwar_tsv):
    df = pd.read_csv(path)
    df = df[df['error'].isna()] if 'error' in df else df
    seen = set()
    for r in pd.read_csv(dwar_tsv, sep='\t', index_col=0)['SMILES']:
        for s in r.replace('>>', ',').split(','):
            m = Chem.MolFromSmiles(s)
            if m is not None:
                seen.add(Chem.MolToSmiles(m))
    df['both_in_dwar_microstates'] = [int(Chem.MolToSmiles(Chem.MolFromSmiles(a)) in seen
                                          and Chem.MolToSmiles(Chem.MolFromSmiles(b)) in seen)
                                      for a, b in zip(df.smiles_t1, df.smiles_t2)]
    res = {}
    for kind, x in df.groupby('label_kind'):
        e = x['logk_A3D-mean'] - x['log_k']
        big = x[x.log_k.abs() >= 1]
        res[kind] = {'n': int(len(x)), 'mae': round(float(e.abs().mean()), 3), 'rmse': round(float(np.sqrt((e ** 2).mean())), 3),
                     'mean_signed': round(float(e.mean()), 3),
                     'sign_acc_|logK|>=1': round(float((np.sign(big['logk_A3D-mean']) == np.sign(big.log_k)).mean()), 3),
                     'n_|logK|>=1': int(len(big)),
                     'mae_2D': round(float((x['logk_A2D'] - x.log_k).abs().mean()), 3),
                     'mae_by_both_in_dwar': {int(k): [int(len(v)), round(float((v['logk_A3D-mean'] - v.log_k).abs().mean()), 3)]
                                             for k, v in x.groupby('both_in_dwar_microstates')},
                     'mae_predict_zero': round(float(x.log_k.abs().mean()), 3),
                     'slope_pred_on_meas': slope(x['logk_A3D-mean'], x.log_k,
                                                 [sc if isinstance(sc, str) and sc else t for sc, t in zip(x.scaffold, x.smiles_t1)])}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage0', required=True, help='taut_zeroshot pairs.csv')
    ap.add_argument('--fit', required=True, help='taut_headfit output dir')
    ap.add_argument('--dwar', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    res = {'stage0_zero_shot_5fold': stage0(a.stage0, a.dwar)}

    oof = pd.read_csv(os.path.join(a.fit, 'tautomer_oof.csv'))
    ext = pd.read_csv(os.path.join(a.fit, 'external_pka.csv'))
    res['n_measured_pairs'] = int(len(oof))
    res['n_groups'] = int(oof.group.nunique())
    res['taut_oof_mae'] = {arm: round(float((oof[f'logk_{arm}'] - oof.log_k).abs().mean()), 3) for arm in TAUT_ARMS}
    res['taut_oof_mae_per_seed'] = {arm: [round(float((oof[c] - oof.log_k).abs().mean()), 3)
                                          for c in oof.columns if c.startswith(f'logk_{arm}_seed')] for arm in TAUT_ARMS}
    res['taut_oof_mae_by_stratum'] = {s: {arm: round(float((x[f'logk_{arm}'] - x.log_k).abs().mean()), 3) for arm in TAUT_ARMS}
                                      | {'n': int(len(x))} for s, x in oof.groupby('stratum_v2')}
    ae = {arm: (oof[f'logk_{arm}'] - oof.log_k).abs().to_numpy() for arm in TAUT_ARMS}
    res['taut_oof_slope'] = {arm: slope(oof[f'logk_{arm}'], oof.log_k, oof.group) for arm in TAUT_ARMS}
    res['taut_oof_mae_predict_zero'] = round(float(oof.log_k.abs().mean()), 3)
    res['delta_taut_tautonly_minus_joint'] = cluster_boot(ae['R-taut'] - ae['R-joint'], oof.group)
    res['delta_taut_pKa_minus_joint'] = cluster_boot(ae['R-pKa'] - ae['R-joint'], oof.group)
    res['delta_taut_R0_minus_joint'] = cluster_boot(ae['R0'] - ae['R-joint'], oof.group)
    res['delta_taut_joint_minus_jointE'] = cluster_boot(ae['R-joint'] - ae['R-joint+E'], oof.group)

    ext['group'] = ext.set.str.replace('_acid', '').str.replace('_base', '')
    pe = {arm: (ext[f'pka_{arm}'] - ext.target).abs().to_numpy() for arm in TAUT_ARMS}
    res['ext_pka_mae'] = {arm: {g: round(float((x[f'pka_{arm}'] - x.target).abs().mean()), 3) for g, x in ext.groupby('group')}
                          | {'all': round(float(pe[arm].mean()), 3)} for arm in TAUT_ARMS}
    keys = ext.set.astype(str) + ':' + ext.i.astype(str)
    res['ext_pka_joint_minus_pKa'] = cluster_boot(pe['R-joint'] - pe['R-pKa'], keys)

    # gate: fold-0 released head reproduces the published external MAEs within 0.15 (single fold, not the 5-fold mean)
    gate = {}
    for g, rep in REPORTED_MAE.items():
        x = ext[ext.group == g] if g == 'novartis' else ext[ext.set == g]
        m =float((x['pka_R0'] - x.target).abs().mean())
        lo, hi = (rep, rep) if isinstance(rep, float) else rep
        gate[g] = {'ours_fold0': round(m, 3), 'reported': rep, 'pass': bool(lo - 0.15 <= m <= hi + 0.15)}
    res['sanity_gate'] = gate
    res['sanity_gate_pass'] = all(v['pass'] for v in gate.values())

    # primary (judge-narrowed): joint vs tautomer-only fit of the same head; pKa non-inferiority vs pKa-only refit
    d, lo, hi = res['delta_taut_tautonly_minus_joint']
    nd, nlo, nhi = res['ext_pka_joint_minus_pKa']
    mae_joint = res['taut_oof_mae']['R-joint']
    if not res['sanity_gate_pass']:
        verdict = 'SANITY GATE FAILED: no conclusion'
    elif lo > 0 and nhi <= 0.05 and mae_joint <= 1.0:
        verdict = 'SUPPORTS T2'
    else:
        why = [w for w, bad in [('delta(taut-only - joint) CI includes/below 0', lo <= 0),
                                ('macro-pKa non-inferiority fails (upper CI > 0.05)', nhi > 0.05),
                                ('joint OOF MAE > 1.0 log K', mae_joint > 1.0)] if bad]
        verdict = 'T2 FALSIFIED: ' + '; '.join(why)
    res['prereg_verdict'] = verdict
    s0 = res['stage0_zero_shot_5fold'].get('M', {})
    b, blo, bhi = s0.get('slope_pred_on_meas', (np.nan, np.nan, np.nan))
    t1_fail = [w for w, bad in [('zero-shot MAE < 1.0', s0.get('mae', 99) < 1.0),
                                ('slope CI includes 1', blo <= 1 <= bhi),
                                ('not compression: MAE not >= predict-zero + 0.3 while slope >= 0.7',
                                 s0.get('mae', 0) < s0.get('mae_predict_zero', 0) + 0.3 and b >= 0.7)] if bad]
    res['T1_verdict'] = 'T1 FALSIFIED: ' + '; '.join(t1_fail) if t1_fail else 'SUPPORTS T1 (compression)'
    with open(a.out, 'w') as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1), flush=True)


if __name__ == '__main__':
    main()
