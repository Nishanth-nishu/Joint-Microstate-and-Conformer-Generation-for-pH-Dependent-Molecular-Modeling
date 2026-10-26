"""Pre-registered analysis for the N1 pilot (research/experiments/N1_pilot/PREREG.md)."""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPORTED_MAE = {'sampl6': 0.493, 'sampl7': 0.554, 'sampl8': (0.619, 0.631), 'novartis': 0.810}
ARMS = ['A3D-mean', 'A2D', 'A3D-only-mean', 'A3D-lse', 'A3D-min', 'A3D-single']


def boot_diff(err_a, err_b, n=10000, seed=0):
    """Paired molecule-level bootstrap of mean(|err_a|) - mean(|err_b|)."""
    rng = np.random.default_rng(seed)
    d = np.abs(err_a) - np.abs(err_b)
    idx = rng.integers(0, len(d), size=(n, len(d)))
    bs = d[idx].mean(axis=1)
    return float(d.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def boot_did(ds, dc, n=10000, seed=1):
    rng = np.random.default_rng(seed)
    bs = (ds[rng.integers(0, len(ds), size=(n, len(ds)))].mean(axis=1)
          - dc[rng.integers(0, len(dc), size=(n, len(dc)))].mean(axis=1))
    return float(ds.mean() - dc.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True)
    ap.add_argument('--stratum_col', default='stratum_v2')
    a = ap.parse_args()
    df = pd.read_csv(os.path.join(a.run, 'transitions.csv'))
    df['stratum'] = df[a.stratum_col]
    n_fail = int(df['error'].notna().sum()) if 'error' in df else 0
    df = df[df['error'].isna()] if 'error' in df else df
    for arm in ARMS:
        df[f'err_{arm}'] = df[arm] - df['target']

    # 2D padding check: transitions where some 3D slot equals the 2D copy (embedding failure)
    pad = 0
    for f in glob.glob(os.path.join(a.run, 'conformer_free_energies_fold0.csv')):
        c = pd.read_csv(f)
        p = c.pivot_table(index='smiles_c', columns='conformer_index', values='conformer_free_energy')
        pad = int((p.iloc[:, :10].eq(p.iloc[:, 10], axis=0)).any(axis=1).sum())

    res = {'n_transitions': int(len(df)), 'n_failed': n_fail, 'microstates_with_2d_padding_fold0': pad}
    df['group'] = df['set'].str.replace('_acid', '').str.replace('_base', '')
    res['mae_by_set'] = {g: {arm: round(float(x[f'err_{arm}'].abs().mean()), 3) for arm in ARMS} | {'n': int(len(x))}
                         for g, x in df.groupby('group')}
    gate = {}
    for g, rep in REPORTED_MAE.items():
        if g in res['mae_by_set']:
            m = res['mae_by_set'][g]['A3D-mean']
            lo, hi = (rep, rep) if isinstance(rep, float) else rep
            gate[g] = {'ours': m, 'reported': rep, 'pass': bool(lo - 0.10 <= m <= hi + 0.10)}
    res['sanity_gate'] = gate
    res['sanity_gate_pass'] = all(v['pass'] for v in gate.values())

    res['stratum_counts'] = df['stratum'].value_counts().to_dict()
    strata = {}
    for s, x in df.groupby('stratum'):
        e = {arm: x[f'err_{arm}'].to_numpy() for arm in ARMS}
        strata[s] = {
            'n': int(len(x)),
            'mae': {arm: round(float(np.abs(e[arm]).mean()), 3) for arm in ARMS},
            'delta_2D_minus_3Donly': boot_diff(e['A2D'], e['A3D-only-mean']),
            'delta_lse_minus_mean': boot_diff(e['A3D-lse'], e['A3D-only-mean']),
            'delta_min_minus_mean': boot_diff(e['A3D-min'], e['A3D-only-mean']),
            'delta_single_minus_mean': boot_diff(e['A3D-single'], e['A3D-only-mean']),
            'spearman_confstd_vs_abserr': list(map(float, spearmanr(x['conf_pka_std'],
                                                                     np.abs(e['A3D-only-mean'])))),
            'median_conf_pka_std': float(x['conf_pka_std'].median()),
        }
    res['strata'] = strata
    if 'sensitive' in strata and 'control' in strata:
        xs, xc = df[df.stratum == 'sensitive'], df[df.stratum == 'control']
        ds = (xs['err_A2D'].abs() - xs['err_A3D-only-mean'].abs()).to_numpy()
        dc = (xc['err_A2D'].abs() - xc['err_A3D-only-mean'].abs()).to_numpy()
        res['did_sensitive_minus_control'] = boot_did(ds, dc)
        d, lo, hi = strata['sensitive']['delta_2D_minus_3Donly']
        dd, dlo, dhi = res['did_sensitive_minus_control']
        if len(xs) < 30:
            verdict = 'UNDERPOWERED (sensitive n<30): extend to stage 2 (IUPAC clean subset)'
        elif lo > 0 and dlo > 0:
            verdict = 'SUPPORTS N1'
        else:
            verdict = 'N1 FALSIFIED on this pilot (CI of delta_s or DiD includes/below 0)'
        res['prereg_verdict'] = verdict if res['sanity_gate_pass'] else 'SANITY GATE FAILED: no conclusion'
    with open(os.path.join(a.run, 'analysis.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1), flush=True)


if __name__ == '__main__':
    main()
