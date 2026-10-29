"""E2/E3: conformer-coupled state free energies (joint prototype) vs sequential baseline on measured Tautobase pairs.
Pre-registered in research/experiments/E2_joint/PREREG.md (read it first; this file implements it).

Run (jmcg venv):
  python -m jmcg.joint.e2_tautomer features --data data_manifests/canonical_v1 --out <dir> [--procs 8]
  python -m jmcg.joint.e2_tautomer fit --data data_manifests/canonical_v1 --features <dir>/features.csv \
         --energies <E1 dir>/energies.csv --out <dir>
Arms: S0 (released, no fit), S1 (a), J1 (a + theta.phi), C1 (a + psi, 2D control), J1-phi_k ablations.
log10 K = (G(T1) - G(T2)) / ln10 with G = a G_U + theta.phi: every column is a T1-minus-T2 difference, no intercept,
so swapping T1/T2 negates the prediction. Columns are scaled by their training-fold RMS (no centring, which would break
that antisymmetry); ridge lambda = 1e-3 on the scaled coefficients. Out-of-fold over the T pilot's 5 scaffold folds.
"""
import argparse
import json
import math
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from jmcg.joint import features as ft

LN10 = math.log(10)
LAM = 1e-3
ARMS = {'S1': ['dGU'], 'J1': ['dGU', 'dphi1', 'dphi2', 'dphi3'], 'C1': ['dGU', 'dpsi1', 'dpsi2', 'dpsi3'],
        'J1-phi1': ['dGU', 'dphi2', 'dphi3'], 'J1-phi2': ['dGU', 'dphi1', 'dphi3'], 'J1-phi3': ['dGU', 'dphi1', 'dphi2']}


def scaffold_folds(groups, n, seed=0):
    """Copied from jmcg.experiments.taut_headfit (T pilot folds)."""
    uniq = np.array(sorted(set(groups)))
    np.random.default_rng(seed).shuffle(uniq)
    fold_of = {g: i % n for i, g in enumerate(uniq)}
    return np.array([fold_of[g] for g in groups])


def measured_pairs(data):
    obs = pd.read_csv(os.path.join(data, 'observations.csv'), keep_default_na=False)
    smi = pd.read_csv(os.path.join(data, 'states.csv'), keep_default_na=False).set_index('state_id').smiles
    t = obs[(obs.task == 'tautomer_logk') & (obs.label_status == 'measured')].copy()
    t['t1'], t['t2'] = smi[t.states_A].to_numpy(), smi[t.states_B].to_numpy()
    return t


def _feat(s):
    try:
        r = ft.ensemble_features(s)
    except Exception as e:  # noqa: BLE001 - recorded, pair becomes ineligible
        r = dict(n_generated=0, n_valid=0, failures=[f'{type(e).__name__}: {e}'[:200]])
    r.update(ft.graph_features(s))
    r['smiles'] = s
    r['failures'] = '; '.join(sorted(set(r['failures'])))
    return r


def cmd_features(a):
    os.makedirs(a.out, exist_ok=True)
    t = measured_pairs(a.data)
    states = sorted(set(t.t1) | set(t.t2))
    with Pool(a.procs) as p:
        rows = p.map(_feat, states, chunksize=4)
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'features.csv'), index=False)
    print(f'{len(states)} states; with phi: {sum("phi1" in r for r in rows)}')


def boot_delta(e_a, e_b, groups, n=10000, seed=0):
    """Mean of |e_a| - |e_b| with a cluster bootstrap over groups."""
    d = pd.DataFrame(dict(d=np.abs(e_a) - np.abs(e_b), g=groups)).groupby('g').d
    s, c = d.sum().to_numpy(), d.size().to_numpy()
    idx = np.random.default_rng(seed).integers(0, len(s), size=(n, len(s)))
    bs = s[idx].sum(1) / c[idx].sum(1)
    return dict(mean=float(np.mean(np.abs(e_a) - np.abs(e_b))), ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])


def design(a):
    """(all measured pairs with eligibility, eligible design frame with T1-minus-T2 columns, scaffold group, cv fold)."""
    t = measured_pairs(a.data)
    F = pd.read_csv(a.features).set_index('smiles')
    E = pd.read_csv(a.energies).groupby('smiles').G.mean()
    elig, why = [], []
    for _, r in t.iterrows():
        bad = [s for s in (r.t1, r.t2) if s not in F.index or pd.isna(F.loc[s].get('phi1'))]
        bad += [f'no G_U: {s}' for s in (r.t1, r.t2) if s not in E.index]
        elig.append(not bad)
        why.append('; '.join(f"{s}: {F.loc[s].failures if s in F.index else 'no features'}" for s in bad))
    t['eligible'], t['why'] = elig, why
    d = t[t.eligible].copy()
    d['dGU'] = [(E[x] - E[y]) / LN10 for x, y in zip(d.t1, d.t2)]
    for k in ('phi1', 'phi2', 'phi3', 'psi1', 'psi2', 'psi3'):
        d['d' + k] = [F.loc[x, k] - F.loc[y, k] for x, y in zip(d.t1, d.t2)]
    sc = []
    for s1, s2 in zip(d.t1, d.t2):
        try:
            sc.append(MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s1)))
        except Exception:  # noqa: BLE001
            sc.append('')
    d['group'] = [g or min(x, y) for g, x, y in zip(sc, d.t1, d.t2)]
    d['cv'] = scaffold_folds(d.group.tolist(), 5)
    return t, d


def fit_predict(X, y, tr, te):
    rms = np.sqrt((X[tr] ** 2).mean(0))
    rms[rms == 0] = 1.0
    Xs = X / rms
    beta = np.linalg.solve(Xs[tr].T @ Xs[tr] + LAM * np.eye(X.shape[1]), Xs[tr].T @ y[tr])
    return Xs[te] @ beta, beta / rms


def cmd_fit(a):
    os.makedirs(a.out, exist_ok=True)
    t, d = design(a)
    t[~t.eligible][['obs_id', 't1', 't2', 'why']].to_csv(os.path.join(a.out, 'ineligible.csv'), index=False)
    y = d.value.to_numpy()
    pred = {'S0': d.dGU.to_numpy()}
    coefs = {k: [] for k in ARMS}
    for arm, cols in ARMS.items():
        p = np.zeros(len(d))
        for f in range(5):
            tr, te = d.cv.to_numpy() != f, d.cv.to_numpy() == f
            p[te], beta = fit_predict(d[cols].to_numpy(), y, tr, te)
            coefs[arm].append(dict(zip(cols, beta.tolist())))
        pred[arm] = p
    for k, v in pred.items():
        d['pred_' + k] = v
    d.to_csv(os.path.join(a.out, 'oof_predictions.csv'), index=False)
    err = {k: v - y for k, v in pred.items()}
    res = dict(n_measured=len(t), n_eligible=len(d), n_ineligible=int((~t.eligible).sum()), n_groups=int(d.group.nunique()),
               fold_sizes=d.cv.value_counts().sort_index().tolist(), predict_zero_MAE=float(np.abs(y).mean()))
    res['arms'] = {k: dict(MAE=float(np.abs(e).mean()), RMSE=float(np.sqrt((e ** 2).mean())),
                           slope=float(np.polyfit(y, pred[k], 1)[0])) for k, e in err.items()}
    res['coefficients_per_fold'] = coefs
    final = {}
    for arm, cols in ARMS.items():  # final model per arm on every eligible pair (for the sampler; not evaluated)
        X = d[cols].to_numpy()
        rms = np.sqrt((X ** 2).mean(0))
        rms[rms == 0] = 1.0
        beta = np.linalg.solve((X / rms).T @ (X / rms) + LAM * np.eye(len(cols)), (X / rms).T @ y) / rms
        final[arm] = dict(zip(cols, beta.tolist()))
    res['final_coefficients_all_eligible'] = final
    res['delta_S'] = boot_delta(err['S1'], err['J1'], d.group)
    res['delta_C'] = boot_delta(err['C1'], err['J1'], d.group)
    res['S0_minus_S1'] = boot_delta(err['S0'], err['S1'], d.group)
    res['ablation_minus_J1'] = {k: boot_delta(err[k], err['J1'], d.group) for k in ARMS if k.startswith('J1-')}
    sup = res['delta_S']['ci95'][0] > 0 and res['delta_C']['ci95'][0] > 0
    fail = [n for n in ('delta_S', 'delta_C') if res[n]['ci95'][0] <= 0]
    res['verdict'] = 'SUPPORTED' if sup else f'NOT SUPPORTED (CI lower bound <= 0 for {", ".join(fail)})'
    json.dump(res, open(os.path.join(a.out, 'e2_result.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != 'coefficients_per_fold'}, indent=1))


def cmd_robust(a):
    """EXPLORATORY, NOT PRE-REGISTERED (added after the E2 verdict): checks of the SUPPORTED result.
    R1 canonical split: fit on canonical train+val eligible pairs, score on canonical test eligible pairs.
    R2 fold-seed sensitivity: the pre-registered 5-fold scaffold CV repeated with fold seeds 1..20.
    R3 one-feature head-to-head: dGU+dphi3 (3D polar SASA) vs dGU+dpsi3 (2D TPSA).
    R4 eligibility: S0 and S1 (5-fold CV, seed 0) on all 221 measured pairs vs the eligible subset.
    R5 phi1 duplicates: N_eff counts identical minimised conformers separately (no RMSD de-duplication)."""
    os.makedirs(a.out, exist_ok=True)
    t, d = design(a)
    obs = pd.read_csv(os.path.join(a.data, 'observations.csv'), keep_default_na=False).set_index('obs_id')
    d['split'] = obs.split[d.obs_id].to_numpy()
    y = d.value.to_numpy()
    arms = dict(ARMS, **{'S1+phi3': ['dGU', 'dphi3'], 'S1+psi3': ['dGU', 'dpsi3']})
    R = {}
    tr, te = d.split.isin(['train', 'val']).to_numpy(), (d.split == 'test').to_numpy()
    e = {'S0': d.dGU.to_numpy()[te] - y[te]}
    for k in ('S1', 'J1', 'C1', 'S1+phi3', 'S1+psi3'):
        e[k] = fit_predict(d[arms[k]].to_numpy(), y, tr, te)[0] - y[te]
    g = d.group[te]
    R['R1_canonical_split'] = dict(n_train=int(tr.sum()), n_test=int(te.sum()), n_test_groups=int(g.nunique()),
                                   MAE={k: float(np.abs(v).mean()) for k, v in e.items()},
                                   delta_S=boot_delta(e['S1'], e['J1'], g), delta_C=boot_delta(e['C1'], e['J1'], g))
    rows = []
    for seed in range(1, 21):
        cv = scaffold_folds(d.group.tolist(), 5, seed=seed)
        err = {}
        for k in ('S1', 'J1', 'C1', 'S1+phi3', 'S1+psi3'):
            p = np.zeros(len(d))
            for f in range(5):
                p[cv == f] = fit_predict(d[arms[k]].to_numpy(), y, cv != f, cv == f)[0]
            err[k] = p - y
        dS, dC = boot_delta(err['S1'], err['J1'], d.group, n=2000), boot_delta(err['C1'], err['J1'], d.group, n=2000)
        d3 = boot_delta(err['S1+psi3'], err['S1+phi3'], d.group, n=2000)
        rows.append(dict(seed=seed, dS=dS['mean'], dS_lo=dS['ci95'][0], dC=dC['mean'], dC_lo=dC['ci95'][0],
                         phi3_vs_psi3=d3['mean'], phi3_vs_psi3_lo=d3['ci95'][0]))
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(a.out, 'fold_seed_sensitivity.csv'), index=False)
    R['R2_fold_seeds_1_20'] = dict(n_seeds=len(S), dS_mean_range=[S.dS.min(), S.dS.max()], dC_mean_range=[S.dC.min(), S.dC.max()],
                                   n_supported=int(((S.dS_lo > 0) & (S.dC_lo > 0)).sum()), n_dC_lo_positive=int((S.dC_lo > 0).sum()),
                                   n_dS_lo_positive=int((S.dS_lo > 0).sum()))
    err = {}
    for k in ('S1+phi3', 'S1+psi3'):
        p = np.zeros(len(d))
        for f in range(5):
            p[d.cv.to_numpy() == f] = fit_predict(d[arms[k]].to_numpy(), y, d.cv.to_numpy() != f, d.cv.to_numpy() == f)[0]
        err[k] = p - y
    R['R3_phi3_vs_psi3_seed0'] = dict(MAE={k: float(np.abs(v).mean()) for k, v in err.items()},
                                      psi3_minus_phi3=boot_delta(err['S1+psi3'], err['S1+phi3'], d.group))
    E = pd.read_csv(a.energies).groupby('smiles').G.mean()
    al = t[t.t1.isin(E.index) & t.t2.isin(E.index)].copy()
    al['dGU'] = [(E[x] - E[z]) / LN10 for x, z in zip(al.t1, al.t2)]
    sc = []
    for s1 in al.t1:
        try:
            sc.append(MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s1)))
        except Exception:  # noqa: BLE001
            sc.append('')
    al['group'] = [g or min(x, z) for g, x, z in zip(sc, al.t1, al.t2)]
    cv = scaffold_folds(al.group.tolist(), 5)
    ya = al.value.to_numpy()
    p = np.zeros(len(al))
    for f in range(5):
        p[cv == f] = fit_predict(al[['dGU']].to_numpy(), ya, cv != f, cv == f)[0]
    el = al.eligible.to_numpy()
    R['R4_eligibility'] = dict(n_all=len(al), n_eligible=int(el.sum()),
                               S0_MAE_all=float(np.abs(al.dGU - ya).mean()), S0_MAE_eligible=float(np.abs(al.dGU - ya)[el].mean()),
                               S0_MAE_ineligible=float(np.abs(al.dGU - ya)[~el].mean()),
                               S1_MAE_all=float(np.abs(p - ya).mean()), S1_MAE_eligible=float(np.abs(p - ya)[el].mean()),
                               S1_MAE_ineligible=float(np.abs(p - ya)[~el].mean()),
                               predict_zero_MAE_eligible=float(np.abs(ya[el]).mean()), predict_zero_MAE_ineligible=float(np.abs(ya[~el]).mean()),
                               ineligible_reasons={k: int(v) for k, v in pd.Series(
                                   [w.split(': ', 2)[-1][:60] for w in t[~t.eligible].why]).value_counts().head(8).items()})
    F = pd.read_csv(a.features)
    R['R5_phi1'] = dict(note='phi1 = -ln N_eff with N_eff over all valid minimised conformers, duplicates not merged: a rigid '
                             'state with 30 identical minima gets N_eff = 30 (phi1 = -3.40), the most "entropic" value',
                        n_states_phi1_at_minus_ln_nvalid=int((np.abs(F.phi1 + np.log(F.n_valid.clip(lower=1))) < 1e-3).sum()),
                        n_states_with_phi=int(F.phi1.notna().sum()))
    json.dump(R, open(os.path.join(a.out, 'e2_robust.json'), 'w'), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float))


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('features')
    s.add_argument('--data', required=True); s.add_argument('--out', required=True); s.add_argument('--procs', type=int, default=8)
    s = sp.add_parser('fit')
    for k in ('data', 'features', 'energies', 'out'):
        s.add_argument(f'--{k}', required=True)
    s = sp.add_parser('robust')
    for k in ('data', 'features', 'energies', 'out'):
        s.add_argument(f'--{k}', required=True)
    a = ap.parse_args()
    dict(features=cmd_features, fit=cmd_fit, robust=cmd_robust)[a.cmd](a)


if __name__ == '__main__':
    main()
