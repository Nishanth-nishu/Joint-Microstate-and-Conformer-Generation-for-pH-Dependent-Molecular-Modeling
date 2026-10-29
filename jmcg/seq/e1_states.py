"""E1 (state component of the sequential baseline), on the tasks that have experimental reference labels.

Run (jmcg venv): python -m jmcg.seq.e1_states --data data_manifests/canonical_v1 --model_dir <unipka> --out <dir>
                    [--n1_transitions results/n1_pilot/full_v2/transitions.csv] [--limit N] [--threads 8]
Tasks (all labels experimental, from canonical_v1; nothing is fitted here, so no split is consumed):
  external macro-pKa  356 transitions (SAMPL6/7/8, Novartis acid/base), role external_test.
     oracle arm      : acid/base microstate sets as given by the dataset (= N1 A3D-mean protocol).
     enumerated arm  : parent = first acid microstate; ensemble from jmcg.seq.enumerate; the predicted macro-pKa
                       between the charge levels (q_acid, q_acid - 1). Missing level -> failure, counted, never dropped.
     state recall    : fraction of the dataset's acid and base microstates present in the enumerated ensemble.
  tautomer log K      Tautobase (curated_eval), measured / estimated / calculated reported separately.
     oracle arm      : log10 K = (G_T1 - G_T2)/ln10 from the two given tautomers (= T-pilot R0, 5-fold mean G).
     enumerator recall: is T2 in the ensemble enumerated from T1 (same charge level)?
Uncertainty: percentile bootstrap (10,000) over split components (obs.split_group), so related molecules resample
together. Reproduction check: oracle external pKa vs the N1 run (same weights, protocol) per transition.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from rdkit import Chem

from jmcg.seq import enumerate as en
from jmcg.seq import score as sc


def boot(err, groups, stat, n=10000, seed=0):
    df = pd.DataFrame(dict(e=err, g=groups))
    gs = [g.e.to_numpy() for _, g in df.groupby('g')]
    rng = np.random.default_rng(seed)
    vals = [stat(np.concatenate([gs[i] for i in rng.integers(0, len(gs), len(gs))])) for _ in range(n)]
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


def metrics(pred, target, groups):
    pred, target = np.asarray(pred, float), np.asarray(target, float)
    e = pred - target
    mae = lambda x: float(np.mean(np.abs(x)))  # noqa: E731
    rmse = lambda x: float(np.sqrt(np.mean(x ** 2)))  # noqa: E731
    slope = float(np.polyfit(target, pred, 1)[0]) if len(e) > 2 else float('nan')
    return dict(n=int(len(e)), n_groups=int(pd.Series(groups).nunique()), MAE=mae(e), MAE_ci95=boot(e, groups, mae),
                RMSE=rmse(e), RMSE_ci95=boot(e, groups, rmse), mean_signed=float(e.mean()), slope_pred_on_target=slope)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--model_dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--n1_transitions', default='')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--templates', default='simple_smarts_pattern.tsv,smarts_pattern.tsv')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    obs = pd.read_csv(os.path.join(a.data, 'observations.csv'), keep_default_na=False)
    smi = pd.read_csv(os.path.join(a.data, 'states.csv'), keep_default_na=False).set_index('state_id').smiles
    ext = obs[obs.role == 'external_test'].copy()
    tau = obs[obs.task == 'tautomer_logk'].copy()
    if a.limit:
        ext, tau = ext.head(a.limit), tau.head(a.limit)
    ext['A'] = ext.states_A.str.split(';')
    ext['B'] = ext.states_B.str.split(';')

    # enumeration: one ensemble per (template, seed). External seeds are scored; tautomer seeds are only used for
    # enumerator recall (no scoring needed).
    templates = a.templates.split(',')
    ext_seeds = sorted(set(smi[x[0]] for x in ext.A))
    tau_seeds = sorted(set(smi[s] for s in tau.states_A))
    seeds = sorted(set(ext_seeds) | set(tau_seeds))
    ens_rows, fail_rows = [], []
    for t in templates:
        for s in seeds:
            try:
                states, f = en.enumerate_ensemble(s, template=t)
            except Exception as e:  # noqa: BLE001
                fail_rows.append(dict(template=t, seed=s, smiles=s, level='', reason=f'{type(e).__name__}: {e}'[:300]))
                continue
            ens_rows += [dict(template=t, seed=s, **x) for x in states]
            fail_rows += [dict(template=t, seed=s, **x) for x in f]
    ens = pd.DataFrame(ens_rows)
    ens.to_csv(os.path.join(a.out, 'ensembles.csv'), index=False)
    pd.DataFrame(fail_rows, columns=['template', 'seed', 'smiles', 'level', 'reason']).to_csv(
        os.path.join(a.out, 'enumeration_failures.csv'), index=False)
    oracle = {smi[x] for col in (ext.A, ext.B) for xs in col for x in xs} | set(smi[tau.states_A]) | set(smi[tau.states_B])
    allsmi = sorted(set(ens[ens.seed.isin(ext_seeds)].smiles) | oracle)
    print(f'{len(seeds)} enumeration seeds, {len(ens)} enumerated states, {len(allsmi)} states to score', flush=True)
    energies = sc.free_energies(allsmi, a.model_dir, threads=a.threads, task='e1')
    energies.to_csv(os.path.join(a.out, 'energies.csv'), index=False)
    E = energies.pivot_table(index='smiles', columns='fold', values='G')
    folds = list(E.columns)

    ens_by_seed = {k: g for k, g in ens.groupby(['template', 'seed'])}
    # external macro-pKa (one row per transition x template; the oracle value repeats)
    rows = []
    for _, o in ext.iterrows():
        A, B = [smi[x] for x in o.A], [smi[x] for x in o.B]
        orc = [(sc.lse(-E.loc[A, f].to_numpy()) - sc.lse(-E.loc[B, f].to_numpy())) / sc.LN10 + sc.PKA_MEAN for f in folds]
        r = dict(obs_id=o.obs_id, source=o.source, source_record=o.source_record, split_group=o.split_group,
                 target=float(o.value), oracle=float(np.mean(orc)), n_acid=len(A), n_base=len(B))
        qa = Chem.GetFormalCharge(Chem.MolFromSmiles(A[0]))
        for t in templates:
            rt = dict(r, template=t)
            g = ens_by_seed.get((t, A[0]))
            if g is None or not ((g.charge == qa).any() and (g.charge == qa - 1).any()):
                rt['enum_failure'] = 'no ensemble' if g is None else f'levels present {sorted(set(g.charge))}, need {qa},{qa - 1}'
            else:
                G = E.loc[g.smiles].to_numpy()
                rt['enumerated'] = float(np.mean([sc.macro_pkas(G[:, k], g.charge.to_numpy())[(qa, qa - 1)]
                                                  for k in range(len(folds))]))
                have = set(g.smiles)
                rt['recall_acid'] = float(np.mean([x in have for x in A]))
                rt['recall_base'] = float(np.mean([x in have for x in B]))
                rt['n_ens_acid_level'] = int((g.charge == qa).sum())
                rt['n_ens_base_level'] = int((g.charge == qa - 1).sum())
                rt['ensemble_size'] = len(g)
            rows.append(rt)
    X = pd.DataFrame(rows)
    X.to_csv(os.path.join(a.out, 'e1_external.csv'), index=False)

    # tautomers
    trows = []
    for _, o in tau.iterrows():
        s1, s2 = smi[o.states_A], smi[o.states_B]
        k = [sc.log10_k(E.loc[s1, f], E.loc[s2, f]) for f in folds]
        rec = dict(obs_id=o.obs_id, label_status=o.label_status, split=o.split, split_group=o.split_group,
                   target=float(o.value), oracle=float(np.mean(k)), oracle_sd_folds=float(np.std(k)))
        for t in templates:
            g = ens_by_seed.get((t, s1))
            have = set(g.smiles) if g is not None else set()
            rec[f't1_in_ensemble|{t}'] = s1 in have
            rec[f't2_in_ensemble|{t}'] = s2 in have
        trows.append(rec)
    T = pd.DataFrame(trows)
    T.to_csv(os.path.join(a.out, 'e1_tautomer.csv'), index=False)

    S = dict(n_seeds=len(seeds), n_enumerated_states=len(ens), n_enumeration_failures=len(fail_rows),
             enumerator=en.provenance(), PKA_MEAN=sc.PKA_MEAN)
    X1 = X[X.template == templates[0]]
    S['external_oracle'] = metrics(X1.oracle, X1.target, X1.split_group)
    if 'enumerated' not in X:
        X['enumerated'] = np.nan
    S['external_enumerated'] = {}
    for t, Xt in X.groupby('template'):
        ok = Xt[Xt.enumerated.notna()]
        S['external_enumerated'][t] = dict(
            n_failed=int(len(Xt) - len(ok)),
            failures=Xt[Xt.enumerated.isna()].reindex(columns=['source_record', 'enum_failure']).to_dict('records')[:50],
            on_successes=metrics(ok.enumerated, ok.target, ok.split_group) if len(ok) > 2 else None,
            oracle_on_same_successes=metrics(ok.oracle, ok.target, ok.split_group) if len(ok) > 2 else None,
            mean_recall_acid=float(ok.recall_acid.mean()) if len(ok) else None,
            mean_recall_base=float(ok.recall_base.mean()) if len(ok) else None,
            full_recall_both=int(((ok.recall_acid == 1) & (ok.recall_base == 1)).sum()) if len(ok) else 0)
    S['external_by_source'] = {s: metrics(g.oracle, g.target, g.split_group) for s, g in X1.groupby('source') if len(g) > 2}
    S['tautomer_oracle'] = {l: metrics(g.oracle, g.target, g.split_group) for l, g in T.groupby('label_status') if len(g) > 2}
    S['tautomer_predict_zero_MAE'] = {l: float(np.abs(g.target).mean()) for l, g in T.groupby('label_status')}
    S['tautomer_enumerator_recall'] = {
        f'{l}|{t}': dict(t1_in_own_ensemble=f"{int(g[f't1_in_ensemble|{t}'].sum())}/{len(g)}",
                         t2_in_ensemble_of_t1=f"{int(g[f't2_in_ensemble|{t}'].sum())}/{len(g)}")
        for l, g in T.groupby('label_status') for t in templates}
    if a.n1_transitions and os.path.exists(a.n1_transitions):
        n1 = pd.read_csv(a.n1_transitions)
        n1['source_record'] = [f'{s}:{i}' for s, i in zip(n1['set'], n1['idx'])]
        # N1's transitions.csv numbers two sampl8 records idx 17 (its record 18 is labelled 17), so (record, target)
        # is the key; a row left unmatched falls back to (set, target) when that is unique in N1.
        n1['tk'] = n1.target.round(3)
        x = X1.assign(tk=X1.target.round(3), set=X1.source_record.str.split(':').str[0])
        m = x.merge(n1[['source_record', 'tk', 'A3D-mean']], on=['source_record', 'tk'], how='inner')
        rest = x[~x.obs_id.isin(m.obs_id)]
        uniq = n1.drop_duplicates(['set', 'tk'], keep=False)[['set', 'tk', 'A3D-mean']]
        m2 = rest.merge(uniq, on=['set', 'tk'], how='inner')
        m = pd.concat([m, m2])
        S['reproduction_vs_N1'] = dict(n_matched=len(m), n_matched_by_set_target=len(m2), n_external=len(X1),
                                       max_abs_diff=float((m.oracle - m['A3D-mean']).abs().max()))
    json.dump(S, open(os.path.join(a.out, 'e1_summary.json'), 'w'), indent=1, default=str)
    print(json.dumps({k: v for k, v in S.items() if k not in ('enumerator',)}, indent=1, default=str)[:8000])


if __name__ == '__main__':
    main()
