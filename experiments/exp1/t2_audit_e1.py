"""EXP1-V2-T2 task 1: independent audit of the saved EXP1-V2-T1 E1 outputs (read only).

For each seed and arm:
  inputs     subset.csv molecules; TD input pickle td{s}.pkl; arm output pickle (same molecules, same conformer count per
             molecule, same atoms in the same order, no duplicated conformers, no NaN coordinates)
  eval       the evaluate_confs.py results pickle: one RMSD matrix per molecule (n_true x n_generated); matrix shape vs
             the conformer counts; failures; NaN matrices
  metrics    MAT-R/P and COV-R/P recomputed from the RMSD matrices (a) over every molecule in the eval pickle and (b) over
             the paired set of the compare JSON, and checked against the compare JSON
  bootstrap  the paired percentile bootstrap re-implemented from scratch (molecule resampling, 5,000 draws) and checked
             against the JSON CIs (same generator order reproduces them exactly; a fresh generator within Monte-Carlo error)
Writes one JSON with every check and, per metric, which molecules contribute.

Run (TD venv): python experiments/exp1/t2_audit_e1.py --in_dir <n200> --t1_dir <T1_e1> --out <json>
"""
from __future__ import annotations

import csv
import json
import pickle
from argparse import ArgumentParser

import numpy as np

THRS = (0.25, 0.5, 0.75)
ARMS = ('Lproj', 'Lnoproj', 'mmff', 'xtb', 'Lprojfix', 'Lnoprojfix')


def metrics(m):
    r, p = np.min(m, axis=1), np.min(m, axis=0)
    out = dict(MAT_R=float(r.mean()), MAT_P=float(p.mean()))
    for t in THRS:
        out[f'COV_R@{t}'] = float(np.mean(r < t) * 100)
        out[f'COV_P@{t}'] = float(np.mean(p < t) * 100)
    return out


def load_eval(path):
    d = pickle.load(open(path, 'rb'))
    per, nan, shapes = {}, [], {}
    for (smi, csmi), r in d['results'].items():
        shapes[csmi] = (smi, tuple(r['rmsd'].shape))
        if np.isnan(r['rmsd']).any():
            nan.append(csmi)
            continue
        per[csmi] = metrics(r['rmsd'])
    return per, nan, shapes, d.get('num_failures')


def check_confs(td, arm):
    """Arm output vs TD input: molecules, counts, atom identity/order, NaN, duplicated conformers."""
    out = dict(missing=sorted(set(td) - set(arm)), extra=sorted(set(arm) - set(td)), count_mismatch=[], atom_mismatch=[],
               none=0, nan_coords=0, duplicated_within_mol=0, conformers=0)
    for smi in td:
        if smi not in arm:
            continue
        a, b = td[smi], arm[smi]
        if len(a) != len(b):
            out['count_mismatch'].append(smi)
        seen = []
        for ma, mb in zip(a, b):
            if mb is None:
                out['none'] += 1
                continue
            out['conformers'] += 1
            za = [x.GetAtomicNum() for x in ma.GetAtoms()]
            zb = [x.GetAtomicNum() for x in mb.GetAtoms()]
            if za != zb:
                out['atom_mismatch'].append(smi)
            x = mb.GetConformer().GetPositions()
            if np.isnan(x).any():
                out['nan_coords'] += 1
            if any(np.allclose(x, y, atol=1e-6) for y in seen):
                out['duplicated_within_mol'] += 1
            seen.append(x)
    out['atom_mismatch'] = sorted(set(out['atom_mismatch']))
    return out


def bootstrap(dlt, rng, n=5000):
    bs = [dlt[rng.integers(0, len(dlt), len(dlt))].mean() for _ in range(n)]
    return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


def main():
    ap = ArgumentParser()
    ap.add_argument('--in_dir', required=True)
    ap.add_argument('--t1_dir', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    subset = list(csv.DictReader(open(f'{a.in_dir}/subset.csv')))
    R = dict(subset_molecules=len(subset), subset_true_conformers=sum(int(r['n_conformers']) for r in subset), seeds={})
    for s in (0, 1):
        S = R['seeds'][s] = {}
        td = pickle.load(open(f'{a.in_dir}/td{s}.pkl', 'rb'))
        n_true = {r['corrected_smiles']: int(r['n_conformers']) for r in subset}
        S['td_input'] = dict(molecules=len(td), conformers=sum(len(v) for v in td.values()),
                             subset_molecules_without_td_output=sorted(set(n_true) - set(td)),
                             count_not_2x_ntrue=sorted(k for k, v in td.items() if k in n_true and len(v) != 2 * n_true[k]))
        S['arm_outputs'] = {}
        for arm in ARMS:
            try:
                S['arm_outputs'][arm] = check_confs(td, pickle.load(open(f'{a.t1_dir}/td{s}_{arm}.pkl', 'rb')))
            except FileNotFoundError as e:
                S['arm_outputs'][arm] = dict(error=str(e))
        evals = {'td': f'{a.in_dir}/td{s}.eval.pkl', 'etkdg': f'{a.in_dir}/etkdg{s}.eval.pkl',
                 **{arm: f'{a.t1_dir}/td{s}_{arm}.eval.pkl' for arm in ARMS}}
        E = {k: load_eval(p) for k, p in evals.items()}
        S['eval'] = {k: dict(molecules_scored=len(v[0]) + len(v[1]), nan_dropped=v[1], num_failures=v[3],
                             all_scored_means={c: float(np.mean([r[c] for r in v[0].values()])) for c in next(iter(v[0].values()))})
                     for k, v in E.items()}
        # RMSD matrix columns must equal the arm's conformer count (generated side) for every molecule
        for k in ('td', *ARMS):
            confs = td if k == 'td' else pickle.load(open(f'{a.t1_dir}/td{s}_{k}.pkl', 'rb'))
            bad = [c for c, (smi, sh) in E[k][2].items() if smi in confs and sh[1] != len(confs[smi])]
            S['eval'][k]['rmsd_cols_ne_conformers'] = bad
        # compare JSONs: paired set, arm means and CIs
        S['compare'] = {}
        for name, ref, arms in (('compare', 'td', ('Lproj', 'Lnoproj', 'mmff', 'etkdg')),
                                ('ablation', 'Lnoproj', ('Lproj', 'mmff', 'td', 'etkdg')),
                                ('xtb_vs_td', 'td', ('Lprojfix', 'Lnoprojfix', 'xtb', 'mmff', 'etkdg')),
                                ('vs_xtb', 'xtb', ('Lprojfix', 'Lnoprojfix', 'Lproj', 'td')),
                                ('ablation_fix', 'Lnoprojfix', ('Lprojfix',))):
            J = json.load(open(f'{a.t1_dir}/{name}_seed{s}.json'))
            common = sorted(set.intersection(*[set(E[k][0]) for k in (ref, *arms)]))
            chk = dict(n_common_json=J['n_common'], n_common_recomputed=len(common), max_abs_mean_diff=0.0,
                       max_abs_ci_diff_same_rng=0.0, max_abs_ci_diff_fresh_rng=0.0)
            rng_same, rng_fresh = np.random.default_rng(0), np.random.default_rng(12345)
            jarms = {p.split('/')[-1].replace('.eval.pkl', ''): v for p, v in J['paired_arm_minus_ref'].items()}
            for k in arms:
                key = f'td{s}_{k}' if k not in ('td', 'etkdg') else f'{k}{s}'
                for c in J['paired_arm_minus_ref'][next(iter(J['paired_arm_minus_ref']))]:
                    dlt = np.array([E[k][0][m][c] - E[ref][0][m][c] for m in common])
                    jm = jarms[key][c]
                    chk['max_abs_mean_diff'] = max(chk['max_abs_mean_diff'], abs(dlt.mean() - jm['mean']))
                    ci = bootstrap(dlt, rng_same)
                    chk['max_abs_ci_diff_same_rng'] = max(chk['max_abs_ci_diff_same_rng'],
                                                          max(abs(x - y) for x, y in zip(ci, jm['ci95'])))
                    ci = bootstrap(dlt, rng_fresh)
                    chk['max_abs_ci_diff_fresh_rng'] = max(chk['max_abs_ci_diff_fresh_rng'],
                                                           max(abs(x - y) for x, y in zip(ci, jm['ci95'])))
            S['compare'][name] = chk
        # who contributes: paired molecules of the main compare, and who is excluded and why
        allm = set(E['td'][2])
        common = sorted(set.intersection(*[set(E[k][0]) for k in ('td', 'Lproj', 'Lnoproj', 'mmff', 'etkdg')]))
        S['contributors'] = dict(
            paired_molecules=common,
            excluded_no_td_output=sorted(set(n_true) - set(td)),
            excluded_nan_rmsd={k: E[k][1] for k in E},
            excluded_scored_by_td_but_missing_in_other_arm=sorted(allm - set(common) - set().union(*[set(E[k][1]) for k in E])),
            rule='each metric: per-molecule value from that molecule\'s full RMSD matrix (all its reference conformers x all '
                 'its generated conformers), then the unweighted mean over the paired molecules; CIs resample molecules')
    json.dump(R, open(a.out, 'w'), indent=1, default=str)
    for s in (0, 1):
        S = R['seeds'][s]
        print('seed', s, 'td_input', {k: (v if not isinstance(v, list) else len(v)) for k, v in S['td_input'].items()})
        for k, v in S['arm_outputs'].items():
            print('  out', k, {kk: (vv if not isinstance(vv, list) else len(vv)) for kk, vv in v.items()})
        for k, v in S['eval'].items():
            print('  eval', k, 'scored', v['molecules_scored'], 'nan', len(v['nan_dropped']), 'fail', v['num_failures'],
                  'cols_mismatch', len(v.get('rmsd_cols_ne_conformers', [])), 'MAT_R', round(v['all_scored_means']['MAT_R'], 4),
                  'COV_R@0.25', round(v['all_scored_means']['COV_R@0.25'], 2))
        for k, v in S['compare'].items():
            print('  cmp', k, v)
        print('  paired', len(S['contributors']['paired_molecules']), 'excluded', {k: v for k, v in S['contributors'].items()
                                                                                  if k.startswith('excluded')})


if __name__ == '__main__':
    main()
