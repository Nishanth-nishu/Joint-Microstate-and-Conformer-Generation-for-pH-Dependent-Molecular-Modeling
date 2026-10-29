"""E6: stratified E1 results (exploratory; strata with n < 20 are flagged as small-sample).

Run (jmcg venv): python -m jmcg.seq.e6_strata --data data_manifests/canonical_v1 --e1 <E1 dir> --out <dir>
External macro-pKa strata (per transition, from the dataset's microstates):
  acid_level_charge  +1 or more (base protonation) / 0 (neutral acid) / -1 or less (anion acid)
  zwitterion         any acid/base microstate carrying both a + and a - atom
  size               heavy atoms <= 9 / 10-20 / > 20 (of the first acid microstate)
  td_qm9_domain      every microstate within H,C,N,O,F and <= 9 heavy atoms
  stereo             any microstate with an unspecified stereo element / specified stereo / none
  in_unipka_train    family also in Dwar-iBonD (released Uni-pKa's fine-tuning set)
Tautomer strata (measured pairs): size, td_qm9_domain, label status, heteroatom of the moving H (N/O/S/C mixtures).
"""
import argparse
import json
import os
import re

import numpy as np
import pandas as pd
from rdkit import Chem

from jmcg.seq.e1_states import metrics


def zwit(smi):
    m = Chem.MolFromSmiles(smi)
    q = [a.GetFormalCharge() for a in m.GetAtoms()]
    return any(x > 0 for x in q) and any(x < 0 for x in q)


def strat_table(df, col, pred='oracle'):
    out = {}
    for k, g in df.groupby(col):
        if len(g) < 3:
            out[str(k)] = dict(n=len(g), note='too few')
            continue
        m = metrics(g[pred], g.target, g.split_group)
        m['small_sample'] = len(g) < 20
        out[str(k)] = {x: m[x] for x in ('n', 'n_groups', 'MAE', 'MAE_ci95', 'mean_signed', 'small_sample')}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--e1', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    obs = pd.read_csv(os.path.join(a.data, 'observations.csv'), keep_default_na=False).set_index('obs_id')
    st = pd.read_csv(os.path.join(a.data, 'states.csv'), keep_default_na=False).set_index('state_id')
    upk = set(obs[obs.role == 'unipka_train'].family_key)
    td_ok = lambda sid: set(re.findall('[A-Z][a-z]?', st.elements[sid])) <= set('HCNOF') and st.n_heavy[sid] <= 9  # noqa: E731

    X = pd.read_csv(os.path.join(a.e1, 'e1_external.csv'))
    X = X[X.template == X.template.iloc[0]].copy()
    A = [obs.states_A[o].split(';') for o in X.obs_id]
    B = [obs.states_B[o].split(';') for o in X.obs_id]
    X['acid_level_charge'] = [('+1 or more' if q >= 1 else '0' if q == 0 else '-1 or less')
                              for q in (st.formal_charge[x[0]] for x in A)]
    X['zwitterion'] = [any(zwit(st.smiles[s]) for s in x + y) for x, y in zip(A, B)]
    X['size'] = [('<=9' if n <= 9 else '10-20' if n <= 20 else '>20') for n in (st.n_heavy[x[0]] for x in A)]
    X['td_qm9_domain'] = [all(td_ok(s) for s in x + y) for x, y in zip(A, B)]
    X['stereo'] = [('unspecified' if any(st.n_unspecified_stereo[s] > 0 for s in x + y) else
                    'specified' if any('@' in st.smiles[s] or '/' in st.smiles[s] or '\\' in st.smiles[s] for s in x + y)
                    else 'none') for x, y in zip(A, B)]
    X['in_unipka_train'] = [obs.family_key[o] in upk for o in X.obs_id]
    R = {'external_oracle': {c: strat_table(X, c) for c in ('acid_level_charge', 'zwitterion', 'size', 'td_qm9_domain',
                                                            'stereo', 'in_unipka_train')}}
    T = pd.read_csv(os.path.join(a.e1, 'e1_tautomer.csv'))
    T['size'] = [('<=9' if st.n_heavy[obs.states_A[o]] <= 9 else '10-20' if st.n_heavy[obs.states_A[o]] <= 20 else '>20')
                 for o in T.obs_id]
    T['td_qm9_domain'] = [td_ok(obs.states_A[o]) and td_ok(obs.states_B[o]) for o in T.obs_id]
    R['tautomer_oracle_measured'] = {c: strat_table(T[T.label_status == 'measured'], c) for c in ('size', 'td_qm9_domain')}
    R['tautomer_oracle_by_label'] = strat_table(T, 'label_status')
    json.dump(R, open(os.path.join(a.out, 'e6_strata.json'), 'w'), indent=1, default=str)
    print(json.dumps(R, indent=1, default=str)[:10000])


if __name__ == '__main__':
    main()
