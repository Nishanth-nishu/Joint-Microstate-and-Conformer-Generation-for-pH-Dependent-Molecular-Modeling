"""N1 pilot: 2D vs 3D-conformer input / aggregation ablation of released Uni-pKa (see research/experiments/N1_pilot/PREREG.md).

One forward pass per microstate (10 ETKDG+MMFF conformers + 1 2D copy, unipkainfer defaults); every arm is computed from
the per-conformer free energies, so arms differ only in which coordinates/aggregation enter the macro-pKa.
"""
import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from unipkainfer.pka_predictor import free_energy as fe
from unipkainfer.pka_predictor.free_energy import UnipkaFreeEnergyConfig, predict_standard_free_energies

from jmcg.data.build_sensitive_set import flags, flags_v2

SETS = ['sampl6', 'sampl7', 'sampl8', 'novartis_acid', 'novartis_base']
N3D = fe.UNIPKA_CONF_SIZE - 1          # conformer indices 0..9 are 3D, index 10 is the 2D copy
_CONF = []


def _capture(conformer_results, n_molecules):
    _CONF.append(conformer_results.copy())
    return _orig_agg(conformer_results, n_molecules)


_orig_agg = fe._aggregate_free_energy_predictions
fe._aggregate_free_energy_predictions = _capture


def lse(x):
    m = max(x)
    return m + math.log(sum(math.exp(v - m) for v in x))


def read_set(path):
    df = pd.read_csv(path, sep='\t', index_col=0)
    rows = []
    for i, r in df.iterrows():
        a, b = r['SMILES'].split('>>')
        rows.append({'idx': i, 'acid': a.split(','), 'base': b.split(','), 'target': float(r['TARGET'])})
    return rows


def canon(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m is not None else None


def state_energies(conf, how):
    g = conf.sort_values('conformer_index')['conformer_free_energy'].to_numpy()
    if how == 'all_mean':
        return float(g.mean())
    if how == '2d':
        return float(g[N3D])
    g3 = g[:N3D]
    if how == '3d_mean':
        return float(g3.mean())
    if how == '3d_lse':
        return -lse([-v for v in g3]) + math.log(len(g3))
    if how == '3d_min':
        return float(g3.min())
    if how == '3d_single':
        return float(g3[0])
    raise ValueError(how)


ARMS = {'A3D-mean': 'all_mean', 'A2D': '2d', 'A3D-only-mean': '3d_mean', 'A3D-lse': '3d_lse',
        'A3D-min': '3d_min', 'A3D-single': '3d_single'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True, help='ext/unipka_orig/dataset')
    ap.add_argument('--model_dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--limit', type=int, default=0, help='smoke test: transitions per set')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    transitions = []
    for s in SETS:
        rows = read_set(os.path.join(a.data, f'{s}.tsv'))
        if a.limit:
            rows = rows[:a.limit]
        for r in rows:
            r['set'] = s
        transitions += rows
    smiles = sorted({c for t in transitions for c in map(canon, t['acid'] + t['base']) if c})
    print(f'{len(transitions)} transitions, {len(smiles)} unique microstates', flush=True)

    per_fold = {}
    for fold in range(5):
        t0 = time.time()
        _CONF.clear()
        cfg = UnipkaFreeEnergyConfig(model_dir=a.model_dir, folds=(fold,), nthreads=a.threads, gpu=False)
        predict_standard_free_energies([Chem.MolFromSmiles(s) for s in smiles], task_name=f'n1_fold{fold}', config=cfg)
        conf = pd.concat(_CONF)
        conf['smiles_c'] = conf['smiles'].map(canon)
        conf.to_csv(os.path.join(a.out, f'conformer_free_energies_fold{fold}.csv'), index=False)
        per_fold[fold] = {smi: grp for smi, grp in conf.groupby('smiles_c')}
        print(f'fold {fold}: {len(conf)} conformer predictions in {time.time() - t0:.0f}s', flush=True)

    uncharger = rdMolStandardize.Uncharger()
    out_rows = []
    for t in transitions:
        acid = [canon(s) for s in t['acid']]
        base = [canon(s) for s in t['base']]
        row = {'set': t['set'], 'idx': t['idx'], 'target': t['target'], 'n_acid': len(acid), 'n_base': len(base)}
        try:
            for arm, how in ARMS.items():
                vals = []
                for fold in range(5):
                    ga = [state_energies(per_fold[fold][s], how) for s in acid]
                    gb = [state_energies(per_fold[fold][s], how) for s in base]
                    vals.append((lse([-g for g in ga]) - lse([-g for g in gb])) / math.log(10))
                row[arm] = float(np.mean(vals))
            # per-conformer pKa spread (fold-averaged), conformer c used for every microstate
            spreads = []
            for fold in range(5):
                pk = []
                for c in range(N3D):
                    ga = [per_fold[fold][s].sort_values('conformer_index')['conformer_free_energy'].iloc[c] for s in acid]
                    gb = [per_fold[fold][s].sort_values('conformer_index')['conformer_free_energy'].iloc[c] for s in base]
                    pk.append((lse([-g for g in ga]) - lse([-g for g in gb])) / math.log(10))
                spreads.append(float(np.std(pk)))
            row['conf_pka_std'] = float(np.mean(spreads))
        except KeyError as e:
            row['error'] = f'missing microstate {e}'
            out_rows.append(row)
            continue
        parent = uncharger.uncharge(Chem.MolFromSmiles(acid[0]))
        heavy = Chem.RemoveHs(parent)
        row.update(flags(heavy)); row.update(flags_v2(heavy))
        row['parent_smiles'] = Chem.MolToSmiles(parent)
        out_rows.append(row)

    with open(os.path.join(a.out, 'transitions.csv'), 'w', newline='') as fh:
        keys = sorted({k for r in out_rows for k in r})
        w = csv.DictWriter(fh, keys)
        w.writeheader(); w.writerows(out_rows)
    print(f'wrote {len(out_rows)} transitions; failures {sum("error" in r for r in out_rows)}', flush=True)


if __name__ == '__main__':
    sys.exit(main())
