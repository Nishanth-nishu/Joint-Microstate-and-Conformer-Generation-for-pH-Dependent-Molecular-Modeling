"""Tautomer pilot stage 0: zero-shot tautomer ratios from released Uni-pKa microstate free energies on Tautobase.

log10 K = log10([T2]/[T1]) = (G_T1 - G_T2) / ln 10, with G the per-microstate free energy (kT units) that Uni-pKa uses in
its macro-pKa log-sum-exp. Same frozen weights, conformer protocol and folds as the N1 pilot (importing n1_pilot applies
its patches: faithful conformer seed loop, per-conformer capture).
"""
import argparse
import csv
import math
import os
import sys
import time

import numpy as np
import pandas as pd
from rdkit import Chem

from jmcg.experiments import n1_pilot as n1
from jmcg.data.build_sensitive_set import MEASURED, build_tautobase, unipka_skeleton_keys
from unipkainfer.pka_predictor.free_energy import UnipkaFreeEnergyConfig, predict_standard_free_energies

ARMS = {'A3D-mean': 'all_mean', 'A2D': '2d', 'A3D-only-mean': '3d_mean'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tautobase', required=True, help='ext/tautobase/Tautobase_SMIRKS.txt')
    ap.add_argument('--dwar', required=True, help='ext/unipka_orig/dataset/dwar-iBond.tsv (leakage flags)')
    ap.add_argument('--model_dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--limit', type=int, default=0, help='smoke test: number of pairs')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    rows, keep = build_tautobase(a.tautobase, kinds=MEASURED + ('e', 'c'))
    print('tautobase filter', dict(keep), flush=True)
    if a.limit:
        rows = rows[:a.limit]
    dwar = unipka_skeleton_keys(a.dwar)
    smiles = sorted({n1.canon(r[k]) for r in rows for k in ('smiles_t1', 'smiles_t2')})
    print(f'{len(rows)} pairs, {len(smiles)} unique tautomers', flush=True)

    per_fold = {}
    for fold in range(5):
        t0 = time.time()
        n1._CONF.clear()
        cfg = UnipkaFreeEnergyConfig(model_dir=a.model_dir, folds=(fold,), nthreads=a.threads, gpu=False)
        predict_standard_free_energies([Chem.MolFromSmiles(s) for s in smiles], task_name=f'taut_fold{fold}', config=cfg)
        conf = pd.concat(n1._CONF)
        conf['smiles_c'] = conf['smiles'].map(n1.canon)
        conf.to_csv(os.path.join(a.out, f'conformer_free_energies_fold{fold}.csv'), index=False)
        per_fold[fold] = {smi: grp for smi, grp in conf.groupby('smiles_c')}
        print(f'fold {fold}: {len(conf)} conformer predictions in {time.time() - t0:.0f}s', flush=True)

    out = []
    for r in rows:
        t1, t2 = n1.canon(r['smiles_t1']), n1.canon(r['smiles_t2'])
        row = dict(r)
        row['in_unipka_finetune'] = int(Chem.MolToInchiKey(Chem.MolFromSmiles(t1))[:14] in dwar
                                        or Chem.MolToInchiKey(Chem.MolFromSmiles(t2))[:14] in dwar)
        try:
            for arm, how in ARMS.items():
                vals = [(n1.state_energies(per_fold[f][t1], how) - n1.state_energies(per_fold[f][t2], how)) / math.log(10)
                        for f in range(5)]
                row[f'logk_{arm}'] = float(np.mean(vals))
        except KeyError as e:
            row['error'] = f'missing tautomer {e}'
        out.append(row)

    with open(os.path.join(a.out, 'pairs.csv'), 'w', newline='') as fh:
        keys = sorted({k for r in out for k in r})
        w = csv.DictWriter(fh, keys)
        w.writeheader(); w.writerows(out)
    print(f'wrote {len(out)} pairs; failures {sum("error" in r for r in out)}', flush=True)


if __name__ == '__main__':
    sys.exit(main())
