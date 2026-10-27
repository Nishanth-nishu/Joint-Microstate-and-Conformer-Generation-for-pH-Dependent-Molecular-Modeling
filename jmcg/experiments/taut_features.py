"""Tautomer pilot stage 1: cache Uni-pKa per-conformer CLS features (input of the free-energy head) for one fold.

Microstates = every SMILES in Dwar-iBonD, the 5 external pKa test sets and Tautobase (all label kinds). The head's own
output is recomputed from the cached features and checked against unipkainfer's conformer free energies, so the cache is
verified to be the exact head input in the exact row order.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
from rdkit import Chem

from jmcg.experiments import n1_pilot as n1
from jmcg.data.build_sensitive_set import MEASURED, build_tautobase
from unipkainfer.unimol.models import unimol as unimol_models
from unipkainfer.pka_predictor.free_energy import UnipkaFreeEnergyConfig, predict_standard_free_energies

_FEATS = []
_HEAD = {}
_orig_head_forward = unimol_models.ClassificationHead.forward


def _capture_head(self, features, **kwargs):
    _FEATS.append(features[:, 0, :].detach().float().cpu().numpy())
    _HEAD['module'] = self
    return _orig_head_forward(self, features, **kwargs)


unimol_models.ClassificationHead.forward = _capture_head


def all_smiles(data_dir, tautobase):
    smi = set()
    for name in ['dwar-iBond'] + n1.SETS:
        for r in pd.read_csv(os.path.join(data_dir, f'{name}.tsv'), sep='\t', index_col=0)['SMILES']:
            smi.update(r.replace('>>', ',').split(','))
    rows, _ = build_tautobase(tautobase, kinds=MEASURED + ('e', 'c'))
    for r in rows:
        smi.update([r['smiles_t1'], r['smiles_t2']])
    return sorted({c for c in map(n1.canon, smi) if c})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True, help='ext/unipka_orig/dataset')
    ap.add_argument('--tautobase', required=True)
    ap.add_argument('--model_dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--fold', type=int, default=0)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--limit', type=int, default=0, help='smoke test: number of microstates')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    smiles = all_smiles(a.data, a.tautobase)
    if a.limit:
        smiles = smiles[:a.limit]
    print(f'{len(smiles)} unique microstates', flush=True)
    t0 = time.time()
    cfg = UnipkaFreeEnergyConfig(model_dir=a.model_dir, folds=(a.fold,), nthreads=a.threads, gpu=False)
    predict_standard_free_energies([Chem.MolFromSmiles(s) for s in smiles], task_name=f'feat_fold{a.fold}', config=cfg)
    conf = pd.concat(n1._CONF).reset_index(drop=True)
    feats = np.concatenate(_FEATS)
    print(f'{len(conf)} conformer rows, {feats.shape} features in {time.time() - t0:.0f}s', flush=True)
    assert len(conf) == len(feats), 'feature/row count mismatch'

    head = _HEAD['module'].eval()
    with torch.no_grad():
        g = head(torch.from_numpy(feats)[:, None, :]).view(-1).numpy()
    err = float(np.abs(g - conf['conformer_free_energy'].to_numpy()).max())
    print(f'head recompute max |dG| = {err:.2e}', flush=True)
    assert err < 1e-3, 'cached features do not reproduce the head output'

    conf['smiles_c'] = conf['smiles'].map(n1.canon)
    np.savez_compressed(os.path.join(a.out, f'features_fold{a.fold}.npz'), feats=feats.astype(np.float32),
                        smiles_c=conf['smiles_c'].to_numpy(dtype=str), conformer_index=conf['conformer_index'].to_numpy(),
                        g_released=conf['conformer_free_energy'].to_numpy(dtype=np.float32))
    torch.save(head, os.path.join(a.out, f'head_released_fold{a.fold}.pt'))  # whole module: keeps the checkpoint's activation
    print('saved', flush=True)


if __name__ == '__main__':
    sys.exit(main())
