"""Checkpoint inventory + training-distribution check for charged states (read-only on the frozen TD checkout).

For every model dir under --workdir: model_parameters.yml essentials and SHA-256 of the checkpoint.
For the QM9 training split (same file list + split.npy indexing as TD's ConformerDataset): how many training molecules
carry any formal charge, and which charged motifs occur. This decides whether ionised states are in-distribution.
"""
import argparse
import glob
import json
import os
import os.path as osp
from collections import Counter

import numpy as np
import yaml
from rdkit import Chem, RDLogger

from jmcg.condgen.td_backend import sha256

RDLogger.DisableLog('rdApp.*')


def charge_stats(smiles):
    n_ok = n_charged = n_net = 0
    motifs, examples = Counter(), []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        n_ok += 1
        q = [a for a in m.GetAtoms() if a.GetFormalCharge()]
        if q:
            n_charged += 1
            n_net += Chem.GetFormalCharge(m) != 0
            for a in q:
                motifs[f'{a.GetSymbol()}{a.GetFormalCharge():+d}(H{a.GetTotalNumHs()})'] += 1
            if len(examples) < 15:
                examples.append(s)
    return dict(n_parsed=n_ok, n_with_formal_charge=n_charged, frac_with_formal_charge=n_charged / max(n_ok, 1),
                n_net_charged=n_net, charged_atom_motifs=dict(motifs.most_common(20)), examples=examples)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--qm9_raw', required=True, help='dir of GEOM-QM9 *.pickle (TD QM9_RAW)')
    ap.add_argument('--qm9_split', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    res = dict(checkpoints=[])
    for yml in sorted(glob.glob(osp.join(a.workdir, '**', 'model_parameters.yml'), recursive=True)):
        d = osp.dirname(yml)
        with open(yml) as f:
            p = yaml.full_load(f)
        ck = osp.join(d, 'best_model.pt')
        res['checkpoints'].append(dict(dir=d, dataset=p.get('dataset'), in_node_features=p.get('in_node_features'),
                                       second_order=p.get('use_second_order_repr'), ns=p.get('ns'), nv=p.get('nv'),
                                       best_model_sha256=sha256(ck) if osp.exists(ck) else None))
    files = np.array(sorted(glob.glob(osp.join(a.qm9_raw, '*.pickle'))))
    split = np.load(a.qm9_split, allow_pickle=True)
    root = a.qm9_raw if a.qm9_raw.endswith('/') else a.qm9_raw + '/'
    for name, idx in zip(['train', 'val', 'test'], split):
        smi = [f[len(root):-7] for f in files[sorted(idx)]]
        res[f'qm9_{name}'] = charge_stats(smi)
        print(name, {k: v for k, v in res[f'qm9_{name}'].items() if k != 'examples'}, flush=True)
    os.makedirs(osp.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res['checkpoints'], indent=1))


if __name__ == '__main__':
    main()
