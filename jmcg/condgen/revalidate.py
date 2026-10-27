"""Re-run validation on already written SDFs (no regeneration) after a validator change.

Run: python -m jmcg.condgen.revalidate --run_dir <dir with <id>.json, <id>__<method>.sdf, summary.csv> --out <new dir>
Writes <out>/summary.csv (same columns; validation counts recomputed, generation fields kept) and <out>/changes.json
(every conformer whose verdict changed, with old and new check values). The original run directory is not touched.
"""
import argparse
import glob
import json
import os

import pandas as pd
from rdkit import Chem

from jmcg.condgen import validate as va


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run_dir', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    old = pd.read_csv(os.path.join(a.run_dir, 'summary.csv'))
    rows, changes = [], []
    for r in old.to_dict('records'):
        sid, method = r['state_id'], r['method']
        rec = json.load(open(os.path.join(a.run_dir, f'{sid}.json')))
        target = rec['state']['state_smiles_canonical']
        mols = list(Chem.SDMolSupplier(os.path.join(a.run_dir, f'{sid}__{method}.sdf'), removeHs=False))
        bounds = va.rdDistGeom.GetMoleculeBoundsMatrix(mols[0])
        new = [va.check_conformer(m, target, bounds) for m in mols]
        for k, (o, n) in enumerate(zip(rec[method]['validation'], new)):
            if o['valid'] != n['valid'] or o.get('clash') != n.get('clash'):
                changes.append(dict(state_id=sid, method=method, conf=k,
                                    old={x: o.get(x) for x in ('valid', 'clash', 'n_clash_pairs', 'min_nonlocal_dist')},
                                    new={x: n.get(x) for x in ('valid', 'clash', 'n_clash_pairs', 'n_hbond_contacts',
                                                               'clash_pairs')}))
        r.update(va.summarise(new))
        r['n_hbond_contacts'] = sum(n['n_hbond_contacts'] for n in new)
        r['clash_pairs_remaining'] = json.dumps([n['clash_pairs'] for n in new if n['clash_pairs']])
        rows.append(r)
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'summary.csv'), index=False)
    json.dump(changes, open(os.path.join(a.out, 'changes.json'), 'w'), indent=1)
    for c in changes:
        print(json.dumps(c))
    print(f'revalidated {len(rows)} (state, method) rows; {len(changes)} conformer verdicts changed')


if __name__ == '__main__':
    main()
