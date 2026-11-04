"""EXP1-V2-T1 leakage check: no E1 test molecule (canonical SMILES, stereo removed) may be in the paired training or
validation cache. Run (TD venv): python experiments/exp1/t1_leak_check.py <subset.csv> <cache.train> <cache.val>"""
import csv
import pickle
import sys

from rdkit import Chem


def canon(m):
    m = Chem.RemoveHs(Chem.Mol(m), sanitize=False)
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)


def main():
    sub, *caches = sys.argv[1:]
    test = set()
    for r in csv.DictReader(open(sub)):
        m = Chem.MolFromSmiles(r['corrected_smiles']) or Chem.MolFromSmiles(r['smiles'])
        test.add(canon(m))
    for c in caches:
        pts = pickle.load(open(c, 'rb'))
        s = {canon(d.mol) for d in pts}
        print(c, 'molecules', len(pts), 'unique', len(s), 'test', len(test), 'overlap', len(s & test), sorted(s & test)[:10])


if __name__ == '__main__':
    main()
