"""Per-molecule torsion strata of a test-format csv: m = number of TD-rotatable bonds (TD skips the model when
``edge_mask.sum() == 0``, ``generate_confs.py:150-153``) and n_rot_heavy (the rigid stratum, DECISION fix 4).

One implementation of both counts (``flexitors.geometry.torsion_types``), checked against ``tools/breakdown.py``'s
``n_rot_heavy`` column on gnode118 (review P1 S3). Used by the re-score's m = 0 / m > 0 split (review P2 S8a).
"""
from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

from flexitors.geometry.torsion_types import heavy_rotatable_bonds, is_structural_stereo, td_rotatable_bonds

RDLogger.DisableLog('rdApp.*')


def _orders(mh: Chem.Mol) -> list[Chem.BondType]:
    return [mh.GetBondWithIdx(k).GetBondType() for k in td_rotatable_bonds(mh)]


def strata(csv: pd.DataFrame) -> pd.DataFrame:
    """Per row: m (TD-rotatable bonds), n_rot_heavy, and the TD-rotatable double / triple / structural-stereo counts
    that define the Exp 3 affected subsets (EXP-3 §9.1: 102 QM9 molecules with a double bond, 211 with a double or
    triple bond for K1-canonical). Unparsable rows get NaN."""
    rows = []
    for r in csv.itertuples(index=False):
        m = Chem.MolFromSmiles(r.corrected_smiles)
        row = dict(smiles=r.smiles, corrected_smiles=r.corrected_smiles)
        if m is None or '.' in r.corrected_smiles:
            rows.append({**row, 'm': float('nan'), 'n_rot_heavy': float('nan')})
            continue
        mh = Chem.AddHs(m)
        orders = _orders(mh)
        rows.append({**row, 'm': len(orders), 'n_rot_heavy': len(heavy_rotatable_bonds(mh)),
                     'n_double': orders.count(Chem.BondType.DOUBLE), 'n_triple': orders.count(Chem.BondType.TRIPLE),
                     'n_structural': sum(is_structural_stereo(mh.GetBondWithIdx(k)) for k in td_rotatable_bonds(mh))})
    return pd.DataFrame(rows)


def main() -> None:
    ap = ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--keys_dir', default=None, help='also write the Exp 3 affected-subset key files here')
    ap.add_argument('--keys_prefix', default='', help="e.g. 'drugs_' for keys_drugs_double.txt")
    a = ap.parse_args()
    S = strata(pd.read_csv(a.csv))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    S.to_csv(a.out, index=False)
    if a.keys_dir:
        d = Path(a.keys_dir)
        d.mkdir(parents=True, exist_ok=True)
        both = (S.n_double > 0) | (S.n_triple > 0)
        for name, sel in (('double', S.n_double > 0), ('double_or_triple', both), ('unaffected', ~both & (S.m > 0))):
            (d / f'keys_{a.keys_prefix}{name}.txt').write_text(''.join(f'{k}\n' for k in S.corrected_smiles[sel]))
    print(f'{len(S)} rows; m = 0: {int((S.m == 0).sum())}; n_rot_heavy = 0: {int((S.n_rot_heavy == 0).sum())}')


if __name__ == '__main__':
    main()
