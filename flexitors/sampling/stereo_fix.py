"""F-flip post-hoc arm (EXP-3 §3.4, §5; E-CHEM-046): on O's output, every specified structural double bond whose
geometric label differs from the SMILES target is rotated to the target and re-planarised (δ_b → 0), i.e. its
dihedral on TD's quadruple is set to φ_b. A rigid side rotation changes no bond length or angle (E-JOINT-026); no NFE.
F-filter, which needs fresh samples, is the ``F_filter`` arm of ``flexitors.sampling.generate``.
"""
from __future__ import annotations

import json
import pickle
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from flexitors.data.ez_targets import load_targets
from flexitors.geometry.internal import dihedral
from flexitors.provenance import write_provenance
from flexitors.sampling.stereo_guard import constrainable_bonds


def f_flip(mols: list[Chem.Mol], targets: dict[int, float]) -> tuple[list[Chem.Mol], int]:
    """Copies of ``mols`` with every wrong specified bond set to its target; returns (mols, number of bonds flipped)."""
    if not targets or not mols:
        return mols, 0
    spec = [b for b in constrainable_bonds(mols[0], targets)[0] if b.kind == 'specified']
    out, n_flip = [], 0
    for m in mols:
        m = Chem.Mol(m)
        conf = m.GetConformer()
        tau = dihedral(conf.GetPositions(), np.array([b.quad for b in spec]))
        for b, t in zip(spec, tau):
            if np.cos(t - b.phi) < 0:
                rdMolTransforms.SetDihedralRad(conf, *map(int, b.quad), float(b.phi))
                n_flip += 1
        out.append(m)
    return out, n_flip


def main() -> None:
    ap = ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--confs', required=True, help="O's confs.pkl")
    ap.add_argument('--ez_table', required=True)
    ap.add_argument('--out_dir', required=True)
    a = ap.parse_args()
    with open(a.confs, 'rb') as f:
        confs = pickle.load(f)
    targets = load_targets(a.ez_table)
    out, flips = {}, {}
    for smi, mols in confs.items():
        out[smi], n = f_flip(mols, targets.get(smi, {}))
        if n:
            flips[smi] = n
    d = Path(a.out_dir)
    d.mkdir(parents=True, exist_ok=True)
    with open(d / 'confs.pkl', 'wb') as f:
        pickle.dump(out, f)
    (d / 'budget.json').write_text(json.dumps(dict(nfe_extra=0, bonds_flipped=sum(flips.values()),
                                                   molecules_flipped=len(flips))) + '\n')
    write_provenance(d)
    print(f'SUMMARY f_flip bonds {sum(flips.values())} molecules {len(flips)}')


if __name__ == '__main__':
    main()
