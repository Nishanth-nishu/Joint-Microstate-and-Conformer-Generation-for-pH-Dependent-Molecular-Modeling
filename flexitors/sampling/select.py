"""Exp 5 Track P selectors: choose exactly n_out = 2K_i of the proposals, GT-free (EXP-5 §4.4, §4.5, §12.2, §12.4;
Q-P1-009 / Q-P2-009 / Q-P3-005 answers; port of round5/ledger/sample_data/select.py:13-41 with the fixes below).

Candidates (EXP-5 §12.4): torsion farthest-point sampling (FPS) from the lowest-energy proposal; FPS from a random
proposal; lowest MMFF94s energy; MMFF filter q = 0.25 (drop the 25 % highest-energy proposals, keep ≥ n_out) then FPS
from either start. ``random`` is the NFE control (P1). "Filter q = 0.5" is not a candidate: at c = 2 it keeps exactly
2K and equals lowest-energy (Q-P2-009). Distance: mean circular distance (1 − cos Δτ) over the molecule's heavy-atom
rotatable torsions (prototype ``select.py:13-14``). Energies: TD's MMFF94s single point
(``diffusion/likelihood.py:19-20``, E-SAMPLE-055); a failed force-field setup gives NaN, treated as +∞ and counted.
The selector signature carries no RMSD, GT or K other than ``n_out`` (EXP-5 §4.5; static test).

CLI ``features``: torsions, MMFF94s energies and their CPU seconds for every proposal of a pool (``confs.pkl``).
"""
from __future__ import annotations

import pickle
import time
from argparse import ArgumentParser

import numpy as np
from rdkit import Chem

from flexitors.geometry.internal import dihedral
from flexitors.geometry.torsion_types import heavy_rotatable_bonds, td_reference_quadruple
from flexitors.td import import_baseline

FILTER_Q = 0.25
SELECTORS = ('random', 'fps_lowE', 'fps_random', 'lowest_E', 'filter25_fps_lowE', 'filter25_fps_random')
CANDIDATES = SELECTORS[1:]  # EXP-5 §12.4: chosen on validation against random


def circular_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Mean (1 − cos Δτ) between torsion rows of ``a`` (p, m) and ``b`` (q, m) -> (p, q); zeros if m = 0."""
    if a.shape[1] == 0:
        return np.zeros((len(a), len(b)))
    return (1 - np.cos(a[:, None, :] - b[None, :, :])).mean(-1)


def fps(tors: np.ndarray, k: int, start: int) -> np.ndarray:
    """Greedy farthest-point sampling of k rows starting at row ``start``."""
    sel = [start]
    d = circular_distance(tors, tors[[start]])[:, 0]
    while len(sel) < k:
        d[sel] = -np.inf
        j = int(np.argmax(d))
        sel.append(j)
        d = np.minimum(d, circular_distance(tors, tors[[j]])[:, 0])
    return np.array(sel)


def select(method: str, tors: np.ndarray, energy: np.ndarray, n_out: int, rng: np.random.Generator) -> np.ndarray:
    """Indices (into the proposals) of the ``n_out`` conformers kept by ``method``."""
    n = len(tors)
    if n_out > n:
        raise ValueError(f'{n_out} outputs from {n} proposals')
    e = np.where(np.isnan(energy), np.inf, energy)
    if method == 'random':
        return np.sort(rng.choice(n, n_out, replace=False))
    if method == 'lowest_E':
        return np.argsort(e, kind='stable')[:n_out]
    keep = np.arange(n)
    if method.startswith('filter25'):
        keep = np.argsort(e, kind='stable')[:max(n_out, int(np.ceil((1 - FILTER_Q) * n)))]
    start = int(np.argmin(e[keep])) if method.endswith('lowE') else int(rng.integers(len(keep)))
    if method not in SELECTORS:
        raise ValueError(method)
    return keep[fps(tors[keep], n_out, start)]


def torsion_quads(mol: Chem.Mol) -> np.ndarray:
    """TD quadruples of the heavy-atom rotatable bonds (H rotors do not move heavy atoms)."""
    return np.array([td_reference_quadruple(mol, k) for k in heavy_rotatable_bonds(mol)], dtype=np.int64).reshape(-1, 4)


def features(mols: list[Chem.Mol]) -> dict:
    """Torsions (rad), MMFF94s single-point energies (NaN on failure) and the CPU seconds spent on them."""
    like = import_baseline('diffusion.likelihood')
    q = torsion_quads(mols[0])
    tors = np.stack([dihedral(m.GetConformer().GetPositions(), q) for m in mols]) if len(q) else np.zeros((len(mols), 0))
    t0 = time.process_time()
    energy = []
    for m in mols:
        try:
            energy.append(like.mmff_energy(m))
        except Exception:  # MMFF setup failure (RDKit raises several types); counted as NaN
            energy.append(np.nan)
    return dict(tors=tors, energy=np.array(energy, dtype=np.float64), mmff_cpu_s=time.process_time() - t0)


def main() -> None:
    ap = ArgumentParser(description='Exp 5 selector features of a pool')
    ap.add_argument('--confs', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    with open(a.confs, 'rb') as f:
        confs = pickle.load(f)
    feats = {smi: features(mols) for smi, mols in confs.items()}
    with open(a.out, 'wb') as f:
        pickle.dump(feats, f)
    n_nan = sum(int(np.isnan(v['energy']).sum()) for v in feats.values())
    print(f'SUMMARY features molecules {len(feats)}; MMFF NaN {n_nan}; '
          f'MMFF cpu s {sum(v["mmff_cpu_s"] for v in feats.values()):.1f}')


if __name__ == '__main__':
    main()
