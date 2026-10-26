"""Per-conformer random streams (round5/DECISION.md graft 1: "per-molecule RNG so arms stay paired within each
molecule"; plan_P3 §3.2; EXP-3 §9.1 condition 1; EXP-5 §8 "hash(global seed, molecule, batch)").

Every conformer j of molecule ``key`` (its ``corrected_smiles``) owns one 63-bit seed h(global_seed, key, j). It drives
that conformer's ETKDG embedding, its initial torsion randomisation and every reverse-step noise draw. Consequences:
* a conformer's trajectory does not depend on the batch it shares, the molecule subset of the run, the pool size or
  the arm, so a 2K run is exactly the prefix of a max(128, 4K) pool (EXP-5 §12.1) and an arm run on the affected
  molecules can be spliced into O (EXP-3 §9.1);
* TD's single global seed (``generate_confs.py:94-95``) would desynchronise every later molecule as soon as an arm
  draws a different number of random values (plan_P2 §3.3).
The key is the SMILES, not the csv row, so subsets reproduce the full run (plan_P3 §3.2).
"""
from __future__ import annotations

import hashlib

import numpy as np

STREAMS = ('etkdg', 'perturb', 'noise')


def conformer_seed(global_seed: int, key: str, index: int, stream: str) -> int:
    """63-bit seed of one conformer's random stream (``stream`` in STREAMS, or an arm-specific name)."""
    h = hashlib.sha256(f'{global_seed}|{key}|{index}|{stream}'.encode()).hexdigest()
    return int(h[:15], 16)


def etkdg_seed(global_seed: int, key: str, index: int) -> int:
    """ETKDG ``randomSeed`` in [1, 2^31 - 1): never 0, which makes RDKit conformers identical
    (``generate_confs.py:64-67``)."""
    return 1 + conformer_seed(global_seed, key, index, 'etkdg') % (2 ** 31 - 2)


def numpy_rng(global_seed: int, key: str, index: int, stream: str) -> np.random.Generator:
    """Independent numpy generator of one conformer stream."""
    return np.random.default_rng(conformer_seed(global_seed, key, index, stream))
