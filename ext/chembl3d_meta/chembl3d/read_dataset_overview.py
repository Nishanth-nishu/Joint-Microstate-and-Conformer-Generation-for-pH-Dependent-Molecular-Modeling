#!/usr/bin/env python3
"""Quick overview of ChEMBL3D Zarr dataset grouped by atom count."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from sgdataset import SizeGroupedDataset


def _unique_molecule_count(mol_ids: np.ndarray) -> int:
    if mol_ids.size == 0:
        return 0
    start_mask = np.empty(mol_ids.size, dtype=bool)
    start_mask[0] = True
    start_mask[1:] = mol_ids[1:] != mol_ids[:-1]
    return int(np.count_nonzero(start_mask))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--zarr_dataset_dir",
        type=str,
        default=str((SCRIPT_DIR.parent / "ChEMBL3D" / "zarr_database").resolve()),
        help="Path to SizeGroupedDataset root",
    )
    parser.add_argument(
        "--show_groups",
        type=int,
        default=10,
        help="How many groups to print from the beginning",
    )
    args = parser.parse_args()

    ds = SizeGroupedDataset(args.zarr_dataset_dir)
    group_keys = sorted(ds.keys())

    print(f"dataset: {args.zarr_dataset_dir}")
    print(f"groups: {len(group_keys)}")

    total_rows = 0
    total_molecules = 0
    for g in group_keys:
        group = ds[g]
        rows = len(group["energy"])
        mols = _unique_molecule_count(group["mol_id"][:])
        total_rows += rows
        total_molecules += mols

    print(f"total_conformers: {total_rows}")
    print(f"total_molecules: {total_molecules}")

    print("\nfirst groups:")
    print("group\trows\tmolecules\tkeys")
    for g in group_keys[: max(0, int(args.show_groups))]:
        group = ds[g]
        rows = len(group["energy"])
        mols = _unique_molecule_count(group["mol_id"][:])
        keys = ",".join(group.keys())
        print(f"{g:03d}\t{rows}\t{mols}\t{keys}")


if __name__ == "__main__":
    main()
