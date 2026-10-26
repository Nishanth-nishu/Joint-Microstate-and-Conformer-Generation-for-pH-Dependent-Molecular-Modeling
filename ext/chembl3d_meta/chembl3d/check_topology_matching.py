#!/usr/bin/env python3
"""Check matching between Zarr mol_id values and SDF topology names."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from sgdataset import SizeGroupedDataset


def _normalize_mol_id(mol_id) -> str:
    if isinstance(mol_id, (bytes, np.bytes_)):
        return bytes(mol_id).decode("utf-8", errors="replace").rstrip(" \x00")
    return str(mol_id).strip()


def _unique_mol_ids_ordered(mol_ids: np.ndarray) -> list[str]:
    if mol_ids.size == 0:
        return []
    start_mask = np.empty(mol_ids.size, dtype=bool)
    start_mask[0] = True
    start_mask[1:] = mol_ids[1:] != mol_ids[:-1]
    starts = np.flatnonzero(start_mask)
    return [_normalize_mol_id(mol_ids[i]) for i in starts]


def _load_sdf_names(sdf_file: Path) -> set[str]:
    names: set[str] = set()
    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=False)
    for mol in supplier:
        if mol is None:
            continue
        names.add(mol.GetProp("_Name").strip())
    return names


def _resolve_groups(ds: SizeGroupedDataset, group_arg: str) -> list[int]:
    if group_arg.lower() == "all":
        return sorted(ds.keys())
    return [int(x.strip()) for x in group_arg.split(",") if x.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--zarr_dataset_dir",
        type=str,
        default=str((SCRIPT_DIR.parent / "ChEMBL3D" / "zarr_database").resolve()),
    )
    parser.add_argument(
        "--topology_folder",
        type=str,
        default=str((SCRIPT_DIR.parent / "ChEMBL3D" / "topologies").resolve()),
    )
    parser.add_argument(
        "--groups",
        type=str,
        default="all",
        help="Comma-separated group ids (e.g. 027,028) or 'all'",
    )
    parser.add_argument(
        "--show_examples",
        type=int,
        default=5,
    )
    args = parser.parse_args()

    ds = SizeGroupedDataset(args.zarr_dataset_dir)
    groups = _resolve_groups(ds, args.groups)

    total_missing_in_sdf = 0
    total_missing_in_zarr = 0

    for g in groups:
        if g not in ds:
            print(f"group {g:03d}: missing in dataset")
            continue

        sdf_file = Path(args.topology_folder) / f"{g:03d}.sdf"
        if not sdf_file.exists():
            print(f"group {g:03d}: missing sdf file {sdf_file}")
            continue

        group = ds[g]
        zarr_ids = set(_unique_mol_ids_ordered(group["mol_id"][:]))
        sdf_ids = _load_sdf_names(sdf_file)

        missing_in_sdf = sorted(zarr_ids - sdf_ids)
        missing_in_zarr = sorted(sdf_ids - zarr_ids)
        total_missing_in_sdf += len(missing_in_sdf)
        total_missing_in_zarr += len(missing_in_zarr)

        print(
            f"group {g:03d}: zarr_mols={len(zarr_ids)} sdf_mols={len(sdf_ids)} "
            f"missing_in_sdf={len(missing_in_sdf)} missing_in_zarr={len(missing_in_zarr)}"
        )

        if missing_in_sdf:
            print("  sample missing_in_sdf:", ", ".join(missing_in_sdf[: args.show_examples]))
        if missing_in_zarr:
            print("  sample missing_in_zarr:", ", ".join(missing_in_zarr[: args.show_examples]))

    print("\nsummary:")
    print(f"total_missing_in_sdf={total_missing_in_sdf}")
    print(f"total_missing_in_zarr={total_missing_in_zarr}")


if __name__ == "__main__":
    main()
