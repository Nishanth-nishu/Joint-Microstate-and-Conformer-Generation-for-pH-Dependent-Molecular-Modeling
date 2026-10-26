#!/usr/bin/env python3
"""Load one molecule from one group and build a single RDKit Mol with all conformers."""

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


def _load_topologies(sdf_file: Path) -> dict[str, Chem.Mol]:
    out = {}
    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=False)
    for mol in supplier:
        if mol is None:
            continue
        out[mol.GetProp("_Name").strip()] = mol
    return out


def _find_molecule_range(mol_ids: np.ndarray, mol_id: str) -> tuple[int, int] | None:
    normalized = np.array([_normalize_mol_id(x) for x in mol_ids], dtype=object)
    idx = np.flatnonzero(normalized == mol_id)
    if idx.size == 0:
        return None
    return int(idx[0]), int(idx[-1]) + 1


def _build_multiconf_mol(
    topology: Chem.Mol,
    coords: np.ndarray,
    energies: np.ndarray | None,
    rel_energies: np.ndarray | None,
    mol_id: str,
) -> Chem.Mol:
    mol = Chem.Mol(topology)
    mol.SetProp("_Name", mol_id)
    mol.RemoveAllConformers()

    natoms = mol.GetNumAtoms()
    for i, xyz in enumerate(coords):
        conf = Chem.Conformer(natoms)
        for a in range(natoms):
            x, y, z = map(float, xyz[a])
            conf.SetAtomPosition(a, Chem.rdGeometry.Point3D(x, y, z))
        cid = mol.AddConformer(conf, assignId=True)
        if energies is not None:
            mol.GetConformer(cid).SetDoubleProp("energy", float(energies[i]))
        if rel_energies is not None:
            mol.GetConformer(cid).SetDoubleProp("relative_energy", float(rel_energies[i]))

    return mol


def _write_conformers_to_sdf(multiconf_mol: Chem.Mol, out_sdf: Path) -> None:
    writer = Chem.SDWriter(str(out_sdf))
    try:
        for conf in multiconf_mol.GetConformers():
            conf_id = conf.GetId()
            mol = Chem.Mol(multiconf_mol)
            single_conf = Chem.Conformer(mol.GetConformer(conf_id))
            mol.RemoveAllConformers()
            mol.AddConformer(single_conf, assignId=True)
            mol.SetProp("conf_id", str(conf_id))
            if conf.HasProp("energy"):
                mol.SetDoubleProp("energy", conf.GetDoubleProp("energy"))
            if conf.HasProp("relative_energy"):
                mol.SetDoubleProp("relative_energy", conf.GetDoubleProp("relative_energy"))
            writer.write(mol)
    finally:
        writer.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr_dataset_dir", type=str, default=str((SCRIPT_DIR.parent / "ChEMBL3D" / "zarr_database").resolve()))
    parser.add_argument("--topology_folder", type=str, default=str((SCRIPT_DIR.parent / "ChEMBL3D" / "topologies").resolve()))
    parser.add_argument("--group", type=int, required=True, help="Group id, e.g. 27")
    parser.add_argument("--mol_id", type=str, default=None, help="Molecule id from mol_id column")
    parser.add_argument("--max_confs", type=int, default=0, help="If >0, truncate to first N conformers")
    parser.add_argument("--out_sdf", type=str, default="", help="Optional output SDF path")
    args = parser.parse_args()

    ds = SizeGroupedDataset(args.zarr_dataset_dir)
    group = ds[args.group]

    mol_ids = group["mol_id"][:]
    if args.mol_id is None:
        args.mol_id = _normalize_mol_id(mol_ids[0])

    span = _find_molecule_range(mol_ids, args.mol_id)
    if span is None:
        raise ValueError(f"mol_id '{args.mol_id}' is not present in group {args.group:03d}")

    start, end = span
    if args.max_confs > 0:
        end = min(end, start + int(args.max_confs))

    sdf_file = Path(args.topology_folder) / f"{args.group:03d}.sdf"
    topologies = _load_topologies(sdf_file)
    if args.mol_id not in topologies:
        raise KeyError(f"mol_id '{args.mol_id}' is missing in topology file {sdf_file}")

    coords = group["coord"][start:end]
    energies = group["energy"][start:end] if "energy" in group.keys() else None
    rel = group["relative_energy"][start:end] if "relative_energy" in group.keys() else None

    multiconf_mol = _build_multiconf_mol(topologies[args.mol_id], coords, energies, rel, args.mol_id)
    nconf = multiconf_mol.GetNumConformers()

    print(f"group={args.group:03d}")
    print(f"mol_id={args.mol_id}")
    print(f"selected_rows=[{start}:{end})")
    print(f"num_conformers={nconf}")

    if args.out_sdf:
        out_sdf = Path(args.out_sdf)
        out_sdf.parent.mkdir(parents=True, exist_ok=True)
        _write_conformers_to_sdf(multiconf_mol, out_sdf)
        print(f"wrote_sdf={out_sdf}")


if __name__ == "__main__":
    main()
