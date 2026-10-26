# ChEMBL3D Dataset

Dataset release: ChEMBL3D

- Dataset posted on: `March 1, 2026`
- Authors: `Filipp Nikitin, Roman Zubatyuk, Olexandr Isayev`

## Brief Description of Dataset Files

- `ChEMBL3D/topologies/`
Contains SDF topology files grouped by atom count (`007.sdf`, `008.sdf`, ..., `142.sdf`).
Each SDF entry name (`_Name`) is the molecule key used to match conformers in Zarr.

- `ChEMBL3D/zarr_database/`
Contains the size-grouped Zarr conformer database. Each top-level group key is an atom-count group (`007`, `008`, ...).

- `scripts/sgdataset.py`
Standalone loader for the size-grouped Zarr format.

- `scripts/read_dataset_overview.py`
Prints global/group statistics and available fields.

- `scripts/check_topology_matching.py`
Checks consistency between Zarr `mol_id` values and SDF `_Name` values.

- `scripts/molecule_to_rdkit_multiconf.py`
Loads one molecule and creates one RDKit `Mol` containing all conformers.
Example:
```bash
python scripts/molecule_to_rdkit_multiconf.py --group 27 --mol_id CHEMBL100101_0 --max_confs 3 --out_sdf /tmp/chembl100101_readme_example.sdf
```

- `scripts/group_to_rdkit_multiconf.py`
Loads multiple molecules from a group and builds one RDKit multi-conformer `Mol` per molecule.

- `scripts/filter_zarr_by_relative_energy.py`
Creates a filtered copy of the Zarr dataset using a relative-energy threshold.

- `scripts/requirements.txt`
Pinned package versions for `sgdataset.py` and all example scripts.

`ChEMBL3D/topologies` and `ChEMBL3D/zarr_database` are release directories.

## Dataset Description

The dataset is curated from ChEMBL v34. The processing workflow includes:

1. Protomer generation (OpenEye FixpKa) and merge with original structures.
2. Enumeration of stereocenters with missing annotation (OpenEye Flipper).
3. Topology generation with hydrogens (Open Babel).
4. Conformer generation (OpenEye Omega Classic).
5. Geometry optimization using FIRE with AIMNet2 potential (toward DFT/CPCM quality).
6. Topology validation, filtering, and deduplication.

Supported elements (AIMNet2): `H, C, N, O, F, S, Cl, Br, I, B, Si, P, As, Se`.

Reported full-pipeline curation statistics:

- Initial conformers: `~505M`
- Processed conformers: `~280M`
- Unique structures: `~2.4M` (`~1.8M` unique ChEMBL compound IDs)
- Conformers with `relative_energy > 6 kcal/mol` were removed.

Current `zarr_database` in this release:

## Zarr Format (Per Group)

Each group stores conformer-level arrays with aligned first dimension `N`:

- `numbers`: `(N, A)` `int8` - atomic numbers for each conformer (group atom count is `A`).
- `coord`: `(N, A, 3)` `float32` - 3D coordinates for each conformer.
- `energy`: `(N,)` `float64` - absolute conformer energy from optimization pipeline.
- `relative_energy`: `(N,)` `float64` - energy shift relative to the lowest conformer in the corresponding `(mol_id, stereo_id)` block.
- `stereo_id`: `(N,)` `int32` - stereochemistry class index enumerated within each `mol_id`.
- `mol_id`: `(N,)` fixed-width bytes string (e.g. `|S15`) - conformer parent identifier; matches SDF `_Name`.
- `charge`: `(N,)` `int8` - molecular charge value stored by the pipeline.

Group key `027` example:

- `coord`: `(632919, 27, 3)`
- `mol_id`: `(632919,)`
- molecules in group: `18442`

## Matching Topology and Conformers

Matching rule:

```python
zarr_mol_id = bytes(mol_id).decode("utf-8", errors="replace").rstrip(" \x00")
sdf_name = sdf_mol.GetProp("_Name").strip()
assert zarr_mol_id == sdf_name
```

Group-level mapping:

- Zarr group `NNN` corresponds to topology file `topologies/NNN.sdf`.

### Concrete Matching Example

For group `027`:

- Zarr rows `[0:4)` -> `mol_id = CHEMBL100101_0` (`4` conformers)
- Zarr rows `[4:18)` -> `mol_id = CHEMBL100123_0` (`14` conformers)
- Zarr rows `[18:19)` -> `mol_id = CHEMBL100203_0` (`1` conformer)

Corresponding SDF entry in `topologies/027.sdf`:

```text
CHEMBL100101_0
     RDKit          3D

 27 30  0  0  0  0  0  0  0  0999 V2000
   -0.8492   -0.7690    2.4340 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1940   -0.6265    3.3517 C   0  0  0  0  0  0  0  0  0  0  0  0
    ...
M  END
$$$$
```

## Installation

Run from the repository root (the folder that contains `ChEMBL3D/` and `scripts/`):

```bash
conda create -y -n chembl3d python=3.10
conda activate chembl3d
pip install -r scripts/requirements.txt
```

## Usage

Run from the same repository root:

```bash
python scripts/read_dataset_overview.py
python scripts/check_topology_matching.py --groups 027,028
python scripts/molecule_to_rdkit_multiconf.py --group 27 --mol_id CHEMBL100101_0 --out_sdf /tmp/chembl100101.sdf
python scripts/group_to_rdkit_multiconf.py --group 27 --limit_molecules 10
```

Quick smoke test on a small group (`010`):

```bash
python scripts/check_topology_matching.py --groups 010
python scripts/molecule_to_rdkit_multiconf.py --group 10 --max_confs 3 --out_sdf /tmp/chembl3d_group010_one_mol.sdf
python scripts/group_to_rdkit_multiconf.py --group 10 --limit_molecules 10
```

## Archive / Unarchive

To package for platforms with file-count limits, run from repository root.

Archive `zarr_database` into 10 GB parts:

```bash
mkdir -p archives
tar -C ChEMBL3D -cf - zarr_database | zstd -19 -T0 | split -d -b 10G - archives/zarr_database.tar.zst.part-
```

Archive `topologies` separately:

```bash
tar -C ChEMBL3D -I 'zstd -19 -T0' -cf archives/topologies.tar.zst topologies
```

To unpack:

```bash
mkdir -p ChEMBL3D
cat archives/zarr_database.tar.zst.part-* | zstd -d | tar -xf - -C ChEMBL3D
tar -I zstd -xf archives/topologies.tar.zst -C ChEMBL3D
```

## License

This dataset is published under the MIT License.

## Citation

```
@article{
doi:10.26434/chemrxiv-2025-k4h7v,
author = {Filipp Nikitin  and Dylan M. Anstine  and Roman Zubatyuk  and Saee Gopal Paliwal  and Olexandr Isayev },
title = {Scalable Low-Energy Molecular Conformer Generation with Quantum Mechanical Accuracy},
journal = {ChemRxiv},
volume = {2025},
number = {0820},
pages = {},
year = {2025},
doi = {10.26434/chemrxiv-2025-k4h7v},
URL = {https://chemrxiv.org/doi/abs/10.26434/chemrxiv-2025-k4h7v},
eprint = {https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2025-k4h7v}}

```
