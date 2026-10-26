This README.txt file was generated on 20260309 by Filipp Nikitin


-------------------
GENERAL INFORMATION
-------------------

1. Title of Dataset:
ChEMBL3D: A Large-Scale Collection of QM-Accurate Conformers for ChEMBL Molecules


2. Author Information
<create a new entry for each additional author>

First Author Contact Information
    Name: Filipp Nikitin
    Institution: Carnegie Mellon University
    Address: N/A
    Email: fnikitin@andrew.cmu.edu
    Phone Number: N/A


Corresponding Author Contact Information
    Name: Olexandr Isayev
    Institution: Carnegie Mellon University
    Address: N/A
    Email: olexandr@olexandrisayev.com 
    Phone Number: N/A

Author Contact Information (if applicable)
    Name: Roman Zubatyuk
    Institution: NVIDIA
    Address: N/A
    Email: N/A
    Phone Number: N/A


---------------------
DATA & FILE OVERVIEW
---------------------


Directory of Files:
   A. Filename: topologies.tar.zst
      Short description: Compressed archive containing ChEMBL3D topology files (SDF/SMI) grouped by atom count.


   B. Filename: zarr_database.tar.zst.part-*
      Short description: Split compressed archive parts for the size-grouped Zarr conformer database.


   C. Filename: ARCHIVE_SHA256SUMS
      Short description: SHA256 checksums for all archive files for integrity verification.


   D. Filename: scripts/
      Short description: Python utility scripts and pinned requirements for loading, validating, and converting ChEMBL3D data.


Additional Notes on File Relationships, Context, or Content
(for example, if a user wants to reuse and/or cite your data,
what information would you want them to know?):
Topologies and Zarr groups are aligned by group ID (e.g., group 027 in Zarr corresponds to topologies for group 027). Molecule identifiers in Zarr mol_id fields match topology molecule names. See README.md for extraction and usage commands.

Documentation Note:
README.md is included and provides installation/setup details, archive pack/unpack commands, and script usage examples.



File Naming Convention:
- zarr_database.tar.zst.part-XX: split archive parts for the full Zarr database.
- topologies.tar.zst: complete topology archive.
- scripts/*.py: helper scripts for loading, validation, filtering, and RDKit conversion.
- scripts/requirements.txt: pinned dependency versions for scripts.
- Inside extracted topologies: NNN.sdf / NNN.smi where NNN is atom-count group.
- Inside extracted zarr_database: NNN/ where NNN is atom-count group.



-----------------------------------------
DATA DESCRIPTION FOR: zarr_database
-----------------------------------------
<create sections for each dataset included>


1. Number of variables:
7 core arrays per group (numbers, coord, energy, relative_energy, stereo_id, mol_id, charge)


2. Number of cases/rows:
283577147 conformers total across all groups


3. Missing data codes:
        Code/symbol        Definition
        N/A                No explicit missing-data code used
        N/A                N/A


4. Variable List



    A. Name: numbers
       Description: Atomic numbers for each conformer, shape (N, A), integer dtype.


    B. Name: coord
       Description: 3D coordinates, shape (N, A, 3), float dtype.

    C. Name: energy
       Description: Absolute conformer energy from optimization pipeline, shape (N,).

    D. Name: relative_energy
       Description: Energy relative to lowest conformer within corresponding (mol_id, stereo_id), shape (N,).

    E. Name: stereo_id
       Description: Stereochemistry class index within each mol_id, shape (N,).

    F. Name: mol_id
       Description: Molecule identifier string matching topology molecule names, shape (N,).

    G. Name: charge
       Description: Molecular charge for each conformer, shape (N,).

-----------------------------------------
DATA DESCRIPTION FOR: topologies
-----------------------------------------
<create sections for each dataset included>


1. Number of variables:
N/A (structure files; fields are embedded in SDF/SMI records)


2. Number of cases/rows:
N/A (depends on group file; one entry per molecule topology)


3. Missing data codes:
        Code/symbol        Definition
        N/A                N/A
        N/A                N/A


4. Variable List



    A. Name: _Name
       Description: Molecule identifier in SDF entries; should match Zarr mol_id.


    B. Name: Atom/Bond block
       Description: Molecular topology with explicit hydrogens and connectivity.

-----------------------------------------
DATA DESCRIPTION FOR: scripts
-----------------------------------------
<create sections for each dataset included>


1. Number of variables:
N/A (code files, not tabular data)


2. Number of cases/rows:
N/A


3. Missing data codes:
        Code/symbol        Definition
        N/A                N/A
        N/A                N/A


4. Variable List


    A. Name: scripts/sgdataset.py
       Description: Core size-grouped Zarr data loader utilities.


    B. Name: scripts/read_dataset_overview.py
       Description: Prints overall and per-group dataset statistics.

    C. Name: scripts/check_topology_matching.py
       Description: Verifies consistency between Zarr mol_id and SDF topology names.

    D. Name: scripts/molecule_to_rdkit_multiconf.py
       Description: Builds one RDKit molecule with multiple conformers for a selected molecule.

    E. Name: scripts/group_to_rdkit_multiconf.py
       Description: Builds RDKit multi-conformer molecules for multiple molecules in a selected group.

    F. Name: scripts/filter_zarr_by_relative_energy.py
       Description: Creates a filtered copy of Zarr data by relative energy threshold.

    G. Name: scripts/requirements.txt
       Description: Pinned package versions required by the scripts.

--------------------------
METHODOLOGICAL INFORMATION
--------------------------


1. Software-specific information:
<create a new entry for each qualifying software program>

Name: Python scripts for dataset access/processing
Version: Python 3.10 environment (project scripts)
System Requirements: Linux/macOS/Windows with Python and required packages
Open Source? (Y/N): Y

(if available and applicable)
Executable URL: N/A
Source Repository URL: N/A
Developer: N/A
Product URL: N/A
Software source components: numpy, zarr, numcodecs, h5py, torch, rdkit, tqdm


Additional Notes(such as, will this software not run on
certain operating systems?):
Core scripts are cross-platform Python, but data processing performance depends on storage throughput and available memory.


2. Equipment-specific information:
<create a new entry for each qualifying piece of equipment>

Manufacturer: N/A
Model: N/A

(if applicable)
Embedded Software / Firmware Name: N/A
Embedded Software / Firmware Version: N/A
Additional Notes: Data were produced using HPC compute systems.


3. Date of data collection (single date, range, approximate date) <suggested format YYYYMMDD>:
Approximate range: 20240101-20260228; release posted 20260301.
