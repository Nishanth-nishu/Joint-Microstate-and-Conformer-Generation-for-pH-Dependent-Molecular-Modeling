This README.txt file was generated on 20260309 by Filipp Nikitin

-------------------
GENERAL INFORMATION
-------------------

1. Title of Dataset:

LoQI data release (model checkpoints and ChEMBL3D stereo processed data)

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
    Email: olexandr@cmu.edu
	Phone Number: N/A

Author Contact Information (if applicable)
    Name: N/A
    Institution: N/A
    Address: N/A
    Email: N/A
	Phone Number: N/A

---------------------
DATA & FILE OVERVIEW
---------------------

Directory of Files:
   A. Filename: README.txt
      Short description: This metadata file in README.txt template format.

   B. Filename: N/A
      Short description: No additional README.md is included in this release.

   C. Filename: loqi.ckpt
      Short description: LoQI checkpoint (diffusion model).

   D. Filename: loqi_flow.ckpt
      Short description: LoQI flow-matching checkpoint.

   E. Filename: sha256sums.txt
      Short description: SHA256 checksums for files in this data release.

   F. Filename: chembl3d_stereo/processed/train_h.pt
      Short description: Processed PyTorch Geometric objects for train split.

   G. Filename: chembl3d_stereo/processed/val_h.pt
      Short description: Processed PyTorch Geometric objects for validation split.

   H. Filename: chembl3d_stereo/processed/test_h.pt
      Short description: Processed PyTorch Geometric objects for test split.

   I. Filename: chembl3d_stereo/processed/test_cremp_h.pt
      Short description: Additional processed test subset (CREMP subset).

   J. Filename: chembl3d_stereo/processed/test_rot_bonds_h.pt
      Short description: Processed test subset emphasizing rotatable-bond handling.

   K. Filename: chembl3d_stereo/processed/test_small_h.pt
      Short description: Processed test subset of small molecules.

   L. Filename: chembl3d_stereo/processed/train_atom_types_h.npy
      Short description: Train split atom-type array.

   M. Filename: chembl3d_stereo/processed/val_atom_types_h.npy
      Short description: Validation split atom-type array.

   N. Filename: chembl3d_stereo/processed/test_atom_types_h.npy
      Short description: Test split atom-type array.

   O. Filename: chembl3d_stereo/processed/train_bond_types_h.npy
      Short description: Train split bond-type array.

   P. Filename: chembl3d_stereo/processed/val_bond_types_h.npy
      Short description: Validation split bond-type array.

   Q. Filename: chembl3d_stereo/processed/test_bond_types_h.npy
      Short description: Test split bond-type array.

   R. Filename: chembl3d_stereo/processed/train_charges_h.npy
      Short description: Train split formal-charge array.

   S. Filename: chembl3d_stereo/processed/val_charges_h.npy
      Short description: Validation split formal-charge array.

   T. Filename: chembl3d_stereo/processed/test_charges_h.npy
      Short description: Test split formal-charge array.

   U. Filename: chembl3d_stereo/processed/train_charges_prior_h.npy
      Short description: Train split prior charge distribution array.

   V. Filename: chembl3d_stereo/processed/val_charges_prior_h.npy
      Short description: Validation split prior charge distribution array.

   W. Filename: chembl3d_stereo/processed/test_charges_prior_h.npy
      Short description: Test split prior charge distribution array.

   X. Filename: chembl3d_stereo/processed/train_hybridization_h.npy
      Short description: Train split hybridization array.

   Y. Filename: chembl3d_stereo/processed/val_hybridization_h.npy
      Short description: Validation split hybridization array.

   Z. Filename: chembl3d_stereo/processed/test_hybridization_h.npy
      Short description: Test split hybridization array.

   AA. Filename: chembl3d_stereo/processed/train_is_aromatic_h.npy
      Short description: Train split aromaticity flags.

   AB. Filename: chembl3d_stereo/processed/val_is_aromatic_h.npy
      Short description: Validation split aromaticity flags.

   AC. Filename: chembl3d_stereo/processed/test_is_aromatic_h.npy
      Short description: Test split aromaticity flags.

   AD. Filename: chembl3d_stereo/processed/train_is_in_ring_h.npy
      Short description: Train split ring-membership flags.

   AE. Filename: chembl3d_stereo/processed/val_is_in_ring_h.npy
      Short description: Validation split ring-membership flags.

   AF. Filename: chembl3d_stereo/processed/test_is_in_ring_h.npy
      Short description: Test split ring-membership flags.

   AG. Filename: chembl3d_stereo/processed/train_angles_h.pickle
      Short description: Train split angle-statistics data.

   AH. Filename: chembl3d_stereo/processed/val_angles_h.pickle
      Short description: Validation split angle-statistics data.

   AI. Filename: chembl3d_stereo/processed/test_angles_h.pickle
      Short description: Test split angle-statistics data.

   AJ. Filename: chembl3d_stereo/processed/train_bond_lengths_h.pickle
      Short description: Train split bond-length statistics.

   AK. Filename: chembl3d_stereo/processed/val_bond_lengths_h.pickle
      Short description: Validation split bond-length statistics.

   AL. Filename: chembl3d_stereo/processed/test_bond_lengths_h.pickle
      Short description: Test split bond-length statistics.

   AM. Filename: chembl3d_stereo/processed/train_dihedrals_h.pickle
      Short description: Train split dihedral statistics.

   AN. Filename: chembl3d_stereo/processed/val_dihedrals_h.pickle
      Short description: Validation split dihedral statistics.

   AO. Filename: chembl3d_stereo/processed/test_dihedrals_h.pickle
      Short description: Test split dihedral statistics.

   AP. Filename: chembl3d_stereo/processed/train_n_h.pickle
      Short description: Train split molecule-size metadata.

   AQ. Filename: chembl3d_stereo/processed/val_n_h.pickle
      Short description: Validation split molecule-size metadata.

   AR. Filename: chembl3d_stereo/processed/test_n_h.pickle
      Short description: Test split molecule-size metadata.

   AS. Filename: chembl3d_stereo/processed/train_smiles.pickle
      Short description: Train split SMILES mapping.

   AT. Filename: chembl3d_stereo/processed/val_smiles.pickle
      Short description: Validation split SMILES mapping.

   AU. Filename: chembl3d_stereo/processed/test_smiles.pickle
      Short description: Test split SMILES mapping.

   AV. Filename: chembl3d_stereo/processed/train_valency_h.pickle
      Short description: Train split valency metadata.

   AW. Filename: chembl3d_stereo/processed/val_valency_h.pickle
      Short description: Validation split valency metadata.

   AX. Filename: chembl3d_stereo/processed/test_valency_h.pickle
      Short description: Test split valency metadata.

Additional Notes on File Relationships, Context, or Content
(for example, if a user wants to reuse and/or cite your data,
what information would you want them to know?):

This release contains pretrained model checkpoints plus a processed ChEMBL3D stereo dataset used by LoQI.
The processed files are intended to be consumed together (train/val/test split-wise). Checksums are provided
in sha256sums.txt to verify file integrity after transfer.

File Naming Convention:

- Split prefix: train_, val_, test_
- Suffix _h indicates hydrogen-inclusive processing.
- Extension meanings:
  - .pt: serialized PyTorch/PyG tensors/objects
  - .npy: NumPy arrays
  - .pickle: Python pickle metadata/statistics
  - .ckpt: PyTorch Lightning checkpoint
  - .txt: plain text checksums/metadata

-----------------------------------------
DATA DESCRIPTION FOR: sha256sums.txt
-----------------------------------------
<create sections for each dataset included>

1. Number of variables:
2

2. Number of cases/rows:
48

3. Missing data codes:
        Code/symbol        Definition
        N/A                No missing values expected

4. Variable List

    A. Name: sha256
       Description: SHA256 hash value for each file in this release.

    B. Name: filepath
       Description: Relative path to file corresponding to sha256 hash.

-----------------------------------------
DATA DESCRIPTION FOR: loqi.ckpt, loqi_flow.ckpt
-----------------------------------------

1. Number of variables:
N/A (binary model checkpoint files)

2. Number of cases/rows:
N/A

3. Missing data codes:
        Code/symbol        Definition
        N/A                Not applicable

4. Variable List

    A. Name: N/A
       Description: Model weights and training-state tensors serialized by PyTorch Lightning.

    B. Name: N/A
       Description: Not a tabular dataset.

-----------------------------------------
DATA DESCRIPTION FOR: chembl3d_stereo/processed/*
-----------------------------------------

1. Number of variables:
N/A (mixed binary arrays/tensors/pickles)

2. Number of cases/rows:
N/A (depends on file/split and object type)

3. Missing data codes:
        Code/symbol        Definition
        N/A                Not encoded as tabular missing codes

4. Variable List

    A. Name: split
       Description: train/val/test partition.

    B. Name: feature family
       Description: atom_types, bond_types, charges, hybridization, aromaticity, ring flags, angles, bond lengths,
                    dihedrals, valency, smiles mapping, or graph objects.

--------------------------
METHODOLOGICAL INFORMATION
--------------------------

1. Software-specific information:
<create a new entry for each qualifying software program>

Name: LoQI
Version: N/A
System Requirements: Python environment with PyTorch and RDKit
Open Source? (Y/N): Y

(if available and applicable)
Executable URL: N/A
Source Repository URL: https://github.com/isayevlab/LoQI
Developer: IsayevLab
Product URL: https://github.com/isayevlab/LoQI
Software source components: N/A

Additional Notes(such as, will this software not run on
certain operating systems?):
N/A

Name: Python
Version: N/A
System Requirements: Linux/Windows/macOS
Open Source? (Y/N): Y

(if available and applicable)
Executable URL: https://www.python.org/
Source Repository URL: https://github.com/python/cpython
Developer: Python Software Foundation
Product URL: https://www.python.org/
Software source components: N/A

Additional Notes(such as, will this software not run on
certain operating systems?):
N/A

Name: PyTorch / PyTorch Lightning
Version: N/A
System Requirements: CPU or CUDA-compatible GPU for training/inference workloads
Open Source? (Y/N): Y

(if available and applicable)
Executable URL: N/A
Source Repository URL: https://github.com/pytorch/pytorch
Developer: PyTorch contributors
Product URL: https://pytorch.org/
Software source components: N/A

Additional Notes(such as, will this software not run on
certain operating systems?):
GPU acceleration depends on CUDA toolchain compatibility.

Name: RDKit
Version: N/A
System Requirements: Standard Python scientific stack
Open Source? (Y/N): Y

(if available and applicable)
Executable URL: N/A
Source Repository URL: https://github.com/rdkit/rdkit
Developer: RDKit contributors
Product URL: https://www.rdkit.org/
Software source components: N/A

Additional Notes(such as, will this software not run on
certain operating systems?):
N/A

2. Equipment-specific information:
<create a new entry for each qualifying piece of equipment>

Manufacturer: N/A
Model: N/A

(if applicable)
Embedded Software / Firmware Name: N/A
Embedded Software / Firmware Version: N/A
Additional Notes: N/A

3. Date of data collection (single date, range, approximate date) <suggested format YYYYMMDD>:

20240101 to 20260309 (authors: Filipp Nikitin, Olexandr Isayev)
