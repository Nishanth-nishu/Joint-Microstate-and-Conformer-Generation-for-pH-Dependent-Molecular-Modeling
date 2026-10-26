# SAMPL8 physical-properties pKa data (repository inspection)

## Citation
samplchallenges/SAMPL8 GitHub repository, physical_properties/ (MIT license on repo).
Inspected directly via the GitHub API and raw files, 2026-10-26. Files copied to
/home2/nishanth.r/jmcg_stage/ext/sampl8/.

## Fulltext status
full_text_read (repository README + data files; the SAMPL8 pKa evaluation paper was NOT read).

## What is actually there
- `experimental_data/pKa.csv`: 23 compounds (SAMPL8-1…23), 1–2 macroscopic pKa each with
  uncertainty; 2 compounds N/A. Measured by Sirius T3 (see `SAMPL8 - SiriusT3 measurements.pdf`).
  **Macroscopic only — no microstate assignment, no populations.**
- `microstates/`: 278 microstate structure files in each of .sdf, .mol2 and .pdb (plus 25 CSVs of
  microstate IDs and SMILES), i.e. **one 3D structure per microstate**, generated with OpenEye
  toolkits from SMILES by `make_mol2_sdf_pdb_files.ipynb`. No energies, no ensembles, no weights.
  "The first microstate in each .CSV file indicated by SAMPL8-X_micro000 is our selected neutral
  reference state."
- Microstate lists are explicitly incomplete: "the microstate lists provided here may not be
  comprehensive", with later additions from Dortmund (Kast/Tielker) and Dassault, and known
  duplicates.

## Relevance
The best small, curated, openly licensed source of (microstate graph + explicit H + one 3D
structure) pairs — a ready-made evaluation set for the *structural* benchmark (task A), and a
reminder that even a curated microstate list is an enumeration, not ground truth.
