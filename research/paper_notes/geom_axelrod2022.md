# GEOM (Axelrod & Gómez-Bombarelli 2022)

## Citation
GEOM: Energy-annotated molecular conformations for property prediction and molecular generation.
Simon Axelrod, Rafael Gómez-Bombarelli. 2022. Scientific Data 9:185. DOI 10.1038/s41597-022-01288-4.
Full text read via arXiv:2006.05531v4 HTML. Code/data: https://github.com/learningmatter-mit/geom. Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v4 HTML; publisher HTML blocked by idp.nature.com redirect).

## Problem
Provide a large energy-annotated 3D conformer dataset for property prediction and generative modelling.

## Method (mechanism level)
CREST metadynamics with GFN2-xTB produces conformer-rotamer ensembles (CREs) per molecule. Statistical weights
p_i ∝ d_i exp(-E_i/k_BT) (Eq. 1). Subsets re-ranked with single-point DFT and xTB quasi-harmonic free energies
(CENSO-style, Eqs. 3-4). Graphs of generated conformers re-identified with xyz2mol.

## Data
133,258 QM9 molecules; 304,466 drug-like (AICures) molecules; plus MoleculeNet subsets and 1,513 BACE species.
Vacuum CREs throughout; BACE additionally has CREs with an implicit water model.

## Inputs/outputs — protonation / charge / H handling (the key section for this project)
- Charged species ARE present and the ionization state is fixed per molecule by preprocessing:
  "3.9% of the drug molecules accessed (11,886 total) were given as clusters, either with a counterbalancing ion
  (e.g. '.[Na+]', '.[Cl-]') or with an acid to represent the protonated salt (e.g. '.Cl'). ... For the acid/base
  SMILES strings, used reaction SMARTS in RDKit to generate the protonated molecule and counterion." (Methods)
- "Not only does de-salting identify the drug-like compound in each cluster and correct its ionization state,
  it also homogenizes the molecular representations in the drug datasets." (Methods)
- For the MoleculeNet portion: "we also selected the heaviest component from each cluster SMILES, but did not
  perform protonation." (Methods) — i.e. protonation handling is NOT uniform across GEOM subsets.
- Charge enters the physics only as a per-molecule total: "Default values were used for all CREST arguments,
  except for the charge of each geometry." (Methods)
- Tautomers explicitly NOT enumerated: "CREST has specific commands to generate tautomers, but they were not
  used here." (Technical Validation / graph re-identification)
- No pH variable anywhere in the dataset; no multiple protomers of the same compound are provided as alternatives.
- Reactive events do occur in simulation: 21% of graph re-attribution failures were tautomerization (~1% of all cases).

## Experiments & baselines
Technical validation: xTB vs r2scan-3c relative energies MAE 1.96 kcal/mol; Spearman ρ per species 0.39±0.35.
Single-point DFT vs fully optimized: MAE 0.54 kcal/mol, mean RMSD 0.36 Å.

## Limitations
"the CREST statistical weights for each conformer are rather inaccurate. Therefore, benchmarks that include
conformer probabilities should use the DFT weights provided in GEOM." (Background/Usage) — DFT weights exist only
for the BACE subset at scale, not for all of GEOM-DRUGS.

## Relevance to p(z, x | G, h)
Decisive for feasibility: GEOM gives exactly ONE chemical graph (one protomer) per compound, chosen by a
deterministic de-salting rule, with no pH label and no tautomer alternatives. Supervision for a latent
protonation-microstate variable z therefore does NOT exist in GEOM and would have to be created (e.g. by
enumerating protomers with a pKa predictor and re-running CREST per protomer). GEOM does supply per-conformer
energies and statistical weights, which is the one ingredient a population/Boltzmann-weight head could use.
Also: vacuum ensembles are the default, which is physically awkward for a pH-conditioned (aqueous) model.
