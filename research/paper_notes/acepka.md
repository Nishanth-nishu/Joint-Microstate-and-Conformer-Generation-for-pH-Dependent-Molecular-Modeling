# AcepKa (PlayMolecule AI)
## Citation
Pesce F., Farr S.E., de Fabritiis G., "AcepKa: Thermodynamics-Informed pKa Prediction and Protonation-State
Generation in PlayMolecule AI", J. Chem. Inf. Model. 2026, DOI 10.1021/acs.jcim.6c03037; preprint arXiv:2604.00841
(v1 2026-04-01, v2 2026-08-31). Code/service: playmolecule.com web app (free account; no API; no licence stated; models not distributed). Checked 2026-10-26.
## Fulltext status
full_text_read — arXiv HTML v1 (earlier pass) and **v2 re-read 2026-10-26 10:40 UTC (16:10 IST)** via WebFetch of
https://arxiv.org/html/2604.00841v2 (Phase VIII retry). JCIM version of record not read (paywalled; not bypassed).
Section headings (v2): 1 Intro; 2 Methods (2.1 Microstate Population and Free Energy; 2.2 Architecture);
3 Results (3.1 Re-training and Model Evaluation; 3.2 App Workflow; 3.3 AceConfgen; 3.4 Implementation); 4 Conclusion; 6 Data and Software Availability.
## Problem
Macroscopic pKa prediction that is thermodynamically consistent across multiple ionizable sites, plus producing
the protonation state(s) of a molecule at a user-specified pH for downstream modelling.
## Method (mechanism level)
Implements the Uni-pKa framework. (1) Microstate enumeration is RULE-BASED: SMARTS templates identify ionizable
sites and the enumerator "combinatorially generate[s] all valid microstates with net charge bound between -2 and +2"
(Sec. 2.2 item 1). (2) A retrained Uni-Mol transformer backbone "predicts its standard Gibbs free energy" per
microstate (Sec. 2.2). (3) An FE2pKa module turns free energies into macro-pKa values and pH-dependent microstate
populations via Boltzmann weighting (Eqs. 1-2, Sec. 2.1). (4) A conformer engine (AceConfgen, GPU, reproduces RDKit
distance geometry + MMFF94, Sec. 3.2-3.3) generates "a conformational ensemble for each of the microstates"
AFTER enumeration (Fig. 2). A 3D-aware mode applies protonation to an existing bound pose without moving heavy atoms.
## Data
Uni-pKa training data (ChEMBL/Novartis-derived pKa sets as in Uni-pKa); conformers generated internally.
## Inputs/outputs
Input: one or more SMILES (or 3D conformers) AND a pH value (Fig. 2). Output: protonation ensemble with
pH-dependent microstate populations, macro-pKa, and 3D conformers per microstate.
## Experiments & baselines (v2, verified 2026-10-26)
- Training/fine-tuning (Sec. 3.1): Uni-pKa recipe re-run; pretraining on ~1M ChEMBL molecules (>3M protonation states),
  fine-tuning on DataWarrior + curated iBonD; "Fine-tuning produced a 5-fold ensemble model." (5-fold = model ensemble, not conformers.)
- pKa evaluation (Sec. 3.1, Fig. 1): macro-pKa **RMSE only** on "public datasets"; the body text does not name them
  (no occurrence of "SAMPL" or "Novartis" in v2); baseline numbers are copied from Uni-pKa [9] and Wagen 2025 (Rowan
  Starling) [11], so the sets are presumably the Uni-pKa/Starling ones (Novartis acid/base, SAMPL6/7/8) — Unknown;
  resolve by reading the Fig. 1 image. Claim: "accuracy on par with the state of the art"; differences "within the
  statistical uncertainty of these relatively small benchmark sets". No MAE, no micro-pKa, no tautomer evaluation.
- Conformer evaluation (Sec. 3.3, Fig. 4): Platinum 2017 (4,548 PDB ligands), 50 conformers/molecule, minimum RMSD to
  the bound pose; AceConfgen ≈ nvMolKit ≈ RDKit in accuracy, ~7x faster than other GPU implementations. This is a
  molecule-level, single-input-state benchmark: **no per-microstate conformer evaluation, no protonation-state
  recall/precision against reference states, no conformational-sensitivity analysis of pKa.**
- Conformers inside the pKa pipeline (Sec. 3.2, Fig. 2): "A conformational ensemble for each microstate is generated
  using AceConfgen" and "The 3D conformations are fed to the Uni-Mol model, which predicts their free energies".
  Number of conformers per microstate and the aggregation rule (mean vs Boltzmann/log-sum-exp vs single) are
  **NOT stated** in v2. Since AcepKa "inherits" Uni-pKa, the Uni-pKa code default (arithmetic mean of per-conformer
  free energies over 10 MMFF + 1 2D conformer; see unipka_2024.md) is the most likely behaviour — Unknown, not verified.
- pH (Sec. 2.1, Eq. 2): Boltzmann populations at query pH; library mode returns "the highest populated protonated
  state" per molecule at the target pH; no population threshold specified.
- Availability (Sec. 6): free web app (account required); "A programmatic API is not currently offered"; trained models
  for the Uni-pKa benchmarks "not distributed". **No software licence stated**; arXiv article licence CC BY-SA 4.0.
## Limitations
Enumeration is combinatorial/rule-based, not learned; charge range limited to [-2,+2]; conformers generated
per-microstate by a classical engine, so geometry does not feed back into state selection; protonation prediction
uses the ligand only, not the protein pocket (no pocket-induced pKa shifts).
## Relevance to p(z,x|G,h)
HIGHEST of all found work on the *interface* conditions: pH IS an explicit input and both z (microstate) and x
(conformers) are produced. But it is exactly the SEQUENTIAL PIPELINE the project proposes to beat: enumerate z
(rules) -> score z (energy model) -> generate x per z (classical conformer generator). No joint p(z,x|G,h), no
generative model over z, no coupling from geometry back to state.
## Exact supporting quotes
- Sec. 3.2: "Accurate pKa prediction requires reasonable ensembles of 3D conformers for the microstates."
- Sec. 3.3: "We generated 50 conformers per molecule for all 4,548 molecules (227,400 conformers in total)".
- Sec. 3.1: "are within the statistical uncertainty of these relatively small benchmark sets".
- Sec. 6: "AcepKa is freely available as a web application within PlayMolecule AI".
- Fig. 2: "takes in input one or multiple SMILEs strings (or 3D conformers) and a pH value".
- Sec. 2.2: "combinatorially generate all valid microstates with net charge bound between -2 and +2".
- Fig. 2: "A conformational ensemble for each of the microstates is generated".
