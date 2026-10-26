# AcepKa (PlayMolecule AI)
## Citation
Pesce F., Farr S.E., de Fabritiis G., "AcepKa: Thermodynamics-Informed pKa Prediction and Protonation-State
Generation in PlayMolecule AI", J. Chem. Inf. Model. 2026, DOI 10.1021/acs.jcim.6c03037; preprint arXiv:2604.00841
(v1 2026-04-01, v2 2026-08-31). Code/service: playmolecule.com (proprietary). Checked 2026-10-26.
## Fulltext status
full_text_read (arXiv HTML v1, via WebFetch; PDF also in research/papers/arxiv_2604.00841.pdf but no local PDF text extractor).
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
## Experiments & baselines
pKa benchmarks vs Uni-pKa/Epik-type baselines (benchmark detail not fully extracted).
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
- Fig. 2: "takes in input one or multiple SMILEs strings (or 3D conformers) and a pH value".
- Sec. 2.2: "combinatorially generate all valid microstates with net charge bound between -2 and +2".
- Fig. 2: "A conformational ensemble for each of the microstates is generated".
