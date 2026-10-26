# Uni-pKa (JACS Au 2024)

## Citation
Zheng H., Luo W., Zhou G., Zhu Z., Yuan Y., Ke G., Wei Z., Gao Z.
"Bridging Machine Learning and Thermodynamics for Accurate pKa Prediction" (Uni-pKa).
JACS Au, 2024. DOI 10.1021/jacsau.4c00271 (PMC11423309). Preprint: ChemRxiv 64e8da38…, Research Square rs-3367941.
Code: https://github.com/dptech-corp/Uni-pKa (Apache-2.0). Date checked 2026-10-26.

## Fulltext status
full_text_read (PMC HTML, first ~100k chars; Supporting Information NOT read).

## Problem
Predict macro-pKa, micro-pKa and pH-dependent microstate populations of small molecules in a
thermodynamically consistent way, rather than regressing a single pKa per site.

## Method (mechanism level)
- Defines a *protonation ensemble*: a microstate is "a well-defined chemical structure of its
  protonation state" (Sec. 2.1); macrostates group microstates of equal total charge.
- A microstate enumerator expands one input structure into the full ensemble (templates adapted from MolGpKa).
- Uni-Mol (3D transformer; geometry enters as "Euclidean distances of all-atom pairs", Sec. 5.4)
  predicts a scalar standard free energy βΔ_fG per microstate.
- Micro-pKa (Eq. 1) and macro-pKa (Eq. 2) are derived from log-sum-exp over microstate free energies;
  distribution fraction vs pH from Eq. 4. Scheme 1: "macro-pKa prediction, micro-pKa prediction,
  and distribution fraction prediction".

## Data
- Pretraining: ChEMBL-derived set, ~1M molecules → >3M enumerated protonation states. Labels are
  **predicted** pKa ("predicted pKa values" / "low fidelity data", Sec. 5.5.1; Conclusions say
  "computationally predicted pKa data") — the Introduction's "empirical" wording is inconsistent.
- Finetuning: Dwar-iBonD = DataWarrior pKaInWater + curated iBonD. "the Dwar-iBonD data set grows
  over 3-fold from 8232 single data points to 27,138 enumerated microstates" (Sec. 2.3). Labels are
  **experimental macro-pKa**, not microstate populations.

## Inputs / outputs (protonation, H, pH)
- Input: SMILES → enumerated microstates (explicit protonation states, explicit H in graph) → 3D conformer(s).
- Output: per-microstate free energy ⇒ macro-pKa, micro-pKa, and pH-dependent distribution fractions.
- pH is not a model input; it enters only in the analytic post-processing (Eq. 4).

## Conformers
- Conformers are used as *data augmentation* only: "swiftly generate multiple random conformations
  allows us to incorporate data augmentation" (Sec. 5.6). **Number not stated** in the main text.
- No Boltzmann/ensemble aggregation over conformers is described; one free energy per microstate.
- No conformer-sensitivity study for Uni-pKa itself. The only conformer evidence is for the DFT
  baseline Jaguar: Table 2 footnote defines "Jaguar++" as "Jaguar with comprehensive conformational
  search, weighting 10 conformers"; on one harder system Jaguar improves "from MAE = 1.07 down to
  0.20" with exhaustive conformational modeling while Uni-pKa (no sampling) gets 0.70 (Sec. 2.5).

## Experiments & baselines
MAE on external macro-pKa sets: Novartis 0.810, SAMPL6 0.493, SAMPL7 0.554, SAMPL8 0.619/0.631
(inconsistent between Conclusions and Table 1). Micro-pKa: 10 SAMPL6 NMR-titration values,
MAE 0.592 / RMSE 0.719. Tautomer ratios (Tautobase) poor: "an MAE of 3.435 for log K".
vs Jaguar on 27 tertiary amines (Table 2): Uni-pKa MAE 0.52 vs Jaguar 0.69.

## Limitations
Conformer dependence of the microstate free energy is never modelled or measured; tautomer
equilibrium constants (which are exactly the microstate-population quantity) are predicted poorly;
pretraining labels are machine-generated, not experimental.

## Relevance to p(z,x|G,h)
Closest existing "state-level" model: it gives p(z | G, pH) through free energies, but it does NOT
give p(x | z) and explicitly does not resolve conformers — exactly the gap the project targets.
Its Tautobase log K failure is the strongest evidence that microstate populations are not solved.

## Exact supporting quotes
- "a well-defined chemical structure of its protonation state" (Sec. 2.1)
- "the Dwar-iBonD data set grows over 3-fold from 8232 single data points to 27,138 enumerated microstates" (Sec. 2.3)
- "swiftly generate multiple random conformations allows us to incorporate data augmentation" (Sec. 5.6)
- "Jaguar with comprehensive conformational search, weighting 10 conformers" (Table 2 footnote)
- "an MAE of 3.435 for log K" (Tautobase evaluation)
