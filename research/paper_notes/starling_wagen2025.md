# Starling (Rowan) — physics-informed ML macro-pKa

## Citation
Wagen C. C., "Physics-Informed Machine Learning Enables Rapid Macroscopic pKa Prediction", 2025-04-25.
ChemRxiv 68388349c1cb1ecda02ba65d (preprint, not peer reviewed); HTML full text at
https://rowansci.com/publications/macroscopic-pka-prediction. Code: Rowan platform (model "starling"); licence not checked.
Date checked: 2026-10-26 10:30 UTC (16:00 IST).

## Fulltext status
full_text_read (rowansci.com HTML via WebFetch; ChemRxiv PDF not read). Details below are from a WebFetch summary of
the page; numbers should be re-checked against the PDF before citing.

## Problem
Retrained, lightweight Uni-pKa-style model (Uni-pKa weights not freely available) for macro-pKa, logD, pI, BBB.

## Method (mechanism level)
- Microstate enumeration: SMARTS templates (RDKit), beam search in charge window −2…+2, candidates scored with AIMNet2,
  beam width 20, ≤10 iterations.
- Conformers: ETKDG + MMFF94 (RDKit); 3 conformers/molecule in pretraining, 10 at fine-tuning and inference.
- **Per-microstate energy = log-sum-exp over conformer predictions, E_micro = −log Σ_i exp(−E_i)** — i.e. a
  Boltzmann-style conformational aggregation (contrast: Uni-pKa code uses arithmetic mean).
- pH: populations over pH 0–14 (0.1 steps); macro-pKa from adjacent charge states; per-charge proton offset.
- No conformer-count ablation or conformational-sensitivity analysis reported.

## Experiments
RMSE on Novartis acid/base, SAMPL6/7/8 (e.g. SAMPL8 1.142 vs Uni-pKa 0.878); amino-acid pI; PROTAC piperazines;
brain Kp,uu classification.

## Relevance to p(z,x|G,h) / I4
Existing ML pKa model that already aggregates a (small, unweighted-by-energy-input) conformer set per microstate via
log-sum-exp. Any I2/I5 claim of "first geometry/ensemble-aware microstate free energies" is NOT available; I4 must
include Starling-style LSE aggregation as a baseline/ablation. It does not score conformer quality or state recall.
