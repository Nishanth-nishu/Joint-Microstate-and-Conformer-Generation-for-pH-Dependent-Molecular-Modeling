# Rowan pKa (Wagen & Wagen 2024) — AIMNet2 + per-state conformer ensembles + CPCM-X + empirical fit

## Citation
Wagen C. C., Wagen A. M. "Efficient and Accurate pKa Prediction Enabled by Pre-Trained Machine-Learned
Interatomic Potentials." ChemRxiv 2024, DOI 10.26434/chemrxiv-2024-8489b (posted 2024-03-08).
Mirror read: https://www.rowansci.com/publications/pka-prediction. Docs: https://docs.rowansci.com/science/workflows/pka.
Successor: Rowan pKa (Wagen 2026, g-xTB), blog https://rowansci.substack.com/p/g-xtb-pka-and-website-redesign (2026-04-15).
Code: proprietary (Rowan platform). Checked 2026-10-26 UTC.

## Fulltext status
full_text_read (via rowansci.com HTML mirror, read through WebFetch summariser; ChemRxiv PDF returned 403).
Exact equations not inspected; numbers below are as reported by the mirror.

## Problem
Fast micro-pKa prediction with physics-level transferability, replacing DFT by a pretrained MLIP.

## Method (mechanism level)
- Protonation states by adding/removing H with RDKit; proton-affinity pre-screen (Methodology).
- Per state: ETKDG conformers (Careful 250 / Rapid 100 / Reckless 50 initial; Table 2), MMFF, de-duplication,
  GFN2-xTB optimisation, AIMNet2 (ωB97M-D3/def2-TZVPP) single points; "the lowest N undergo full optimization
  and frequency calculations", N = 10 by default (Methodology). Docs (2026): "AIMNet2 searches each state
  independently"; "Conformer contributions are combined into ensemble free energies at 298.15 K" following Pracht & Grimme.
- Solvation: "∆Gsolv is computed from a single-point GFN2-xTB calculation with the CPCM-X implicit water model."
- "we combine energies from each conformer to generate a single hybrid ∆G for the conjugate acid/conjugate base pair"
  (exact weighting formula not given in the mirror text → Unknown whether Boltzmann sum or min; docs say "ensemble free energies").
- ΔG → pKa via quadratic free-energy relationship plus element/valence-specific constants fitted by Levenberg–Marquardt
  on "a slightly modified version of the TR224 dataset" (Table 1).
- 2026 g-xTB version: ETKDG + ReSCoSS-style clustering, optimisation in implicit solvent, per-functional-group LFERs fitted
  to a subset of Dwar-iBond (experimental compilation).

## Data / results
SAMPL6 (24 kinase-like): MAE 0.94, RMSE 1.21. Trend studies (spiro-oxetanes, folate sulfonamides, Yudin macrocycles where
IMHB lowers pKa — "Rowan pKa was able to recapitulate the observed trend"). Ertl ring systems used for mode-vs-mode
consistency. g-xTB 2026: Rombouts BACE1 amidines MAE 1.06 → 0.59.

## Inputs/outputs
SMILES → micro-pKa per site from the neutral state ("We computed microscopic pKa values starting from each molecule's
neutral state"). Macro-pKa and tautomer ratios: not produced (SAMPL6 doubly ionised states excluded).

## Limitations (stated)
Rapid modes "may prune conformers too aggressively and often fail to find the global minimum" (fn 42); errors attributed to
"insufficient consideration of conformational effects or analysis of the wrong conformations" (Discussion). No ablation of
ensemble vs single conformer reported in the mirror text (Unknown; would be resolved by the SI/PDF).

## Relevance to I2/I5
I2: same objective + same mechanism except the conformer generator (ETKDG/CREST-screening vs torsional diffusion) and that
Rowan does micro-pKa from the neutral state rather than full macro-pKa/tautomer-ratio populations. I5: the empirical
ΔG→pKa map is a fitted correction but it is per-element/functional-group, not geometry-conditioned.

## Exact supporting quotes
- "we combine energies from each conformer to generate a single hybrid ∆G for the conjugate acid/conjugate base pair" (Methodology)
- "AIMNet2 searches each state independently" (Rowan docs, pKa workflow page)
- "the rapid workflows here may prune conformers too aggressively and often fail to find the global minimum" (fn 42)
