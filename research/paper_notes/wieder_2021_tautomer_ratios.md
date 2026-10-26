# Wieder, Fass, Chodera 2021 — tautomer ratios in solution (KEY evidence on conformer coupling)

## Citation
Wieder M., Fass J., Chodera J. D. "Fitting quantum machine learning potentials to experimental free
energy data: predicting tautomer ratios in solution." Chem. Sci. 12:11364 (2021).
DOI 10.1039/D1SC01185E (PMC8409483). Code: github.com/choderalab/neutromeratio. Checked 2026-10-26.

## Fulltext status
full_text_read (PMC, first ~100k chars).

## Data
Experimental aqueous tautomer ratios from **Tautobase** (DataWarrior file, Wahl/Sander, from
P. W. Kenny's tautomer codex): 1680 deposited pairs → 476 usable + 2 SAMPL2 = 478; 468 pairs used
for QM, 354 for QML (212/71/71 train/val/test). Labels are experimental ΔG = −RT ln K. Calculations
at 298.15 K (QM) / 300 K (MD). Experimental temperature not stated in the paper.

## Method
Two routes for the same ΔG between tautomers: (a) RRHO thermochemistry on minima — RDKit conformers
(20 generated, pruned to 10 or 5 by RMSD); (b) alchemical relative free energies with full
configurational sampling (single topology, 11 λ, 200 ps/state, 5 repeats, MBAR).

## The magnitude result (directly relevant to conformation-dependence)
- RRHO vs full configurational sampling: "The mean absolute deviation of 0.9 kcal mol−1 should not
  be underestimated" and "the RRHO approximation introduces a mean absolute error of ∼1 kcal mol−1"
  (Discussion). ≈0.7 log-units of population.
- Multiple conformers vs the global minimum alone: "Much more relevant than including multiple
  conformations is locating the global minimum conformation" — only 12 of the pairs changed by
  >1 kcal/mol when multiple minima were used.

## Accuracy
Best QM (B3LYP/aug-cc-pVTZ//B3LYP/6-31G(d)/SMD, 468 pairs) RMSE 3.1 kcal/mol, wrong dominant
tautomer for 13% of pairs. ANI-1ccx vacuum alchemical RMSE 6.4; after fitting to experiment,
2.8 kcal/mol on the 71-pair test set.

## Relevance to p(z,x|G,h)
Two-edged and important:
(+) Microstate (tautomer) free energies in water are badly predicted even by DFT (RMSE ~3 kcal/mol
    ≈ 2 pKa units), so task C is genuinely unsolved and worth modelling.
(−) The paper's own ablation says *which* conformers you include matters far less than finding the
    global minimum (only ~2.5% of pairs shift >1 kcal/mol), and the full-sampling-vs-RRHO gap is
    ~1 kcal/mol. That is the strongest published evidence that explicit joint conformer-level
    modelling of microstate populations buys ≲1 kcal/mol for typical drug-like molecules.

## Exact supporting quotes
- "The mean absolute deviation of 0.9 kcal mol−1 should not be underestimated." (RRHO vs alchemical)
- "the RRHO approximation introduces a mean absolute error of ∼1 kcal mol−1 in the calculation" (Discussion)
- "Much more relevant than including multiple conformations is locating the global minimum conformation"
