# E2/E3: conformer-coupled state free energies vs sequential baseline. Results

Pre-registration: `PREREG.md` (written 2026-10-28 ~18:00 UTC / 23:30 IST; amendment 2 records a feature bug fixed before
any fit). Everything below is on measured aqueous Tautobase pairs from canonical_v1, out-of-fold.

## Runs
| job | when (UTC / IST) | what |
|---|---|---|
| 45307 | 20:57–20:58 / 02:27–02:28 | features v1. **Superseded**: the MMFF aromaticity bug left only 246/402 states with φ |
| 45568, 45570 | 21:46–21:48 / 03:16–03:18 | bug diagnosis (MMFF typing rewrites bond types in place) and fix check |
| 45572 | 21:48–21:49 / 03:18–03:19 | features v2 (fixed), pre-registered fit, E6, N1 reproduction, E1 conformer comparison |
| 45574 | 21:50–21:51 / 03:20–03:21 | exploratory robustness checks (post hoc); refactored fit reproduces 45572 byte for byte |

PREREG sha256 at fit time: `19d22f83512604de`. Outputs: `features_v2/features.csv`, `fit/` (`e2_result.json`,
`oof_predictions.csv`, `ineligible.csv`), `robust/`.

## Eligibility (identical for every arm)
- 402 states; 314 have φ (≥ 1 validated MMFF94s conformer).
- Remaining rejections:
  - `bonds_angles`: 98 states (MMFF minima outside the ±0.10/0.20 Å bounds);
  - `rt_isomeric`/`stereo`: about 14;
  - `connectivity`: 1;
  - `rt_bonds`: 0 after the fix.
- 150 of 221 measured pairs are eligible (105 scaffold groups; fold sizes 29/26/28/31/36). The 71 ineligible pairs and
  their reasons are in `fit/ineligible.csv`.
- The eligible pairs are not a random subset:
  - Their |log K| is larger (predict-0 MAE 3.75 vs 3.05 on the ineligible pairs).
  - Released Uni-pKa does worse on them (S0 MAE 3.42 vs 2.78).
  - Calibrated S1 does slightly better on them (2.67 vs 2.90).

## Pre-registered result (job 45572)
| arm | parameters | MAE | RMSE | slope (pred on measured) |
|---|---|---|---|---|
| predict 0 | 0 | 3.75 | | |
| S0 released Uni-pKa | 0 | 3.42 | 4.17 | 0.10 |
| S1 calibrated (a) | 1 | 2.64 | 3.29 | 0.47 |
| C1 a + 2D ψ | 4 | 2.48 | 3.21 | 0.49 |
| **J1 a + conformer φ** | 4 | **2.24** | 2.91 | 0.55 |
| J1 − φ1 | 3 | 2.27 | 2.98 | 0.52 |
| J1 − φ2 | 3 | 2.24 | 2.91 | 0.55 |
| J1 − φ3 | 3 | 2.63 | 3.25 | 0.50 |

Paired scaffold-cluster bootstrap (10,000 resamples), positive = J1 better:
- **Δ_S = MAE(S1) − MAE(J1) = 0.39 [0.14, 0.64]**
- **Δ_C = MAE(C1) − MAE(J1) = 0.23 [0.025, 0.49]**
- S0 − S1 = 0.78 [0.35, 1.21] (calibration alone)
- E3 ablations, increase in MAE when the term is dropped: φ1 0.03 [−0.09, 0.15]; φ2 0.00 [−0.04, 0.04];
  **φ3 0.39 [0.15, 0.63]**

**Pre-registered verdict: SUPPORTED** (both lower bounds > 0).

Coefficient signs are stable across all 5 folds:
- a = 3.7–4.5;
- θ3 = +24 to +32 log10 units per 100 Å² of Boltzmann-weighted polar SASA difference (T1 − T2);
- θ1, θ2 > 0.

The C1 TPSA coefficient is likewise +28 to +42.

## What the result does and does not show (interpretation)
- **Confirmed:** on the eligible measured pairs, adding the three conformer-ensemble descriptors to the calibrated
  Uni-pKa score lowers out-of-fold MAE more than adding three 2D descriptors does, by the pre-registered rule.
- **The gain is carried by φ3 alone** (Boltzmann-averaged polar SASA). φ1 (the conformational-entropy proxy) and φ2
  (intramolecular H-bonds) contribute nothing measurable. So this is not evidence that conformational entropy or IMHB
  coupling helps.
- **Exploratory and post hoc (not pre-registered; `robust/e2_robust.json`):**
  - **R1, canonical split.** Fit on canonical train+val eligible pairs (130), score on canonical test eligible pairs
    (19 pairs, 15 groups):
    - MAE: S0 4.35, S1 3.43, C1 3.14, J1 2.68;
    - Δ_S 0.75 [0.33, 1.31]; Δ_C 0.46 [0.15, 0.79].
    - Same direction as the pre-registered result, but n = 19.
  - **R2, fold-seed sensitivity.** With fold seeds 1–20 instead of 0:
    - Δ_S mean 0.28–0.45, lower bound > 0 in 19 of 20 seeds;
    - Δ_C mean 0.19–0.39, lower bound > 0 in only **15 of 20**;
    - both criteria met in 14 of 20.
    - **Δ_C is fragile:** the verdict depends on the fold assignment about a quarter of the time.
  - **R3, one term each.** dGU + φ3 (MAE 2.27) vs dGU + ψ3 (2D TPSA, 2.43): difference 0.16 [−0.02, 0.39], CI
    includes 0. 3D polar SASA is not shown to beat 2D TPSA head to head.
    - J1's margin over C1 may partly come from C1's ψ1/ψ2 adding noise rather than from 3D information. This is a
      hypothesis, not tested.
  - **R5, φ1 defect.** φ1 = −ln N_eff counts identical minimised conformers separately. 128 of 314 states sit at
    −ln n_valid; a rigid state with 30 identical minima gets the most "entropic" value. φ1 is therefore not a valid
    conformational-entropy proxy as implemented. It carries no weight in the fit, but it must be fixed (RMSD
    de-duplication) before any entropy claim.
- **Physics caution:**
  - a ≈ 4.2–4.8 means released Uni-pKa compresses tautomer free-energy differences about fivefold.
  - The large positive θ3 says the tautomer with more exposed polar surface is favoured in water: a solvation-like
    correction learned from 150 pairs. MMFF94s is a vacuum force field, and nothing here is a free energy. G_θ is a
    fitted score.
- **Selection caution:** the conclusion holds only for pairs where both tautomers yield a validated MMFF ensemble. The
  validator rejects many MMFF minima (see limits). Generalisation to the 71 ineligible pairs is unknown.

## Limits found while running (not fixed; thresholds kept)
1. **Validator vs force-field minima.** The 1–2/1–3 bounds come from RDKit distance geometry (UFF-based) with
   ±0.10/0.20 Å tolerance. MMFF94s minima often fall outside them:
   - 2-hydroxypyridine with the OH syn to the ring N (C–O–H ≈ 99°);
   - glycine zwitterion;
   - ibuprofen and its carboxylate;
   - phenolate.

   These states lose their whole ensemble. The cause for the anions is a hypothesis (MMFF carboxylate/phenolate
   angles against UFF-based bounds), not verified.
2. φ uses validated conformers only. That can drop the IMHB-forming rotamer and bias φ2 low.
3. **NFS attribute caching:** a job started about 30 s after an edit ran a stale copy of the source (job 45577).
   Jobs now `rsync --checksum` and log md5s of the source they ran. Affected runs: 45577 only (failed, superseded).
