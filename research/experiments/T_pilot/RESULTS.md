# Tautomer pilot: result (job 43158, finished 2026-10-27 21:06:26 UTC / 2026-10-28 02:36:26 IST)

Pre-registration: `PREREG.md`. Numbers: `run1/analysis.json` (identical to the STAGE 6 JSON in `logs/jmcg_tpilot_s4_43158.log`).

- **Sanity gate:** passed. Fold-0 macro-pKa MAE was within tolerance of Uni-pKa Table 1 on Novartis acid and base,
  SAMPL6, SAMPL7 and SAMPL8.
- **T1 (zero-shot compression): supported.** On 221 measured aqueous Tautobase pairs, MAE was 3.21 log K. The
  regression slope of predicted on measured was 0.10 [0.08, 0.12], and predicting zero gives an MAE of 3.53.
- **T2 (joint head refit beats tautomer-only): falsified.** All three support conditions failed:
  - Δ = MAE(R-taut) − MAE(R-joint) = −0.007 log K [−0.118, 0.112]; the CI includes 0.
  - Macro-pKa non-inferiority failed: R-joint − R-pKa = +0.22 [+0.16, +0.28] on the external sets.
  - MAE(R-joint) = 1.91 > 1.0.
- **Exploratory, not pre-registered as primary:** adding experimental tautomer energies (R-joint+E) gave 1.78, a
  difference of +0.13 [+0.04, +0.22] relative to R-joint.
