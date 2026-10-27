# Tautomer pilot (T1 diagnostic + T2 head refit) — pre-registration
Written 2026-10-27 ~10:50 UTC (16:20 IST), before any tautomer model output exists. Claims are the judge-narrowed
versions (`scholar_runs/judge_T/verdict.md`); the intent is in `architecture_candidates/T_intent.md`.

## Claims under test
- **T1 (diagnostic):** released Uni-pKa's aqueous Tautobase errors are dominated by compression of within-macrostate
  tautomer free-energy gaps, so the slope of predicted on measured log K_T is well below 1. This is expected because an LSE macro-pKa
  loss is invariant to within-macrostate redistributions that keep each LSE fixed.
- **T2 (method):** on held-out (scaffold-grouped) measured aqueous Tautobase pairs, the Uni-pKa free-energy head fitted
  jointly to Dwar-iBonD macro-pKa and Tautobase log K_T beats the same head fitted to Tautobase alone, while macro-pKa MAE on
  SAMPL6/7/8 + Novartis rises by ≤ 0.05 against a pKa-only refit.

## Model and data (frozen)
- Released Uni-pKa weights (the same SHA-verified 5 folds as N1) through unipkainfer 0.2.0, with the N1 driver patches:
  the original seed-loop conformers (10 ETKDG+MMFF + 1 2D copy) and PKA_MEAN. CPU only.
- **Tautobase** (`ext/tautobase/Tautobase_SMIRKS.txt`, aqueous rows, numeric log K).
  - Built by `build_tautobase(kinds=MEASURED+('e','c'))`. M = measured (blank tag; 221 rows in sensitive_set_v2), E = estimated (166 numeric), C = calculated.
  - Convention: log10 K = log10([T2]/[T1]) = (G_T1 − G_T2)/ln 10.
  - **Only M pairs are evaluated.** E pairs are used only in the secondary R-joint+E training arm. C pairs are never used, except for the descriptive stage-0 table.
- **pKa:** Dwar-iBonD (8,232 transitions, 23,795 microstates) for refits. External tests are the 356 N1 transitions.
- **Leakage:** the released fold-0 head was trained on Dwar-iBonD, which contains pKa (not tautomer-ratio) data on some
  Tautobase molecules. All arms start from that same head, so the effect is common to all of them. Pairs whose tautomers both occur as
  Dwar microstates are reported separately (stage 0).

## Stage 0 — T1, zero-shot (5-fold ensemble, default A3D-mean)
Report for M, E and C separately: MAE, RMSE, mean signed error, sign accuracy on |log K| ≥ 1, the predict-zero MAE, the 2D-only MAE, MAE split by
"both tautomers present as Dwar microstates", and the OLS slope of predicted on measured log K with a 95 % scaffold-bootstrap
CI (2,000 resamples).
- **T1 falsified** (on M pairs) if any of these hold: MAE < 1.0; slope CI includes 1; or MAE < predict-zero MAE + 0.3 while slope ≥ 0.7.
- Uni-pKa reports MAE 3.435 on an unspecified Tautobase subset. Ours is compared to it for information only, not as a gate.

## Stage 1–2 — T2, head refit on cached features (fold-0 backbone frozen)
- Features: the CLS input of the free-energy head for every conformer of every microstate (Dwar-iBonD, the 5 test sets, Tautobase). The cache is
  verified by recomputing the released head's output (max |ΔG| < 1e-3, otherwise the job aborts).
- G(s) = mean over the 11 coordinate sets of head(feature), the released inference protocol.
- Arms. All are initialised from the released fold-0 head and use identical optimisation: Adam, lr 1e-4, 20 epochs × 129 steps, pKa
  batch 64, tautomer batch 16 (sampled with replacement), λ = 1 (pKa units ≡ log K units), and seeds 0/1/2.
  - R0: released head, no training.
  - R-pKa: Dwar-iBonD macro-pKa MSE only.
  - R-taut: tautomer log K MSE only. Same number of steps and tautomer batches as R-joint.
  - R-joint: pKa MSE + λ·tautomer MSE.
  - R-joint+E (secondary): R-joint plus all E pairs in training.
- Tautomer evaluation: 5-fold out-of-fold over M pairs, grouped by Murcko scaffold. Acyclic molecules form one group per pair, so duplicate
  pairs never straddle folds. Predictions are averaged over the 3 seeds. Macro-pKa evaluation averages the 15 fold×seed models
  (or the 3 seed models for R-pKa) on the 356 external transitions.
- No hyperparameter is tuned. The values above are fixed now.

## Endpoints and decision rules (fixed now)
- **Sanity gate:** R0 (fold 0 alone) external macro-pKa MAE is within 0.15 of Uni-pKa Table 1 for each of Novartis acid, Novartis base and Novartis total,
  and SAMPL6/7/8. The tolerance is wider than N1's 0.10 because this is one fold, not the 5-fold mean. If the gate fails: no conclusion.
- **Primary T2:** Δ = MAE(R-taut) − MAE(R-joint) on out-of-fold M pairs, with a scaffold-cluster bootstrap 95 % CI (10,000 resamples).
- **Co-primary (non-inferiority):** external macro-pKa MAE(R-joint) − MAE(R-pKa), with a paired bootstrap upper 95 % bound ≤ 0.05.
- **Support for T2:** Δ CI lower bound > 0 AND non-inferiority holds AND MAE(R-joint) ≤ 1.0 log K.
- **Falsification of T2:** anything else, with the failing condition(s) reported.
- **Secondary:** Δ(R-pKa − R-joint), Δ(R0 − R-joint), Δ(R-joint − R-joint+E); the slope per arm; per-stratum (stratum_v2) MAE; per-seed
  MAE; external pKa for R-taut (to show what tautomer-only fitting does to pKa).
- **Context only:** sPhysNet-Taut reports MAE 1.16 kcal/mol ≈ 0.85 log K on a different fixed 100-pair test. This is not comparable and not a gate.

## Not run here (declared stage 2, conditional)
- sPhysNet-Taut retrained and evaluated on the same folds. This is required before any paper claim ("within 0.2 log K"). It needs a
  separate approval for third-party code.
- SAMPL6 micro-pKa non-inferiority (needs the NMR micro-pKa set curated).
- Backbone fine-tuning (GPU). This is only triggered if R-taut and R-joint both stay > 1.5 log K, meaning the frozen features cannot separate tautomers.

## Expected result (stated before running)
- **T1:** slope well below 1 (~0.2–0.5) and MAE ≥ 1.5. About 75 % that T1 holds.
- **T2:** R-joint and R-taut both improve strongly over R0. Whether pKa data adds to tautomer fitting is uncertain, and I put about 55 % on T2 being falsified because the
  Δ CI includes 0. Non-inferiority is likely to hold with λ = 1.

## Compute
CPU only, 8 cores / 32 GB, gnode118, one job:
- smoke
- zero-shot (~1 h)
- feature cache (~5 h; 262k conformer passes)
- refits (~0.5–1 h)
- analysis

The released weights and code are the same as N1. This run executes them again in a new context and needs user approval before submission.

## Amendment log
(none)
