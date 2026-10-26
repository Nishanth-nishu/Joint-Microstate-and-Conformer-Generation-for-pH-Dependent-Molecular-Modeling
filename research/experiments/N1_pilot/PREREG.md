# N1 pilot — pre-registration (written 2026-10-26, before any model output exists)

Claim under test (verbatim, `architecture_candidates/N1_intent.md`): for thermodynamically consistent microstate pKa
models (Uni-pKa class), whether the input carries 3D conformers — and how they are generated and aggregated — changes
experimental pKa error by more than resampling noise only on conformationally sensitive molecules.

## Model and data (frozen)
- Model: released Uni-pKa weights via `bindresearch/pkasso_unipka` @ e28a591 (Apache-2.0), 5 folds, SHA-256 verified
  against unipkainfer 0.2.0's registry. Inference code: unipkainfer 0.2.0 (Apache-2.0/MIT). No training, no tuning.
- Pilot test sets: Uni-pKa's own external test sets (`ext/unipka_orig/dataset/`): SAMPL6 (31), SAMPL7 (20), SAMPL8 (25),
  Novartis acid (112), Novartis base (168) macro-pKa transitions given as acid>>base microstate sets. These were NOT in
  Uni-pKa fine-tuning (per the paper), so arms are compared out of sample.
- Strata: **`stratum_v2`** from `jmcg/data/build_sensitive_set.py::flags_v2` (6/7-membered flexible IMHB to a strong
  partner, flexible acid–base pair 3–8 bonds, or flexible site pair 4–8 bonds), computed on the uncharged parent of each
  transition's first acid microstate. v1 (`flags`) was replaced on 2026-10-26 BEFORE any model output because it put ~75 %
  of IUPAC entries in 'sensitive' (non-discriminative); v1 columns are still written for transparency.

## Arms (all from ONE forward pass: 10 ETKDG+MMFF conformers + 1 2D-coordinate copy per microstate, seed 0)
- A3D-mean (default Uni-pKa/unipkainfer): G = mean over the 11 coordinate sets.
- A2D: G = model output on the 2D-coordinate copy only.
- A3D-only-mean: mean over the 10 3D conformers (2D copy excluded).
- A3D-lse: G = −ln( (1/10) Σ_c exp(−G_c) ) over the 10 3D conformers (Boltzmann-style, β-scaled units).
- A3D-min: min over the 10 3D conformers.
- A3D-single: conformer 0 only (resampling-noise reference).
pKa per transition = (LSE_A(−G) − LSE_B(−G)) / ln 10; fold-averaged over 5 folds.

## Endpoints
- Primary: Δ_s = MAE(A2D) − MAE(A3D-only-mean) within the **sensitive** stratum; molecule-level paired bootstrap
  (10,000 resamples), 95 % CI. And difference-in-differences DiD = Δ_s − Δ_c (control stratum).
- Secondary: MAE of A3D-lse and A3D-min vs A3D-only-mean (sensitive stratum); Spearman ρ between per-transition
  conformer spread of predicted pKa and |error|; per-dataset MAE of A3D-mean.

## Gates and decision rules (fixed now)
- **Sanity gate:** A3D-mean MAE per set must be within 0.10 of Uni-pKa's reported values
  (SAMPL6 0.493, SAMPL7 0.554, SAMPL8 0.619–0.631, Novartis 0.810, `paper_notes/unipka_2024.md`). If not, the pipeline
  is not faithful → stop and debug; no scientific conclusion.
- **Support for N1:** Δ_s > 0 with CI excluding 0 AND DiD CI excluding 0.
- **Falsification of N1 (on this pilot):** Δ_s CI includes 0, or DiD CI includes 0 (3D helps everywhere equally, so
  the effect isn't sensitivity-specific).
- **Underpowered:** if the sensitive stratum has < 30 transitions, report as underpowered and extend to the
  leakage-clean IUPAC subset (`data_manifests/sensitive_set_v1`, `in_unipka_finetune == 0`) as pre-declared stage 2.
- Conformer-source arm (torsional-diffusion conformers) is stage 2; not run here.

## Expected result (stated before running)
Effect small: Δ_s ≈ 0.05–0.2 pKa units, likely with a CI touching 0 at pilot n. I put roughly 60 % on N1 being
falsified or underpowered at this stage. A3D-lse/min likely no better than mean (cf. Starling vs Uni-pKa).

## Compute
CPU only (torch 2.5.1+cpu), 8 cores / 32 GB (node had 2 free CPUs at submission time), gnode118. No GPU requested.

## Amendment log
- 2026-10-26 ~10:52 UTC (before any model output): primary stratum switched v1→v2 (see Strata). Analysis default `--stratum_col stratum_v2`.
- 2026-10-26 22:24 UTC — **Run 1 (job 41159) failed the pre-registered sanity gate** (A3D-mean MAE ≈ 6.4 vs reported
  0.49–0.81) → no scientific conclusion drawn, per the gate. Root causes (diagnosed jobs 41187/41189/41192):
  (a) the pKa conversion omitted Uni-pKa's target-normalisation mean 6.5049 (`unimol/tasks/unimol_pka.py:184`): mean
  signed error was −6.47 ± 0.87; (b) unipkainfer 0.2.0's free-energy path generates conformers with
  `EmbedMultipleConfs(randomSeed=0)`, which returns 10 identical conformers (RDKit 2024.03 and 2026.09; seeds 1/42/−1 are
  diverse), unlike the original Uni-pKa loop `EmbedMolecule(randomSeed=seed) for seed in range(10)`; so run 1's "3D
  ensemble" arms were a single conformer. Fixes (driver only, before run 2): add PKA_MEAN; replace the conformer function
  with a verbatim port of the original Uni-pKa protocol. Arms, strata, endpoints and decision rules are unchanged.
  Run 1 outputs kept in `run1_failed_sanity_gate/` and must not be interpreted.
