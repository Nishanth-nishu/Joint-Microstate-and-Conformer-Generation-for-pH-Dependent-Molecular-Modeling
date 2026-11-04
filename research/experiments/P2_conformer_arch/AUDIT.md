# P2: audit of the "modified conformer architecture" before any training (2026-10-29)

Sources:
- read-only inspection jobs 46545, 46552, 46565 and 46569 (gnode118, 2026-10-29 10:21–10:30 UTC / 15:51–16:00 IST);
- the TorDiff-NAX handoff `prior_work/TORDIFF_NAX_HANDOFF.md`, `round5/specs/EXP-1.md` rev 2/3, `round5/IMPLEMENTATION.md`
  and `EXPERIMENTS.md`, all in the jmcg repo;
- the flexitors package (`$P/flexitors`, identical to `/scratch/nishanth.r/tordiff/flexitors` at commit
  `1c4b7a54…` except `__pycache__`).

## 1. What actually differs from the original Torsional Diffusion model

**No modified score-model architecture exists in code.**
- `flexitors/models/__init__.py` and `flexitors/training/__init__.py` are docstring-only stubs ("Empty in milestone
  1 by design").
- The EXP-1 geometry head (V2: a per-atom `1o` L-velocity head on TD's `TensorProductScoreModel`, a second clock
  through the λ embedding, a P_W-projected endpoint loss, and DAgger) is specified in EXP-1 §2.2/§10 but was never
  built. The handoff lists it as "Milestone 3 — not started".
- So there is no modified checkpoint to fine-tune or audit for compatibility.

What does exist:

| Component | Original TD | What exists beyond it | Trained checkpoints |
|---|---|---|---|
| Score network | `TensorProductScoreModel`: e3nn, 4 layers, ns 32, nv 8, `0o` torsion head | round-2 S4: `--lambda_embed_dim 32` concatenates a λ embedding to node and edge inputs (`score_model.py:70,76,194-198`); `--no_parity` hook. Otherwise unchanged | S4 × 3 seeds (`workdir/qm9_S4_lamcond_e100_s{0,1,2}`) |
| Noise / score parameterisation | wrapped normal on T^m, σ ∈ [0.0314, 3.14] log-uniform, DSM loss divided by E‖∇log p‖² | unchanged in every trained variant | — |
| Torsion representation | rotatable = cut leaves > 1 atom on both sides; one pseudoscalar per bond; L frozen from ETKDG | unchanged. Inference-time only: flexitors stereo guard (K1/K2, Exp 3) and the jmcg E/Z lock remove or freeze bonds, with no retraining | — |
| Training data | GEOM-QM9 standardized pickles (conformer matching) | recipe variants: rematch, raw/true L (B1), DE-only (B2), MMFF-matched (B3/S5), heavy-atom matching (B5), jitter (B6), 50/50 mix (S3), paired λ blends (S4) | 3 seeds each, `workdir/qm9_*_e100_s*` |
| Loss | TD DSM | unchanged everywhere | — |
| Checkpoint compatibility | `qm9_default` (sha256 `6a9a6ee9…`) | S4 checkpoints are not strictly loadable into the default architecture (wider first linear layers). All other variants are loadable | — |

Earlier round-2 findings, quoted from `EXPERIMENTS.md`, not re-run here:
- No training-recipe variant beat the control in real use; only better geometry at inference (MMFF seeds) did.
- S4 on RDKit L: +0.0033 Å [−0.0009, +0.0078], not significant.

## 2. Reference for any comparison: frozen TD in E1 (already measured, not re-run)

- **Molecules:** 200 GEOM-QM9 test molecules (`research/experiments/M3_seq/e1_confs/subset.csv`, seed 0 sample of
  TD's `test_smiles.csv`), 183 common to all arms after model failures and RMSD NaN.
- **Evaluation policy:** budget 2 × n_true, scored with TD's `evaluate_confs.py`, paired bootstrap over molecules
  (`jmcg/condgen/e1_confs.py compare`).

Results on the 183 molecules (job 45308, `compare_v2.json`):

| arm | MAT-R (Å) | COV-R@0.5 | COV-R@0.25 |
|---|---|---|---|
| TD `qm9_default`, seed 0 / 1 | 0.170 / 0.164 | 94.8 / 96.2 | 82.7 / 81.3 |
| ETKDGv3, seed 0 / 1 | 0.230 / 0.216 | 89.1 / 91.0 | 69.3 / 68.6 |

This subset is drawn from the **test** split, so it may be read once for a finished model. It must not be used to
choose checkpoints, σ_b, w or any other setting; EXP-1 §10.3 assigns those to the validation split.

## 3. What the binding spec requires before any training (EXP-1 §5.2, §10.3)

1. **M0, oracle-velocity integration test (CPU, no learning).** Does the sampler's projected, torsion-restored,
   constraint-guarded L update follow an oracle field to the pure-L target?
   - **Kill** if recovery < 50 % or > 5 % of primary pairs end with a skipped step or a flipped centre.
   - Run on validation pairs only.
2. **1-epoch training smoke (≈ 0.12–0.18 GPU-h).** Measures the real overhead of the geometry head and the
   projector before the Gate-0 ask.
3. **M1 = Gate 0 (≈ 10.5 GPU-h, ceiling 16).** Geometry head only (w_τ = 0), 3 σ_b × 30 epochs, decided on
   validation (rigid stratum AMR-R ≤ 0.10 Å; ≥ 4-ring > 10° share ≤ 20 %).
   - M1 does not produce a full conformer generator: torsions stay at the seed in the read-out. So no MAT/COV
     comparison with TD and ETKDG on the E1 molecules is meaningful before M2 (full V2 with DAgger, ≈ 35–53 GPU-h
     including controls).
4. The handoff also records unresolved milestone-2 review blockers in the Exp 3/Exp 5 sampler code
   (`round5/review_P1_m2.md`, `review_P2_m2.md`). They do not touch the Exp 1 geometry path built here.

## 4. What was done in this step
- **M0 implemented and run on the validation cache: PASS by the pre-declared rule** (`research/experiments/EXP1_M0/RESULTS.md`).
  - Primary population (30,420 pairs): recovery 99.6 %; 4.74 % of pairs end with a skipped step (limit 5 %); 0 flipped centres.
  - Code: `flexitors/geometry/projector.py`, `flexitors/geometry/l_step.py`, `experiments/exp1/m0_oracle.py`;
    tests in `flexitors/tests/test_projector.py`.
- **Not done, and why:** no model code, training configuration, checkpoint path, training smoke or E1 comparison.
  - The modified model (the V2 geometry head, transform, projected loss and batched torch projector) does not exist
    yet. Building it is the next milestone.
  - Its first GPU spend (the 1-epoch smoke, then Gate 0 at ≈ 10.5–16 GPU-h) goes against the GPU budget rule in the
    handoff ("ask the user before GPU-heavy runs").
  - A like-for-like E1 comparison (MAT/COV against TD and ETKDG on the 200 E1 molecules) needs M2: full V2 with a
    torsion head and DAgger, ≈ 35–53 GPU-h including controls.
- **Before M1:** record the reason for every skipped step (bond window / clash / chirality) per step. M0's 2.7 %
  same-stereo skips are not attributed.
