# EXP-1 V2 training smoke: build check, CPU smoke, overfit diagnostic (exploratory)

**Verdict: smoke FAIL on S2 (learning). Plumbing passes (S1, S3, S4, S5).** By the brief's rule ("if the smoke passes,
run a controlled training experiment"), no controlled V2-vs-A3 training or E1 evaluation was run. Gate 0 was not started.
Nothing is committed.

## What was built (staging, not committed)
- `flexitors/geometry/projector_torch.py`: batched, differentiable P_W (padded C, float64 pinv); equals the numpy projector (test).
- `flexitors/training/transforms.py`: `PairedJointTransform` (§10.1): same Δτ on both endpoints, pure-L target x_G*′,
  bridge x_t, chirality resample (≤ 5) or mask, pair_ok-false samples get the τ loss only.
- `flexitors/models/flexitors_score.py`: TD `TensorProductScoreModel` + `o3.Linear → 1x1o` geometry head on the last atom
  layer (forward hook); geometry clock via TD's λ embedding (32). Golden test: torsion output bit-equal to TD's with the
  head present; parity (1o negates under x → −x) and rotation equivariance tested.
- `flexitors/training/losses.py` (TD torsion loss; projected endpoint loss / rank), `flexitors/training/loop.py`
  (own `--log_dir`, atomic checkpoints with RNG states, reproducible validation), `flexitors/sampling/l_readout.py`
  (τ-frozen 20-step read-out with torsion restoration and guard), `flexitors/tests/test_v2.py`, `experiments/exp1/v2_smoke.py`.

## Compute
GPU job 46725 never started: plafnet2 account GPU-minute quota exhausted (`AssocGrpGRESMinutes`, GrpTRESMins gres/gpu=345600).
Everything below ran on CPU on gnode118. GPU time per epoch is unmeasured.

## Jobs
| job | what | result |
|---|---|---|
| 46738 | CPU smoke | cancelled (`nvidia-smi` exit under `set -e`) |
| 46739 | CPU smoke | test failure: CLI flags declared in loop.py, read in flexitors_score.py (moved to `add_model_args`) |
| 46744 | tests, CPU resume check, 400-iteration training (12:39:58–14:00:56 UTC / 18:09:58–19:30:56 IST) | tests 110 passed / 3 skipped; RESUME_OK; training finished; checks crashed reading `model_parameters.yml` (`torch.__version__` dumped as a python object; fixed in loop.py, reader uses `yaml.Loader`) |
| 46848 | smoke checks on the 46744 run (no retraining) | report below |
| 46855 | overfit diagnostic, 32 training molecules (14:04:57–14:09:36 UTC / 19:34:57–19:39:36 IST) | below |

Configuration (46744): w_τ 1, w_geo 1, σ_b 0.15, seed 0, batch 32, lr 1e-3, 1 epoch capped at 400 iterations,
validation capped at 60 batches; round-2 paired cache (`/scratch/nishanth.r/tordiff/data/QM9/cache_paired.{train,val}`, read only).
Checkpoints: `$P/workdir/exp1/v2_smoke_s0/{best,last}_model.pt` (sha256 in `checkpoints.sha256`).

## Smoke criteria (declared in `v2_smoke.py` before the run)
| criterion | result | detail |
|---|---|---|
| S1 finite | PASS | all logged and validation losses finite |
| S2 learning | **FAIL** | val geo 0.01302 vs Δ̂ = 0 base 0.01288 (ratio 1.010; required < 0.9). Train first-5 0.196 → last-5 0.0108 passes, but only because iteration 0 starts at 117× base |
| S3 gradients | PASS | all groups finite and > 0; l_head grad 1.05 → ~0.005 |
| S4 checkpoint | PASS | strict loads of best/last; bit-identical round trip; 400 optimizer steps; CPU resume continued at epoch 1 |
| S5 determinism | PASS | CPU read-out repeat max diff 0.0 (no GPU available: the "GPU" repeat also ran on CPU) |

Torsion loss after 400 iterations: val 1.011 vs base 1.003 (no learning yet either; TD itself learns slowly early on).

Read-out (informative, 200 validation pairs, seed = cache-matched RDKit conformer): median heavy-atom RMSD to x_G*
0.106 Å (seed) → 0.127 Å (read-out); read-out closer in 15 %; ring dihedral > 10° 56 % → 54 %; ≥ 1 guard skip in 72.5 %.
The untrained-in-effect head makes geometry worse; this is expected given S2.

Training trajectory (train geo / base, every 20 iterations): 116.8, 1.62, 1.30, 1.24, 1.13, 1.25, 1.28, 1.11, 1.09, 1.04,
1.13, 1.04, 1.04, 1.04, 1.06, 1.09, 1.06, 1.01, 0.98, 1.01. The head decays from a large random-init output towards
Δ̂ ≈ 0; it does not get below the zero-prediction base.

## Overfit diagnostic (job 46855; rule declared in the sbatch header before the run)
32 training molecules, batch 8, 100 epochs = 400 iterations, CPU. Pass iff mean train geo/base over the last 10 epochs < 0.9.

| arm | ratio first 10 epochs | ratio last 10 | verdict | td/base last 10 |
|---|---|---|---|---|
| w_τ 0, w_geo 1 | 3.33 | 0.856 | pass (narrowly) | 1.05 (untrained) |
| w_τ 1, w_geo 1 | 3.45 | 0.913 | fail | 0.998 |

## Interpretation (hypotheses, not confirmed)
1. The geometry objective can be fitted (w_τ 0 overfit gets below base), but slowly and only to ~0.86 × base on 32 molecules.
2. With the torsion loss on, the geometry loss is ~80× smaller (0.013 vs 1.0) and its gradient share of the shared trunk is
   small. A loss-scale imbalance is the leading hypothesis for S2. Untested: normalising L_geo by its base, or a larger w_geo.
3. Part of the base is irreducible bridge noise (σ_b = 0.15 Å vs RMS signal ≈ 0.11 Å per coordinate), so the attainable
   floor is above 0. The size of that floor is not measured.
4. Timing anomaly: 46744 took 11.7 s/iteration (batch 32, 6 loader workers, OMP 10 on 16 CPUs); 46855 took ~0.3 s/iteration
   (batch 8, no workers). Likely CPU oversubscription in 46744; not investigated. CPU training may be more feasible than
   46744 suggested.

## Not done / blockers
- Controlled training (V2 vs A3 ablation, w_geo = 0) and E1 evaluation vs frozen TD and ETKDGv3: not run, smoke failed.
  A joint sampler (torsion SDE + geometry step) for E1 is also not built.
- GPU quota exhausted for plafnet2.
- Per-step guard skip reasons are still not logged (TODO before Gate 0).
- No S2 threshold was changed; any fix (loss normalisation, w_geo) needs a new pre-declared smoke.

## Re-smoke with L_geo normalised by its base (job 46998, 16:43–17:11 UTC / 22:13–22:41 IST): FAIL on S2
Declared before the run (sbatch header): one objective change, `--geo_norm base` (optimise L_geo / its Δ̂ = 0 base,
base detached, per batch); all else as 46744; S1–S5 unchanged, S2 on the raw losses. Infrastructure: loader workers
use 1 torch thread. Tests: 111 passed, 3 skipped (new `test_geo_norm_base_objective`).

| | 46744 (raw) | 46998 (geo_norm base) |
|---|---|---|
| val geo / base | 0.01302 / 0.01288 = 1.010 | 0.01171 / 0.01183 = **0.989** (S2 needs < 0.9) |
| val td / base | 1.011 / 1.003 | 0.976 / 0.973 |
| l_head grad, iterations 20–380 | 0.066 → ~0.005 | 0.15–1.24 (stays O(1)) |
| read-out closer to x_G* (200 val pairs) | 15 % | 20.5 % |
| median read-out RMSD to x_G* (seed 0.106 Å) | 0.127 Å | 0.128 Å |
| read-outs with ≥ 1 guard skip | 72.5 % | 4 % |
| ring dihedral > 10° (57 ring pairs; seed 56 %) | 54 % | 49 % |
| s / iteration (CPU) | 11.7 | 2.9 |
| overfit diagnostic, w_τ 1, last-10 ratio | 0.913 | 0.952 |

S1, S3, S4, S5 PASS again; overall FAIL.

Train geo/base every 20 iterations: 116.8, 1.74, 1.32, 1.20, 1.11, 1.05, 1.14, 1.12, 1.13, 1.03, 1.07, 1.01, 1.01, 0.97,
1.07, 1.00, 1.08, 0.99, 1.00, 0.98.

**Reading.** Normalisation fixed the gradient share (head gradient stays O(1) instead of decaying to ~0.005) and the
oversubscription fix made CPU iterations 4× faster, but the geometry loss still sits at the zero-prediction base. So the
loss-scale hypothesis is **not supported** as the main cause. Remaining hypotheses (untested):
1. Random-init head output is ~117× base; most of 400 iterations go to unlearning it. Zero-initialising `l_head`
   would remove that phase.
2. Much of the per-batch target is unpredictable bridge noise (σ_b = 0.15 Å vs ~0.11 Å RMS signal); a smaller σ_b
   (spec sweep includes 0.05) would raise the learnable share.
3. 400 iterations (~0.1 epoch) is too short; TD's torsion loss is also still at base.

## Re-smoke 3 with a zero-initialised geometry head (job 47102, 17:22–17:48 UTC / 22:52–23:18 IST): **PASS**
Declared before the run (sbatch header): one change vs 46998, `--l_head_init zero`; criteria unchanged. Job 47095
(first submission) failed at the tests (dead-code check: flag read via `getattr`), so nothing trained; fixed and resubmitted.
Tests: 112 passed, 3 skipped.

| | 46998 | 47102 |
|---|---|---|
| val geo / base (S2 < 0.9) | 0.989 | **0.703** (0.00831 / 0.01183) |
| read-out closer to x_G* (200 val pairs) | 20.5 % | **61 %** |
| median read-out RMSD to x_G* (seed 0.106 Å) | 0.128 Å | **0.0996 Å** |
| ring dihedral > 10° (seed 56 %) | 49 % | **25 %** |
| read-outs with ≥ 1 guard skip | 4 % | 13.5 % |

S1, S3, S4, S5 PASS → **smoke PASS**. Train geo/base every 20 iterations: 1.00, 0.99, 0.99, 0.99, 0.99, 0.99, 0.995,
0.98, 0.97, 0.975, 0.93, 0.92, 0.92, 0.85, 0.82, 0.79, 0.85, 0.77, 0.75, 0.71 (still falling at iteration 400).
Reading: with the default init, the head's ~117× base initial output had to be unlearned first, and the run ended
before learning started. Three configurations were tried in sequence (one change each); the gate criteria never changed.
The secondary 2,000-iteration run in the same job was cancelled after the gate passed (superseded by EXP1-V2-T1).
The focused diagnosis (`experiments/exp1/v2_diag.py`) was written but not run, as S2 passed.
