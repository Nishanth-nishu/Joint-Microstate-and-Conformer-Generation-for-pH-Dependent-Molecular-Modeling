# EXP1-V2-T1 results (design pre-declared in `DESIGN.md`; exploratory, CPU, reduced budget)

## Jobs
| job | what | start → end (UTC / IST) | state |
|---|---|---|---|
| 47246 | T1-proj training, 10 epochs, CPU 14 cores | 2026-10-29 18:56 / 10-30 00:26 → 10-30 21:47 / 10-31 03:17 | COMPLETED (26.9 h) |
| 47247 | T1-noproj training, same | 2026-10-29 18:56 / 10-30 00:26 → 10-31 00:17 / 10-31 05:47 | COMPLETED (29.4 h) |
| 47249 | E1 read-outs + evaluate_confs + paired bootstrap | 10-31 00:17 / 05:47 → 00:50 / 06:20 | all outputs written; Slurm marks OUT_OF_MEMORY (see Failures) |
| 49875 | leakage check, E1 test molecules vs paired caches | 10-31 01:09 / 06:39 → 01:11 / 06:41 | COMPLETED: 1 shared graph (see Failures) |

Checkpoints (validation-selected `best_model.pt`, epoch 9 for both arms): proj `09a0d876…82d0`, noproj `c647742f…9c36`
(full hashes in `T1_*_s0/checkpoints.sha256`).

## Training curves (geometry loss / its Δ̂ = 0 base; val = first 100 validation batches, fixed noise)
| epoch | proj train | proj val | noproj train | noproj val |
|---|---|---|---|---|
| 0 | 0.450 | 0.343 | 0.499 | 0.397 |
| 2 | 0.289 | 0.293 | 0.339 | 0.344 |
| 4 | 0.261 | 0.252 | 0.307 | 0.300 |
| 6 | 0.242 | 0.238 | 0.285 | 0.293 |
| 8 | 0.233 | 0.232 | 0.277 | 0.281 |
| 9 | 0.232 | 0.232 | 0.278 | 0.274 |

Still decreasing at epoch 9 (budget-limited, not converged); lr never reduced; no non-finite loss. The two arms
optimise differently defined losses (projected/rank vs plain/3N), so their loss values are not comparable to each other.

## E1 (200-molecule GEOM-QM9 test subset, 184 molecules common to all arms, 2 × n_true, TD evaluate_confs)
Absolute (seed 0 / seed 1), Å or %:

| arm | MAT-R | MAT-P | COV-R@0.25 | COV-P@0.25 | COV-R@0.5 | COV-P@0.5 |
|---|---|---|---|---|---|---|
| TD | 0.170 / 0.165 | 0.215 / 0.206 | 82.8 / 80.9 | 71.4 / 71.4 | 94.9 / 96.2 | 91.0 / 92.9 |
| TD + L(proj) | **0.115 / 0.110** | **0.149 / 0.138** | 87.8 / 87.9 | 81.4 / 82.5 | 95.7 / 96.8 | 93.8 / 94.7 |
| TD + L(noproj) | 0.116 / 0.111 | 0.150 / 0.140 | 87.8 / 88.3 | 81.2 / 82.4 | 95.7 / 96.8 | 93.6 / 94.7 |
| TD + MMFF94s (control) | 0.163 / 0.152 | 0.182 / 0.170 | 82.5 / 83.0 | 78.8 / 79.7 | 94.2 / 96.3 | 92.9 / 93.8 |
| ETKDGv3 | 0.229 / 0.217 | 0.227 / 0.221 | 69.5 / 68.2 | 67.9 / 67.6 | 89.2 / 91.1 | 90.0 / 91.8 |

Paired differences, mean [95 % bootstrap CI over molecules]:

| contrast | seed | MAT-R | MAT-P | COV-R@0.25 | COV-R@0.5 |
|---|---|---|---|---|---|
| **L(proj) − TD (primary)** | 0 | **−0.0547 [−0.0628, −0.0473]** | −0.0659 [−0.0761, −0.0568] | +4.99 [2.76, 7.63] | +0.86 [−0.04, 2.27] |
| L(proj) − TD | 1 | −0.0557 [−0.0638, −0.0486] | −0.0679 [−0.0777, −0.0589] | +7.03 [3.95, 10.69] | +0.53 [0.08, 1.22] |
| L(proj) − L(noproj) | 0 | −0.0009 [−0.0020, 0.0002] | −0.0016 [−0.0026, −0.0006] | −0.02 [−0.37, 0.35] | +0.03 [−0.03, 0.10] |
| L(proj) − L(noproj) | 1 | −0.0012 [−0.0022, −0.0003] | −0.0017 [−0.0027, −0.0007] | −0.36 [−1.04, 0.09] | −0.03 [−0.08, 0.00] |
| MMFF − L(noproj) | 0 | +0.0470 [0.0369, 0.0573] | +0.0312 [0.0217, 0.0405] | −5.28 [−7.95, −2.93] | −1.49 [−3.29, −0.18] |
| MMFF − L(noproj) | 1 | +0.0415 [0.0331, 0.0503] | +0.0302 [0.0217, 0.0388] | −5.30 [−7.93, −2.96] | −0.51 [−1.70, 0.17] |
| MMFF − TD | 0 | −0.0068 [−0.0184, 0.0043] | −0.0332 [−0.0444, −0.0221] | −0.27 [−3.66, 3.26] | −0.66 [−2.13, 0.49] |
| ETKDG − TD | 0 | +0.0593 [0.0421, 0.0768] | +0.0125 [−0.0012, 0.0274] | −13.27 [−17.19, −9.48] | −5.67 [−8.54, −2.91] |

Full tables: `e1/compare_seed{0,1}.json` (ref TD), `e1/ablation_seed{0,1}.json` (ref L(noproj)).

## Verdicts under the pre-declared rules
- **Improves on TD: yes.** Primary MAT-R CI excludes 0 (−0.055 Å, ≈ −32 %), seed 1 same sign and size.
- **Projector matters: not shown.** proj − noproj on MAT-R, seed 0, CI includes 0. Seed 1 excludes 0, but the effect is
  ≈ 1 mÅ against a 55 mÅ gain over TD: the learned correction, not the projector, carries the result.
- **Beats the MMFF control: yes** (by ≈ 0.04–0.05 Å MAT-R, CI excludes 0, both seeds); MMFF alone barely moves MAT-R.

## Read-out diagnostics (per seed, 5,390 / 5,386 conformers, 185 molecules)
| arm | W-dihedral drift (max) | conformers with ≥ 1 guard skip | mean halvings | chirality rejections | final geometry / chirality fail | median heavy-atom displacement | CPU s / conformer | errors |
|---|---|---|---|---|---|---|---|---|
| L(proj) | 2e-13° | 4.8 % / 5.1 % | 1.07 / 1.13 | 0 | 1 / 0 ; 0 / 0 | 0.101 Å | 1.35 | 46 |
| L(noproj) | 2e-13° | 4.5 % / 4.4 % | 0.85 / 0.90 | 0 | 1 / 0 ; 0 / 0 | 0.114 Å | 1.33 | 46 |
| MMFF | 1e-13° | – | – | – | 107 / 16 ; 105 / 16 | 0.33 Å | 0.02 | 0 |

Torsions are unchanged in every arm, so all differences come from local geometry (bonds, angles, ring puckers).

## Failures and deviations
- **Read-out errors (46 conformers, 12 molecules, every L arm):** all molecules with no rotatable bond; TD's
  `score_model.forward` fails on an empty torsion-edge set (`score_model.py:164`). Those conformers were kept as TD's,
  as the script specifies, so the reported gain is if anything slightly understated. Fix before the next run.
- **47249 OUT_OF_MEMORY (18 oom_kill events, MaxRSS 12 GB of 40 GB):** every step finished, no `eval FAILED`, every
  evaluation log has the same failure counts as TD's own evaluation (15 model failures, 1 additional, 1 RMSD-NaN drop),
  and all arms pair on the same 184 molecules. The killed processes were not identified; outputs are complete.
- Earlier: 47158/47159 TIMEOUT at the 1 h partition default (no epoch finished; logs in `failed_47158_47159/`).

- **Leakage check (49875):** comparing canonical SMILES with stereochemistry removed, 1 of the 200 test molecules
  (`N#CC12C3C4CN1C2C43`) also occurs in the training cache (likely a stereoisomer or a duplicate graph in TD's split);
  none in the validation cache. One molecule out of 184 cannot account for a 0.055 Å mean shift; noted, not removed.

## Limits
Single training seed, one σ_b, 10 epochs on CPU, not converged. The GEOM references are xTB geometries and TD's
conformers have RDKit local geometry, so the head is learning a systematic RDKit → xTB local-structure correction; the
missing control is a torsion-constrained xTB (GFN2) relaxation of the same TD conformers.

# Addendum A results: constrained GFN2-xTB control and fixed read-out (job 49899, CPU, 32 cores)
Job 49899: 2026-10-31 01:16 → 01:42 UTC (06:46 → 07:12 IST). 114 tests passed, 3 skipped (incl. the new
rigid-molecule test). xTB dry run (3 molecules) 0 errors before the full run. GPU: probe 49881 still held with
`AssocGrpGRESMinutes` (cancelled); everything here is CPU.

Absolute E1 (seed 0 / seed 1):

| arm | MAT-R | MAT-P | COV-R@0.25 | COV-P@0.25 | COV-R@0.5 | COV-P@0.5 |
|---|---|---|---|---|---|---|
| TD | 0.170 / 0.165 | 0.215 / 0.206 | 82.8 / 80.9 | 71.4 / 71.4 | 94.9 / 96.2 | 91.0 / 92.9 |
| TD + L(proj)fix | 0.111 / 0.106 | 0.145 / 0.134 | 87.8 / 88.0 | 81.5 / 83.0 | 95.7 / 96.8 | 93.8 / 94.8 |
| TD + L(noproj)fix | 0.112 / 0.107 | 0.146 / 0.136 | 87.8 / 88.4 | 81.2 / 82.9 | 95.7 / 96.8 | 93.7 / 94.7 |
| **TD + xTB (constrained)** | **0.097 / 0.089** | **0.123 / 0.113** | 88.3 / 88.9 | 83.3 / 84.5 | 95.6 / 96.6 | 93.8 / 94.8 |
| TD + MMFF94s | 0.163 / 0.152 | 0.182 / 0.170 | 82.5 / 83.0 | 78.8 / 79.7 | 94.2 / 96.3 | 92.9 / 93.8 |

Paired differences, mean [95 % CI]:

| contrast | seed | MAT-R | MAT-P | COV-R@0.25 | COV-P@0.25 |
|---|---|---|---|---|---|
| **L(proj)fix − xTB (addendum primary)** | 0 | **+0.0138 [0.0083, 0.0187]** | +0.0220 [0.0154, 0.0299] | −0.56 [−1.39, 0.15] | −1.81 [−3.33, −0.67] |
| L(proj)fix − xTB | 1 | +0.0174 [0.0122, 0.0241] | +0.0208 [0.0154, 0.0272] | −0.85 [−2.39, 0.42] | −1.47 [−2.83, −0.25] |
| xTB − TD | 0 | −0.0726 [−0.0828, −0.0634] | −0.0922 [−0.1057, −0.0794] | +5.55 [3.25, 8.14] | +11.90 [8.74, 15.39] |
| xTB − TD | 1 | −0.0767 [−0.0876, −0.0665] | −0.0931 [−0.1057, −0.0813] | +7.99 [4.61, 11.76] | +13.08 [9.70, 16.79] |
| L(proj)fix − TD | 0 | −0.0589 [−0.0670, −0.0514] | −0.0702 [−0.0805, −0.0609] | +4.99 [2.76, 7.63] | +10.10 [7.10, 13.24] |
| L(proj)fix − L(noproj)fix | 0 | −0.0009 [−0.0021, 0.0002] | −0.0018 [−0.0029, −0.0008] | −0.02 [−0.37, 0.35] | +0.22 [−0.02, 0.48] |
| L(proj)fix − L(noproj)fix | 1 | −0.0012 [−0.0023, −0.0003] | −0.0019 [−0.0029, −0.0009] | −0.39 [−1.15, 0.20] | +0.13 [−0.33, 0.53] |

Read-out diagnostics: xTB 0 errors, 5 / 5 non-converged optimisations (last geometry kept), W drift ~0, final
geometry-guard fail 103 / 109 conformers (bond or clash check relative to the TD seed; chirality 0), median heavy-atom
displacement 0.12 Å, **0.36 CPU s / conformer**. L(proj)fix: 0 errors (fix works), 1 / 0 final geometry fail,
**1.67 CPU s / conformer** (32 workers; 1.35 at 16). Same 184 paired molecules and the same evaluate_confs failure counts
as TD in every arm. Slurm again reported oom_kill events (80) with every output complete and no step failed.

## Addendum A verdicts (pre-declared rules)
- **Learned head beats xTB: no.** Constrained GFN2-xTB is better on MAT-R by 0.014 Å (seed 0) and 0.017 Å (seed 1),
  CIs exclude 0, and on MAT-P and COV-P@0.25; COV-R differences include 0.
- **Speed: no advantage either** on CPU: xTB costs 0.36 s per QM9 conformer vs 1.67 s for the 20-step read-out.
- So on GEOM-QM9 the V2 geometry head is a less accurate, slower approximation of a torsion-constrained xTB
  relaxation. The T1 "improves on TD" result stands, but it is explained by the reference geometries being xTB.
- The projector conclusion is unchanged with the fix (≈ 1 mÅ, seed-0 CI includes 0).
- Side result worth keeping: TD + constrained xTB (MAT-R 0.097 / 0.089) is a much stronger post-hoc baseline than TD
  (0.170 / 0.165) on this metric; any future geometry method on GEOM must be compared against it.

## Code-review follow-up (2026-10-31, job 49940, 02:00 UTC / 07:30 IST; diagnostic only, no reported number replaced)
An independent review of the code diff flagged that `w_drift_deg` is measured after the exact dihedral restore, so it
could not show whether the methods kept the torsions themselves (e.g. an xtb vs RDKit dihedral sign mismatch).
Re-run of seed 0 with the drift also logged on the raw output (`e1_signcheck/`):
- **xTB:** raw W drift median 0.25°, 95th percentile 0.64°; 5,327 / 5,329 conformers within 5° (1 above 30°, a
  non-converged optimisation). The constraint holds, so xtb's sign convention matches RDKit's. The rerun is
  bit-identical to the reported `td0_xtb` (max coordinate difference 0.0).
- **MMFF control:** unconstrained MMFF moves the W dihedrals a lot (raw drift median 38°, 95th percentile 114°), and
  restoring them afterwards undoes most of the relaxation and causes its guard failures. MMFF is therefore a weak
  control; the constrained xTB arm is the meaningful non-learned control.
Fixes applied after the review (none affects a reported number): `td_loss` returns 0 instead of NaN for an all-rigid
batch (new test; 115 tests pass); a failing Pool initializer becomes per-conformer errors instead of a hang; the
training script copies logs and returns the real exit code when training fails; `jmcg_tests.sbatch` now fails when a
suite fails. Known and left as is: the original L arms in `compare_seed*.json` contain the 46 unrefined rigid
conformers (disclosed above; the `fix` arms supersede them); `heavy_displacement_rms` is not alignment-invariant.

## Pointer (added 2026-10-31, nothing above changed)
Verification of these outputs, the OOM explanation, the rigid-molecule root cause and the geometry decomposition are in
`../EXP1_V2_T2/RESULTS.md`. The audit reproduced every number in this file exactly; the OOM events did not truncate any
output.
