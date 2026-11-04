# EXP1-V2-T2 results (design: `DESIGN.md`; T1's pre-registered results in `../EXP1_V2_train/RESULTS.md` unchanged)

Times UTC / IST. All CPU (plafnet2 GPU-minute quota exhausted). Parts A–C are complete; part D (training) is running.

## A. Verification of the saved T1 E1 outputs (task 1)
**Audit** (`t2_audit_e1.py`, job 49955, 2026-10-31 02:22 / 07:52; `audit_e1.json`):
- Inputs: subset 200 molecules. TD produced conformers for 185 per sampling seed (15 TD failures each, partly different
  molecules between seeds), 5,390 / 5,386 conformers, every molecule exactly 2 × n_true.
- Arm outputs (L(proj), L(noproj), MMFF, xTB, both fix arms, both seeds): no missing or extra molecule, no count
  mismatch, same atoms in the same order, no NaN coordinate, no duplicated conformer.
- RMSD matrices: one per scored molecule, columns = that molecule's conformer count in every arm.
- Paired set: 184 per seed = 185 − 1 molecule (`[H]N=c1nc[nH]cn1`) whose alignment fails in every arm (NaN RMSD).
  Every arm pairs on exactly those 184; no arm lost a molecule the others kept.
- Metrics recomputed from the RMSD matrices: every paired mean in all 10 compare JSONs reproduced exactly (max
  difference 0.0). Bootstrap re-implemented: identical CIs with the same generator stream; with a fresh generator the
  CI ends move by ≤ 0.33 (coverage-percent scale) / Monte-Carlo error, as expected for 5,000 draws.
- **Which observations enter each metric:** per molecule, MAT-R = mean over its reference conformers of the minimum
  RMSD to its generated conformers, MAT-P the transpose, COV-R/P@t the share of those minima below t (strict <); then
  the unweighted mean over the 184 paired molecules; CIs resample molecules (5,000 draws, percentile).
- **Why TD's own printed coverage differs** (e.g. COV-R@0.25 80.74 vs 87.76 in our table): `evaluate_confs.py` counts the
  15 molecules with no TD output as 0 % coverage and averages over 200; our table averages over the 184 paired
  molecules. MAT agrees because TD uses `nanmean` over scored molecules. No error, two denominators.

**OOM status** (job 49957, 02:27–02:36 / 07:57–08:06; `oom_probe/`): the seed-0 steps of 47249 re-run with the same
resources and the cgroup oom_kill counter logged every 2 s. All 5 kills happened in the first 15 s of the learned
read-out step, while 16 freshly forked workers each loaded torch, TD and the model (2–3 GB each; cgroup peak 37.9 of
40 GB). The Pool replaced the killed workers (process count constant at 22) and returned every result. No kill during
MMFF, evaluate_confs or the comparison. Every output is bit-identical to the saved one (max coordinate difference 0.0,
max RMSD difference 0.0, same molecules), and MAT-R(L − TD) reproduces (−0.0589). 47249's 18 and 49899's 80 events are
the same start-up mechanism (49899 used 32 workers). Conclusion: no output or metric was truncated.

## B. Molecules without rotatable bonds (task 2)
- **Root cause** (TD `diffusion/score_model.py`, `TensorProductConvLayer.forward`): `out_nodes = out_nodes or
  node_attr.shape[0]`. With no rotatable bond `edge_mask.sum()` is 0, which is falsy, so the bond convolution silently
  returns one row per atom; the σ-scaling then multiplies 18 rows by 0 σ values and raises.
- **Fix** (`FlexiTorsScoreModel.forward`): for a batch with no rotatable bond, one bond is marked rotatable only to run
  TD's forward, then the mask is restored and the torsion output set to the empty tensor. TD's atom layers build their
  graph from `edge_index` and never read `edge_mask`, so Δ̂ is exactly the model's; the torsion output is correctly empty.
  Batches with any rotatable bond take TD's unchanged path. `td_loss` returns 0 (not NaN) for an all-rigid batch.
- **Tests:** `test_rigid_molecule_geometry_head` (Δ̂ equals the one with any bond marked, empty torsion output, mask
  restored), `test_rigid_molecule_in_mixed_batch` (a rigid and a flexible molecule batched give the same torsion and
  geometry outputs as alone), `test_td_loss_all_rigid_batch_is_zero_not_nan`. 116 passed, 3 skipped.
- **Re-run:** the "fix" arms of job 49899 (0 read-out errors; all 12 rigid molecules refined). The complete evaluation is
  L(proj)fix − TD MAT-R −0.0589 [−0.0670, −0.0514] (seed 0), −0.0593 [−0.0672, −0.0522] (seed 1). The original arms
  (46 unrefined conformers) stay as the historical T1 result.

## C. Constrained xTB and where the gain is (tasks 3 and 5)
Arm: GFN2-xTB 6.7.1, every W dihedral harmonically constrained (1.0 Eh) to TD's value, then restored exactly
(pre-declared in T1 addendum A; job 49899). Constraint check (49940): raw W drift median 0.25°, 5,327 / 5,329 within 5°.
Not an equivalence claim: the head and xTB share neither objective nor constraints.

E1 (184 paired molecules; seed 0 / seed 1):

| arm | MAT-R | MAT-P | COV-R@0.25 | COV-P@0.25 | CPU s / conformer |
|---|---|---|---|---|---|
| TD | 0.170 / 0.165 | 0.215 / 0.206 | 82.8 / 80.9 | 71.4 / 71.4 | – |
| TD + L(proj)fix | 0.111 / 0.106 | 0.145 / 0.134 | 87.8 / 88.0 | 81.5 / 83.0 | 1.67 |
| TD + L(noproj)fix | 0.112 / 0.107 | 0.146 / 0.136 | 87.8 / 88.4 | 81.2 / 82.9 | 1.65 |
| TD + constrained xTB | **0.097 / 0.089** | **0.123 / 0.113** | 88.3 / 88.9 | 83.3 / 84.5 | **0.36** |
| TD + MMFF94s | 0.163 / 0.152 | 0.182 / 0.170 | 82.5 / 83.0 | 78.8 / 79.7 | 0.02 |
| ETKDGv3 | 0.229 / 0.217 | 0.227 / 0.221 | 69.5 / 68.2 | 67.9 / 67.6 | – |

L(proj)fix − xTB MAT-R +0.0138 [0.0083, 0.0187] / +0.0174 [0.0122, 0.0241]: xTB is better on both seeds.

**Geometry validation and error decomposition** (`t2_geom_validate.py`, job 49960; `geom_validate_seed{0,1}.json`).
Errors vs the reference conformer nearest to each generated conformer, heavy atoms, mean over 184 molecules, seed 0
(seed 1 within ±0.05 Å / ±4° of these, same ordering):

| arm | bond MAE (Å) | angle MAE (°) | ring-dihedral MAE (°) | acyclic-dihedral MAE (°) | graph re-perceived wrong | stereo changed vs TD input | conformers with a clash < 2 Å |
|---|---|---|---|---|---|---|---|
| TD | 0.0319 | 3.30 | 14.3 | 40.1 | 136 | – | 22 |
| L(proj)fix | 0.0128 | 1.25 | 9.1 | 33.6 | 29 | 2 | 5 |
| L(noproj)fix | 0.0121 | 1.09 | 9.0 | 31.4 | 23 | 2 | 5 |
| constrained xTB | **0.0027** | **0.68** | **7.1** | **25.8** | 3 | 4 | 1 |
| MMFF94s | 0.0194 | 2.19 | 10.9 | 40.2 | 6 | 79 | 8 |
| ETKDGv3 | 0.0319 | 3.25 | 14.3 | 41.4 | 108 | – | 1 |

(5,388 conformers per arm. "Graph re-perceived wrong": RDKit `DetermineConnectivity` from the 3D coordinates does not
reproduce the bond graph, i.e. a bond stretched or two atoms too close. Stereo assigned from 3D vs the SMILES fails for
~270 conformers in every TD-derived arm and 191 for ETKDG, mostly centres the SMILES leaves unspecified; the comparison
with the TD input conformer is the meaningful one for refinement. The read-out's own guard flags ~105 xTB conformers per
seed because it measures bonds relative to TD's distorted seed geometry; the absolute checks above show xTB's
geometries are the most valid.)

Paired differences (seed 0, 95 % CI): L(proj)fix − xTB bond MAE +0.0101 [0.0094, 0.0107] Å, angle +0.57° [0.51, 0.64],
ring dihedral +2.0° [1.5, 2.7], acyclic +7.8° [4.2, 12.0]; seed 1 the same signs and sizes.

**Reading.** xTB reproduces GEOM's bond lengths almost exactly (0.0027 Å MAE), as expected if the references are
GFN2-xTB geometries. The trained head removes ≈ 60 % of TD's bond-length error and ≈ 62 % of its angle error, xTB ≈ 92 %
and ≈ 79 %. The head's gain is local structure (bonds, angles) and ring puckering, plus the acyclic dihedrals that are not
TD's W torsions; coverage moves at 0.25 Å (+5–7 points) but barely at 0.5 Å. Constrained xTB does better on every
component, is ≈ 4.7× cheaper per conformer on CPU, and yields more valid geometries. The evidence supports the confound:
on GEOM-QM9 the head mostly learns an approximation of an xTB-style local relaxation.

## D. Longer training, three seeds (task 4): running
Jobs (`train_jobs.txt`): proj / noproj × seeds 0, 1, 2 = 49962, 49964, 49966, 49968, 49970, 49972, each with an
`afternotok` resume job (49974–49979). Seed 0 resumed T1 at epoch 10 (verified: "RESUMED at epoch 10", T1 best values
0.2294 / 0.2778); seeds 1 and 2 from scratch. 116 tests passed in each job before training. Seeds 0 and 1 started
2026-10-31 02:32–02:37 UTC (08:02–08:07 IST); seed 2 waits for CPUs (node at 48/48; the per-user QOS caps CPUs and
4-day wall time). Logs mirror every 10 min into `train/T2_*_s*/`. E1 evaluation of the six checkpoints runs once after
all six finish, as pre-declared; no T2 number exists yet.

**Status 2026-10-31 19:35 UTC / 2026-11-01 01:05 IST:** all six running, no failure, no resume triggered.
Seed 2 started 09:31 UTC (15:01 IST) once CPUs freed. Epochs completed (of 20): s0 14 / 14 (proj / noproj),
s1 4 / 4, s2 2 / 2. ≈ 3.9–4.0 h of training per epoch on 8 CPUs, plus validation. Projected finish:
- s0: ≈ 2026-11-01 20:00 UTC;
- s1: ≈ 2026-11-03 16:00 UTC;
- s2: needs ≈ 100 h, more than the 4-day limit, so its `afternotok` resume job is expected to finish the last epoch
  or two (≈ 2026-11-04/05).

Final evaluation and decision rule: DESIGN.md addendum B (job 51354, waiting).

**Final (2026-11-04).** All six runs completed 20 epochs; the final evaluation ran once (57155, 05:05–08:27 UTC /
10:35–13:57 IST). Full tables: `FINAL_RESULTS.md`. Seed-averaged MAT-R (TD seed 0), 184 paired molecules:
- proj − TD −0.0634 [−0.0724, −0.0552], all three training seeds and TD seed 1 agree → endpoint 1 **yes**;
- proj − noproj −0.0002 [−0.0013, +0.0011] → endpoint 3 **not shown**;
- proj − xTB +0.0093 [+0.0036, +0.0141] (seed 1 +0.0142) → endpoint 4 **worse**;
- T2 − T1 seed 0 −0.0039 [−0.0061, −0.0021] → endpoint 5 **helped** (small);
- cost, like for like on 20 molecules at 6 workers: 1.23 vs 0.30 CPU s / conformer (head 4.1× xTB).
Decision rule (addendum B) → **C**; see `FINAL_SCIENTIFIC_ASSESSMENT.md` and `RESEARCH_DECISION.md`.

## Answers, final (task 5; supersedes the table below)
| question | answer | evidence (FINAL_RESULTS.md) |
|---|---|---|
| Outperforms frozen TD after longer training? | yes | −0.063 Å MAT-R, 3 / 3 training seeds, both TD sampling seeds |
| Persists across training seeds? | yes | training-seed SD 0.5 mÅ vs effect 63 mÅ |
| Projector improves independently? | not shown, not meaningful | −0.2 mÅ [−1.3, +1.1]; ≤ 1 mÅ on TD seed 1; local geometry slightly worse with proj |
| Constrained xTB reproduces most of the gain? | yes, and exceeds it | xTB − TD −0.073 / −0.077 Å vs head −0.063; every local component better |
| Gain in local, ring or coverage? | local first, ring dihedrals second; coverage at 0.25 Å only | §5 |
| Generalises beyond the training distribution? | not testable here | no GEOM-DRUGS data or model |
| Cost justified vs physical optimisation? | no, on QM9 CPU | 1.23 vs 0.30 CPU s / conformer, like for like |

## Answers so far (task 5)
| question | answer now | evidence |
|---|---|---|
| Outperforms frozen TD after longer training? | pending (part D) | T1 (10 epochs): yes, −0.059 Å MAT-R, both TD sampling seeds |
| Persists across training seeds? | pending (part D) | one training seed so far |
| Projector improves independently? | not shown | ≈ 1 mÅ MAT-R, seed-0 CI includes 0; noproj has slightly lower bond / angle error |
| Constrained xTB reproduces most of the gain? | yes, and exceeds it | xTB − TD −0.073 / −0.077 Å vs L − TD −0.059; every local component better |
| Gain in local, ring or coverage? | local (bonds, angles) first, ring dihedrals second; coverage only at tight thresholds | part C |
| Generalises beyond the training distribution? | not testable here | no GEOM-DRUGS data or DRUGS model on the cluster; E1 is in-distribution QM9 |
| Cost justified vs physical optimisation? | no, on QM9 CPU | 1.67 vs 0.36 CPU s / conformer, and less accurate; GPU read-out cost unmeasured (no GPU) |

## Limitations and open questions
Single TD model (qm9_default) and E1 subset of 200; constraint force constant for xTB fixed at 1.0 Eh (not varied);
nearest-reference matching uses the precision-side RMSD minimum; GPU timing of the read-out unmeasured; the head might
behave differently on larger molecules where xTB is costlier, which cannot be tested with the data here.
