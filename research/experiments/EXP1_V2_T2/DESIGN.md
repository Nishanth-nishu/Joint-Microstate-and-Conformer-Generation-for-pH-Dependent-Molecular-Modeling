# EXP1-V2-T2: verification, constrained-xTB analysis, longer training with three seeds (pre-declared 2026-10-31)

Separate from EXP1-V2-T1 (`../EXP1_V2_train/`), whose pre-registered results are not changed. T1 checkpoints stay in
`$P/workdir/exp1/T1_{proj,noproj}_s0/` and are only read or copied.

## Part A: verification of the saved T1 outputs (task 1), done before any new training
- `t2_audit_e1.py` (job 49955): inputs, outputs, RMSD matrices, recomputed metrics and bootstrap vs the saved JSONs.
- `exp1_t2_oom_probe.sbatch`: seed-0 re-run of the 47249 steps with the cgroup oom_kill counter logged per phase, and
  every output compared with the saved one.

## Part B: rigid molecules (task 2)
Root cause, fix, tests and the re-run are documented in `RESULTS.md`; the corrected complete evaluation is the T1
"fix" arms (job 49899). The original arms stay as the historical result.

## Part C: constrained xTB and geometry decomposition (task 3, task 5)
Constrained GFN2-xTB arm: pre-declared in `../EXP1_V2_train/DESIGN.md` addendum A, run in 49899, sign convention checked
in 49940. Added here (`t2_geom_validate.py`): per arm, connectivity re-perceived from 3D, stereo assigned from 3D vs the
SMILES and vs the TD input, non-bonded clashes, and the local error against the nearest reference conformer split into
bond lengths, bond angles, ring dihedrals and acyclic dihedrals; paired bootstrap of arm − TD and arm − xTB.
Not an equivalence claim: the head and xTB do not share constraints or objective (the head is trained against
GEOM's xTB-optimised references through the pure-L target; xTB minimises GFN2 energy with harmonic W constraints).

## Part D: longer training, three seeds (task 4)
**Arms × seeds:** {proj, noproj} × training seeds {0, 1, 2}; six runs.
**Unchanged from T1:** every flag (`--w_tau 0 --w_geo 1 --geo_norm base --l_head_init zero --sigma_b 0.15 --lr 1e-3
--batch_size 32 --limit_val_iters 100 --num_workers 4`), data split (round-2 paired cache: 82,895 train / 10,376 val;
E1 test molecules untouched), validation set and its fixed noise, checkpoint rule (`best_model.pt` = lowest validation
geometry objective over all epochs). Code: the T1 code plus the review fixes (all-rigid batch loss 0 instead of NaN,
rigid-molecule forward path); neither changes a batch that has a rotatable bond.
**Budget:** 20 epochs per run (2 × T1). Seed 0 continues from T1's `last_model.pt` (copied into its own T2 directory,
`--resume`, which restores model, optimiser, scheduler and RNG states), so seed 0 is one 20-epoch run whose first 10
epochs are T1. Seeds 1 and 2 train from scratch. The budget is fixed now and not extended after looking at E1.
**Compute:** CPU (plafnet2 GPU-minute quota exhausted: probe 49881 held with `AssocGrpGRESMinutes`). 8 CPUs per run
(OMP 4 threads + 4 loader workers), at most what the node's free CPUs allow at a time. Expected ≈ 3–5 days of wall time.
If the GPU quota is restored mid-run, a run may be resumed on a GPU; the device used per epoch is recorded.
**Recorded per run:** train/val geometry loss and base, torsion loss (w_τ = 0, so recorded only), per-module gradient
norms (train_log.jsonl), lr, epoch time, CPU count, peak memory (sacct MaxRSS), failures. After training, on validation
molecules only: the geometry-loss decomposition of `v2_diag.py` (model / base by t quintile, 1,000 samples).
**E1 (run once, after all six runs finish):** the fixed read-out of each `best_model.pt` on the frozen-TD E1 samples
(TD sampling seeds 0 and 1), TD `evaluate_confs.py`, and `t2_geom_validate.py`. Read-out CPU s / conformer recorded.

**Endpoints and rules (fixed now):**
1. *Improves on TD after longer training:* MAT-R (L(proj) − TD), TD sampling seed 0, paired bootstrap 95 % CI, per
   training seed: "yes" only if all three CIs exclude 0 favourably and TD sampling seed 1 has the same sign.
2. *Persists across seeds:* the between-training-seed SD of MAT-R and the seed-averaged effect (per-molecule metrics
   averaged over the three training seeds, then paired bootstrap over molecules).
3. *Projector:* seed-averaged MAT-R (proj − noproj): "improves" only if the CI excludes 0 in favour of proj and at least
   2 of 3 training seeds have the same sign.
4. *vs constrained xTB:* seed-averaged MAT-R (L(proj) − xTB): "beats" if the CI excludes 0 favourably; "indistinguishable"
   if it includes 0; otherwise "worse".
5. *Longer training helped:* MAT-R of T2 seed 0 (20 epochs) − T1 seed 0 (10 epochs), proj, paired CI.
6. *Where the gain is:* bond / angle / ring-dihedral / acyclic-dihedral MAE differences (Part C) and COV at 0.25 vs 0.5 Å.
7. *Cost:* CPU s per conformer, read-out vs constrained xTB, same node and worker count.

**Not testable with the data here:** generalisation beyond the training distribution. Only GEOM-QM9 data and TD's QM9
models exist on the cluster (no GEOM-DRUGS data or DRUGS model). E1 test molecules are QM9 molecules from the same
distribution; this is reported as an open question, not answered by stratifying QM9.

## Addendum B (2026-10-31 19:50 UTC / 2026-11-01 01:20 IST, before any T2 E1 number): final evaluation and decision rule
Written while the six runs train (seed 0 at epoch 14, seeds 1–2 at epochs 2–4); no T2 checkpoint has been read out.
**Procedure** (`slurm_jmcg/exp1_v2_t2_e1.sbatch`, submitted as 51354 with `afterany` on the six training jobs; if a run is
incomplete and its resume job is queued, it resubmits itself and does nothing else):
1. each run has epochs 0–19 in `epochs.jsonl` and `best_model.pt`; both checkpoints match the `checkpoints.sha256` the
   training job wrote; E1 inputs match `e1_manifest.sha256` (written by `exp1_t2_manifest.sbatch`, files unmodified since
   the part-A audit);
2. run-once guard: refuses if `$P/results/exp1/T2_e1` exists (a crash after this point is reported, not silently re-run);
3. the fixed L read-out (T1 "fix" code path, guard, torsion restoration, 32 workers, CPU) of each `best_model.pt` on TD
   sampling seeds 0 and 1, wall time / cgroup peak memory / node CPU load per read-out; `evaluate_confs.py`;
4. `t2_geom_validate.py` with every arm (TD, ETKDGv3, MMFF, constrained xTB, T1 L(proj)fix / L(noproj)fix, the six T2
   read-outs) on the molecules paired across all arms; `t2_final_analysis.py` (endpoints 1–7 mechanically, plus
   training curves, gradient statistics by epoch, read-out guard skips and cost); `v2_diag.py` on 1,000 validation
   samples; `sacct` of the training jobs (elapsed, CPU time, MaxRSS).
Predictions (`*.pkl`), RMSD matrices (`*.eval.pkl`) and logs stay in `$P/results/exp1/T2_e1`; JSON/log copies in `e1/`.
If the paired set differs from T1's 184 molecules, the difference is named and both are reported.

**Decision rule** (fixed now; the options of the 2026-10-31 brief):
- **A, continue developing the head:** only if endpoint 1 is "yes" *and* endpoint 4 is "beats" or "indistinguishable"
  *and* the read-out costs no more CPU time per conformer than constrained xTB.
- **B, keep as an optional refinement component:** endpoint 1 "yes" and endpoint 4 "worse", but the head is cheaper than
  constrained xTB per conformer or has a demonstrated use xTB cannot serve.
- **C, stop the geometry-head architecture, redirect to the joint microstate–conformer work:** endpoint 1 "yes" or "no",
  endpoint 4 "worse", and the head is not cheaper than xTB (it is then dominated on accuracy and cost), or endpoint 1 "no".
- **D, one more narrow experiment:** only if a specific hypothesis remains whose test is feasible with the data and
  compute here and could change the choice among A–C. Generalisation to GEOM-DRUGS is not feasible here (no data or
  model) and does not qualify. Longer training is not a default; it qualifies only if endpoint 5 shows a clear gain *and*
  validation loss still falls at epoch 20 *and* closing the measured gap to xTB is plausible from the slope.

## Addendum C (2026-11-02 09:30 UTC / 15:00 IST, before any T2 E1 number): exploratory analyses requested by the final brief
Added to `t2_final_analysis.py`. Nothing in addenda A–B changes: endpoints 1–7 and the decision rule stay as written.
None of the analyses below can change an endpoint verdict, and all are reported as **exploratory**.
- **Strata (graph only, RDKit):**
  - rotatable bonds (`CalcNumRotatableBonds`): 0 / 1 / ≥ 2;
  - ring count: 0 / 1 / ≥ 2;
  - presence of a 3- or 4-membered ring.

  Per stratum: n, and seed-averaged MAT-R for proj − TD, proj − xTB and proj − noproj, with bootstrap CIs.
  No multiplicity correction; strata with n < 20 are not interpreted.
- **Near-duplicate sensitivity:** endpoints 1, 3 and 4 recomputed without the one E1 molecule whose
  stereo-free graph occurs in the training cache (`N#CC12C3C4CN1C2C43`, `t1_leak_check.py`).
- **Variability, reported separately:**
  - training-seed SD of the mean MAT-R (3 seeds);
  - TD sampling seed (seed 0 vs seed 1);
  - molecule-level bootstrap CI width.

  These are reported side by side and not combined.
- **Denominators:** 200 attempted, molecules with TD output per sampling seed (185), paired molecules
  (expected 184).
- **Bootstrap draws:** the addendum-C bootstraps use their own generator (`default_rng(1)`), so the endpoint CIs
  are the same as without addendum C.
- **Training log:** if an (epoch, iteration) record occurs twice (noproj s2 after its cancel and resume), the
  later record is kept.

## Addendum D (2026-11-04 05:05 UTC / 10:35 IST, after training, before any T2 E1 number): evaluation resources
The six runs finished at 2026-11-04 03:39 UTC. The queued evaluation (51354, resubmitted as 56777: 32 CPUs, 90 G)
had an estimated start of 2026-11-08 because other users share gnode118 (10 CPUs and about 25 GB free). At the user's
request it was cancelled unstarted and resubmitted with 8 CPUs and 22 G, leaving 2 CPUs free for others.
- **Workers:** 6 read-out / evaluate_confs / geometry workers instead of 32 / 16 / 16. The read-out is deterministic
  per conformer, and `evaluate_confs` and the geometry checks are deterministic per molecule, so the predictions,
  RMSD matrices and every accuracy metric are unaffected. Only wall time changes.
- **Cost (endpoint 7):** time per conformer depends on the worker count and node load, so the T2 read-out timings are
  not comparable with the T1 xTB timing at 32 workers. Step 4b re-times the head (T2 proj s0) and constrained xTB on
  the same first 20 molecules of TD seed 0, with the same 6 workers, in the same job. Endpoint 7 uses that pair;
  the full-run timings are reported as secondary.
- Nothing else changes: endpoints, rules, checkpoints and molecules are as in addenda B and C.
