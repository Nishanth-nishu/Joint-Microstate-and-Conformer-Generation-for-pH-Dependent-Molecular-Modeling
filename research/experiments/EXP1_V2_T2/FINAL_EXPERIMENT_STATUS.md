# EXP1-V2-T2: experiment status (audit 2026-11-02 09:08 UTC / 14:38 IST; final update 2026-11-04 08:40 UTC / 14:10 IST)

Read from `squeue`, `sacct`, the mirrored training logs in `train/T2_*/` and the job logs in `jmcg_stage/logs/`.
Nothing was re-run for this audit. Updated again when the final evaluation finishes.

## 1. Training jobs (6 configurations: projector on/off × training seeds 0, 1, 2; 20 epochs each)

| run | job(s) | state | epochs logged (of 20) | notes |
|---|---|---|---|---|
| proj s0 | 49962 | COMPLETED, exit 0, 2026-11-01 19:04 UTC (02 Nov 00:34 IST), elapsed 1-16:31 | 0–19 | resumed T1 at epoch 10 ("RESUMED at epoch 10, best val 0.2294") |
| noproj s0 | 49964 | COMPLETED, exit 0, 2026-11-01 17:34 UTC (23:04 IST), elapsed 1-15:01 | 0–19 | resumed T1 at epoch 10 (best val 0.2778) |
| proj s1 | 49966 | RUNNING (2-06:31) | 0–13 | from scratch |
| noproj s1 | 49968 | RUNNING (2-06:31) | 0–13 | from scratch |
| proj s2 | 49970 | RUNNING (1-23:37) | 0–11 | from scratch |
| noproj s2 | 49972 → 52355 | 49972 CANCELLED by user request 2026-11-01 17:23 UTC (22:53 IST), right after epoch 7 (8th) was saved, to free 2 CPUs; 52355 RUNNING (15:44) on 6 CPUs | 0–10 | 52355: RESUME=1, "RESUMED at epoch 8 (best val 0.2870)", 116 tests passed |

All six configurations are present. Settings are identical to T1 for every run (same flags; `train_jobs.txt`).
The only deviation is noproj s2's CPU allocation (6 instead of 8), which changes wall time, not the computation
(same thread and worker counts).

**Seed-0 history:** epochs 0–9 of both T2 seed-0 runs are byte-identical to T1's `epochs.jsonl`. The copied
T1 checkpoints had the same sha256 in T1 and T2 at the time of copying (job log 49962:
T1/T2 `last_model.pt` cc877b58…, `best_model.pt` 09a0d876…). Final T2 hashes are in `train/T2_*_s0/checkpoints.sha256`.

**Consistency:** no duplicated (epoch, iteration) record in any `train_log.jsonl`, except one in noproj s2
(epoch 8, iteration 0), logged by 49972 before it was cancelled and again by 52355. The analysis keeps the
later record. `epochs.jsonl` has no duplicate epoch.

**Validation (selection rule: lowest validation objective):** both seed-0 runs pick epoch 19, and the
objective was still falling at epoch 20: proj 0.2294 (epoch 10) → 0.2006, noproj 0.2778 → 0.2483. Seeds 1 and 2
are also still improving at their latest epochs.

## 2. Resume mechanism
- Each training job has an `afternotok` resume job with RESUME=1, which refuses to start without
  `last_model.pt`: 49976 (proj s1), 49977 (noproj s1), 49978 (proj s2), 52356 (noproj s2, 6 CPUs).
- 49974 and 49975, the seed-0 resume jobs, became DependencyNeverSatisfied once seed 0 completed and were
  cancelled at 09:10 UTC.
- 49979, the original noproj-s2 resume job, was cancelled 2026-11-01 before 49972 was stopped, so that it
  could not start a second copy; 52356 replaces it.
- A watcher cancels each remaining resume job when its training job completes, so the evaluation's
  dependency chain cannot hang.
- The resume path was exercised twice: in T1→T2 seed 0, and in 52355.

## 3. Final evaluation readiness
- **51354** (`exp1_v2_t2_e1.sbatch`, 32 CPUs, 90 G, 16 h), pending on `afterany` of the six original training
  jobs.
  - It checks all six runs: epochs 0–19, checkpoints present and matching the hashes the training job wrote.
  - If a run is incomplete while a T2 job is still queued, it resubmits itself after that job.
  - It refuses to run twice.
- **Manifest (51353, done):** `e1_manifest.sha256` covers the 24 E1 inputs. All were last modified
  2026-10-28/29, before the part-A audit.
- **Dry run (51353):** with the T1 outputs standing in for the six read-outs, every expected value was
  reproduced:
  - stand-ins equal T1 (MAT-R 0.111 / 0.106 Å);
  - "20 vs 10 epochs" exactly 0;
  - between-seed SD 0;
  - geometry as in part C.

## 4. Expected completion (at the rates so far; UTC / IST)
| run | rate | expected end |
|---|---|---|
| proj s1, noproj s1 | ≈ 3.9 h / epoch | ≈ 2026-11-03 09:00 / 14:30 |
| proj s2 | ≈ 4.0 h / epoch | ≈ 2026-11-03 17:00 / 22:30 |
| noproj s2 (6 CPUs) | ≈ 5.2 h / epoch | ≈ 2026-11-04 08:00 / 13:30 |

Each run stays within its 4-day limit. The evaluation then needs about 2–4 h.

## 5. Missing or duplicated outputs
None found. No T2 E1 output exists yet (`$P/results/exp1/T2_e1` is created only by 51354).

## 6. Final status (2026-11-04 08:40 UTC / 14:10 IST)
Sections 1–5 are the 2026-11-02 audit and are kept as written.

**Training.** All six runs COMPLETED, exit 0, epochs 0–19, best epoch 19 in every run, no non-finite loss.
| run | job(s) | end (UTC / IST) |
|---|---|---|
| proj s0 | 49962 | 2026-11-01 19:04 / 02 Nov 00:34 |
| noproj s0 | 49964 | 2026-11-01 17:34 / 23:04 |
| proj s1 | 49966 | 2026-11-03 06:08 / 11:38 |
| noproj s1 | 49968 | 2026-11-03 05:38 / 11:08 |
| proj s2 | 49970 | 2026-11-03 16:23 / 21:53 |
| noproj s2 | 49972 → 52355 | 2026-11-04 03:36 / 09:06 |

Every `afternotok` resume job (49974–49979, 52356) was cancelled unstarted once its training job completed. The
resume path ran only by design (seed 0 from T1; noproj s2 after the CPU release).

**Final evaluation.**
- **Submissions:** 51354 found noproj s2 incomplete and resubmitted itself as 56777 (32 CPUs, 90 G, estimated start
  2026-11-08). At the user's request 56777 was cancelled unstarted, and the evaluation was resubmitted as **57155**:
  8 CPUs, 22 G, 6 workers, leaving 2 CPUs free (DESIGN.md addendum D, written before any result).
- **57155:** started 05:05 UTC / 10:35 IST, ended 08:27 UTC / 13:57 IST, COMPLETED, exit 0, elapsed 3 h 29 min.
- **Checks:** 12 / 12 checkpoint hashes OK; 24 / 24 E1 inputs OK; 116 tests passed (3 skipped); one run only
  (`STARTED`, `DONE`).
- **Read-outs:** 12 of 12, 0 errors.
- **Pairing:** 184 paired molecules on both TD sampling seeds, the same set as T1.

**Outputs present** (`e1/`):
- `final_analysis.json` / `.log`;
- `geom_validate_seed{0,1}.json` / `.log`;
- 12 × `td{s}_T2_*.{readout.log,eval.log,pkl.log.json,time}`;
- `timing_{Lproj,xtb}.*`;
- `v2_diag_val.json` / `.log`;
- `train_sacct.txt`, `checkpoints.sha256`, `code.md5`.

Predictions and RMSD matrices are in `/scratch/nishanth.r/jmcg/results/exp1/T2_e1/`. Nothing is missing.
Nothing was deleted: every checkpoint, log and earlier output remains.

Results: `FINAL_RESULTS.md`. Interpretation: `FINAL_SCIENTIFIC_ASSESSMENT.md`. Decision: `RESEARCH_DECISION.md`.
