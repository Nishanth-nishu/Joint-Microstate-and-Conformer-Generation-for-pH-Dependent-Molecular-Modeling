# Code review P1: Stage A milestone 2 (structure, reverse loop, job machinery, budget)

Reviewer: coding planner P1, 2026-10-17. Scope: commits `0f70dc5..53029b6` and `round5/IMPLEMENTATION.md`
(milestone 2). Checked against `round5/DECISION.md`, EXP-3 §9, EXP-5 §12 and EVAL §1/§3.1. No code was edited and
nothing was committed.

**Re-run locally (Windows, Python 3.11, no torch_geometric):**
- `pytest flexitors/tests`: 89 passed, 13 skipped (cluster-only);
- `pytest tools/tests/{test_lsa_floor,test_leakfree,test_leakfree_budget}.py`: 15 passed.

Both match IMPLEMENTATION §2. The cluster goldens (`flexitors/tests/test_reverse.py`) are recorded in
`cluster_sync/round5/tests_m2/`.

**Milestone-1 review items.** All my milestone-1 items are addressed (commits `a5b5d32`, `76785ae`, `03a5d31`,
`72980b4`, `53290e9`, `02e39f6`, `b6e4aa8`, `b2e9ed4`):
- `td.py` now runs every TD import with cwd = TD root;
- `TEST_ONLY_API` is empty (`test_no_dead_code.py:19`).

**Citation spot-check.** Every E-ID cited in the new code exists in the ledgers with the claimed content:
E-CHEM-043/046/048/064, E-JOINT-026, E-SAMPLE-033/053/054/055/057/058. The README PG values in
`flexitors/sampling/generate.py:45-64` match `torsional-diffusion/README.md:78-91` digit for digit.

**Overall.** The code is clean and readable:
- one reverse loop, golden-tested bitwise against `sample()` on CPU;
- PG through TD's own `sample()` with no copy;
- per-conformer streams that make prefixes and subsets exact;
- operators that only produce torsion-update vectors applied by TD's `modify_conformer`;
- pre-declared analyses committed before the runs.

The problems are in how the job machinery charges and protects the budget, in two DRUGS endpoint details, and in the
loop's hook interface, which does not yet admit the milestone-3 geometry step.

## BLOCKER

**B1. DRUGS recall is scored at 0.5 Å.**
- `experiments/exp3/analysis.py:46-47` hard-codes `COV-R@0.5`.
- `experiments/common/infer.sbatch:149-150` runs every `eval` with `--dataset qm9`, so the evaluator's summary and the
  `report_threshold` also default to 0.5 Å (`evaluate_confs.py:163`).
- EVAL §1 and EXP-3 §4 (ii) set δ = 0.75 Å for GEOM-DRUGS.
- Failure scenario: wave w3 runs, and the pre-declared co-primary DRUGS non-inferiority test (COV-R margin −0.5 pt,
  EXP-3 §9.3) is computed on the wrong coverage threshold. Because commit time is pre-registration, fixing it after
  the run is a post hoc change.
- Fix before w3:
  - make the threshold a per-dataset constant (`{'qm9': 0.5, 'drugs': 0.75}`) in `analysis.py`;
  - pass `--dataset drugs` for `drugstest` rows (`dataset_files` already knows the dataset);
  - add a synthetic test where the two thresholds give different COV values.

**B2. The 24.95 GPU-h DRY_RUN counts generation only, while the ledger charges whole GPU allocations, so the 25 GPU-h
stop will trip with work still unrun.**
- `status.py:30-43, 64` takes `used` from `sacct` (elapsed × allocated GPUs). That is the right source, but every
  post step runs inside the GPU allocation (`infer.sbatch:147-161`): `evaluate_confs.py` (`--n_workers 8`), MMFF
  `features`, the validity `screen`, `stereo_eval`.
- `GPUH_EST` (`make_tasks.py:1-6`) is conformers × 0.034 GPU-s, i.e. generation only. The smoke cannot reveal the gap,
  because its 20 molecules have small K. For td_s0 the ledger row is 101 s against 80 s of generation, so evaluation was
  about 20 % there (`cluster_sync/round5/smoke_m2/smoke_m2_summary.txt`, ROW 1).
- Evaluation cost is Σ K_i × n_i, quadratic in K. For the K = 727 molecule alone the seed-0 pool needs 727 × 2,908 ≈
  2.1 M `GetBestRMS` calls. Round 2 measured ≈ 0.26 h unpacked per plain 2K QM9 evaluation (plan_P3 §1 table,
  `round3/grounding.md:56-63`).
- Inside GPU jobs, Stage A runs about 14 QM9 2K/4K/pool evaluations plus six DRUGS evaluations. DRUGS has K ≈ 100 and
  ≈ 25 heavy atoms; plan_P3 §12 estimated 50–200 CPU-h for DRUGS RMSD work.
- My estimate is + 3–8 GPU-h of allocation that the plan does not show, mostly in w3.
- Failure scenario: w1 and w2 overrun their estimates. The per-wave check (B3) then refuses w3 (DRUGS, the dataset that
  carries Exp 3's primary endpoint, EXP-3 §9.3) or later refuses Exp 1's M1, the main experiment.
- Fix (as plan_P3 §9 `fx_dreval`): generation rows in GPU jobs; `eval`/`features`/`validity`/`stereo` in CPU jobs
  submitted `afterok` on the generation job (no `--gres`). Re-run the DRY_RUN and state the estimate as GPU-h of
  allocation.

**B3. The budget stop checks one wave at a time and reserves nothing for Exp 1** (`submit_stageA.sh:123-128`).
- `--planned "$total"` is that wave's estimate only. w1 + w2 + w3 may therefore spend up to 25 GPU-h, and Exp 1
  (milestone 3, the reason for Stage A) is the experiment that gets refused.
- Pass `planned = this wave + every later wave + EXP1_GPU_H` (`EXP1_GPU_H` is already defined at `:99`), so the stop
  fires *before* the wave that would make the total exceed 25, not after.
- Also print the M1 ceiling case: with EXP-1 §5.2's 16 GPU-h ceiling, the plan is 13.76 + 16.5 ≈ 30.3 GPU-h.

## SHOULD-FIX

**S1. The reverse loop's hook interface is Exp 3-specific and does not admit the milestone-3 geometry step**
(`flexitors/sampling/reverse.py:25-65`).
- The `hooks` object is duck-typed to `mask_update / fold_updates / verify`. It receives numpy *copies* of the
  positions and may only return torsion corrections, which are applied through `apply_torsion_and_update_pos`.
- EXP-1 §10.2 needs:
  - (i) the geometry clock in `data_gpu.node_lambda` before each model call;
  - (ii) a fixed σ input (σ_min) in τ-frozen mode;
  - (iii) no torsion step;
  - (iv) access to `data_gpu.l_pred` after the call;
  - (v) arbitrary position updates written into both `data.pos` and `conf_dataset.data[i].pos`;
  - (vi) exact torsion restoration after each step.
- None of these has a hook point. Milestone 3 would have to edit the loop, which risks the golden equality and could
  fork a second loop (DECISION graft 1: one loop).
- Do it now, while only Exp 3 uses the loop:
  - define a `typing.Protocol` with explicit points: `before_model(data_gpu, k)`, `torsion_update(perturb, k)`,
    `after_step(state, k)`;
  - `state` should expose `conf_dataset`, `data`, `data_gpu` and a `set_positions(idx, pos)` helper;
  - add `sigma_input` and `torsion_step` parameters;
  - make `StereoHooks` implement the protocol;
  - keep `test_reverse_without_hooks_equals_td_sample` and add a no-op-protocol golden.

**S2. The derived all-molecule difference uses an arbitrary N** (`experiments/exp3/analysis.py:106, :113`).
- `n_all = len(next(iter(rec.values())))` takes the first frame in `sorted(glob('*_s*'))`. For QM9 that is
  `F_filter_s0`, whose `eval.pkl` covers only the 102 affected molecules, because `infer.sbatch:149` evaluates the
  arm's own subset `confs.pkl`, not `confs_spliced.pkl`.
- So `(n_affected / N) × diff` divides 102 by 102 and reports the undiluted subset difference as the all-molecule
  value. That is the opposite of EXP-3 §9.1's dilution formula.
- For DRUGS no full O exists (`stageA_tasks.tsv` w3: O runs on `keys_drugs_double_or_triple`), so N must come from the
  csv: 1000, or the released-pickle universe per §9.1 condition 4.
- Take N from the csv of the dataset. Also fix the module docstring (`:4`, "holds the spliced eval.pkl"), which does
  not match what `infer.sbatch` scores.

**S3. The DRUGS timing row is missing.** `make_tasks.py:6` says "measured first by the timing row", and EXP-3 §9.4
requires "a measured 20-molecule timing run first" under the 4.5–5 GPU-h cap. `stageA_tasks.tsv` w3 has no such row,
and w3 is one job of 4.1 GPU-h estimated. Add a 20-molecule `exp3/drugs/timing_s0` row (generation + evaluation) as its
own job, and make `submit_stageA.sh w3` re-derive the arm estimates from its measured rate before submitting.

**S4. Re-running a job regenerates rows that already finished** (`infer.sbatch:114-165`).
- `run_row` has no skip-if-done. `generate.py` resumes only an unfinished `confs.partial.pkl` and deletes it at the end
  (`generate.py:291-292`).
- Failure scenario: a packed job (e.g. w1 seed 0: 0.52 + 4 × 0.47 + 1.43 ≈ 3.8 GPU-h of generation plus evaluations)
  hits `--time=06:00:00` (`infer.sbatch:8`). Resubmitting it reruns every finished row from scratch and pays their
  GPU-h twice.
- Write a `ROW_DONE` marker (the outputs of `fx_strict_check` plus the post-step files) and skip rows that have it.
  Size `--time` from the summed estimate with margin, and add `--signal=B:USR1@600` plus requeue as the training array
  does.

**S5. Pre-flight does not check DRUGS inputs** (`infer.sbatch:102`). It always checks `$MODEL_QM9`, never the
dataset's own model and csv. IMPLEMENTATION §5 item 5 says "the preflight refuses the rows until those files exist",
but `dataset_files` (`:104-112`) only sets variables, so a w3 row fails later inside `generate.py`. Call `fx_preflight`
per row with `$MODEL/best_model.pt "$CSV" "$MOLS"`.

**S6. The identity condition of EXP-3 §9.1 is reported, not enforced** (`splice.py:118-122`; IMPLEMENTATION §5 item 1).
- On GPU the arm-vs-O deviation is 0.265 Å, and the job still succeeds. That is correct for the reason stated (CUDA
  reductions), but then condition 3 holds only on CPU.
- I support option (b): run the QM9 Exp 3 arms and their O (the `O_prefix` source) on CPU. That is ≈ 20 min per arm and
  seed on 8 cores and 0 GPU-h, and it keeps exact pairing.
- Make `--identity_check` fail (exit 3) above 1e-6 Å when the device is CPU.
- R-CHEM must confirm (b) before w2. If the Exp 3 O is the 2K prefix of a GPU pool, the CPU arms are again unpaired
  with it, so O for Exp 3 must then be its own CPU run on the affected keys.

**S7. The OOM retry sleeps inside the allocation** (`fx_common.sh:52`: 300 · 2^(n−1) s, up to 35 min over three
retries). That time is charged to the 25 GPU-h cap via `sacct`. Either cap the total wait (e.g. ≤ 10 min) or exit with
a requeue request so the wait is not billed. Also print the waited seconds in `STATUS.md`.

**S8. The cut order needs concrete, measured levers.** DECISION item 10 says: pack M1, then Exp 5 P6, then reduced PG
tuning, then stop and ask. In this implementation P6 has **no GPU rows**: it is the selector applied to the `pg4k`
pools, which P2b (primary) needs anyway. Cutting P6 therefore saves only CPU. Recommended pre-declared order, written
into `submit_stageA.sh plan` output with the GPU-h each step saves:
1. Move post-processing to CPU jobs (B2): saves the evaluation share, about 3–8 GPU-h. This is not a scientific cut.
2. Run Exp 3 QM9 arms + O on CPU (S6): about 0.25 GPU-h, and restores exact pairing.
3. Pack the three M1 runs if the 1-epoch smoke's packing test passes (plan_P3 §4.3): saves 0 to ≈ 5 GPU-h,
   unmeasured.
4. Reduce PG validation tuning from 4 to 2 settings (one invariant, one non-invariant, by README recall): −0.94 GPU-h.
5. Drop DRUGS `K1_canonical` (−0.75), the first cut in EXP-3 §9.4's own order. Keep DRUGS K1 / K1_wrong
   (DECISION 10, EXP-3 §9.4 caveat).
6. Stop and ask the user (DECISION 10).

Record that "drop P6" saves 0 GPU-h, so nobody relies on it.

## NIT
- **N1.** `generate.py:146-164`. `run_conformers` mixes the PG branch (global-RNG seeding, TD `sample`) and the loop
  branch in one function. A `_run_pg` helper would let milestone 3 add the geometry arm without growing this function
  further.
- **N2.** `generate.py:43-44`, `make_tasks.py:15`. The PG setting names live in two places (`PG_SETTINGS` dict keys and
  the `make_tasks` list). Import `sorted(PG_SETTINGS)` in `make_tasks.py` so a renamed setting cannot silently desync
  the task table.
- **N3.** `select.py:63-65`. The unknown-method check runs after the method has been used to pick `keep` and `start`.
  Move `if method not in SELECTORS: raise` to the top of `select`.
- **N4.** Coverage units differ between modules. `pool_metrics.subset_metrics` and `exp5/analysis.py` use fractions
  (`COVP_MARGIN = 0.005`), while `exp3/analysis.py` uses percent (`100 *`, `COV_MARGIN = 0.5`). Both are consistent
  internally; state the unit in each constant's comment or use one convention.
- **N5.** `infer.sbatch:84-86`. The header calls sequential rows in one allocation "packing". plan_P3's packing is
  PACK = 3 *concurrent* processes per GPU. Rename (e.g. "batched rows") so readers do not assume concurrency in the
  timing.
- **N6.** `exp5/analysis.py:154-172`. The validation decision is hashed (`inputs_sha256`), but nothing records that the
  test pools (generated and scored in w1, `stageA_tasks.tsv`) were not read before it. Write the decision file's mtime
  and the test `eval.pkl` mtimes into `exp5_summary.txt` as a procedural audit (the rule itself is mechanical, so this
  is documentation only).
- **N7.** `reverse.py:57-63`. Each numpy ↔ torch round-trip per step reshapes with `n_conf`. Add an assertion that
  `perturb.numel() == n_conf * m` with a molecule-key message, because `edge_pred.squeeze()` in TD
  (`score_model.py:158`) yields a 0-d tensor for one conformer with one rotatable bond.
- **N8.** `fx_common.sh:34-35`. `df -BG … | tr -dc 0-9` on a GNU `df` that prints "50G" works, but a missing `$RES`
  makes `free_gb` empty and `(( > 50 ))` a syntax error under `set -e`. Guard with `${free_gb:-0}`.
- **N9.** `stageA_tasks.tsv`. QM9 `F_flip` uses `keys_double.txt` while DRUGS `F_flip` uses
  `keys_drugs_double_or_triple.txt` (its O's subset). Both are right for their base, but a one-line comment in
  `make_tasks.py:41-45` would save the next reader the comparison.

## Answers to the coordinator's specific questions
- **One loop, extensible for milestone 3?** It is the one loop, golden-equal to `sample()`, and PG is not copied. It is
  not yet extensible without editing; S1 gives the minimal protocol change.
- **Duplication of TD code outside `td.py`?** Every TD symbol is imported through `import_baseline`. The remaining
  re-implementations are each justified and tested:
  - `embed` and `seed_stage` (per-conformer seeds, which TD's single global seed cannot give);
  - the reverse loop (golden-tested);
  - `torsion_types` (ported in milestone 1 with cluster goldens).

  No unjustified copy was found.
- **DRY_RUN 24.95 GPU-h: right?** The arithmetic is right: w1 4.87 + w2 4.785 + w3 4.10 = 13.76, plus 11.2 = 24.96.
  But it is generation only (B2), it reserves nothing for Exp 1 at submission (B3), it uses the M1 TD-speed figure
  rather than the 16 GPU-h ceiling (≈ 30.3 GPU-h at the ceiling), and the DRUGS rows are unmeasured (S3). As an
  allocation estimate it is low by several GPU-h.
- **Cut order if real timings run over?** DECISION item 10 as written (pack M1 → drop P6 → reduce PG tuning → stop and
  ask), with the caveat that "drop P6" saves 0 GPU-h here. S8 gives the concrete order with savings.
