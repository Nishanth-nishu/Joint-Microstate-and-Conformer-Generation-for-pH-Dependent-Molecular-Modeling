# TorDiff-NAX: FlexiTors-Diffusion research programme

Handoff document for a new chat session. State as of **24 October 2026**. Read this first, then `EXPERIMENTS.md`.

The goal is a publishable paper extending **Torsional Diffusion** (TD; Jing, Corso et al., NeurIPS 2022) beyond
torsion-only sampling. TD learns torsions τ but copies local structure L (bond lengths, angles, ring geometry) from
RDKit. We are testing whether learning L, coupled to τ, fixes TD's error ("FlexiTors-Diffusion"). Phase 1 is QM9;
GEOM-DRUGS comes next. The professor asked for this direction; the user (GitHub `Nishanth-nishu`) runs the project.

---

## 1. Ground rules (user decisions; follow them)

- **Repo:** GitHub `Nishanth-nishu/TorDiff-NAX`, branch `main`. The user is the **sole author**: commit as
  `Nishanth-nishu <71981689+Nishanth-nishu@users.noreply.github.com>` (repo-local git config is set). **Never add
  Co-Authored-By or any co-author line.** History was rewritten on 7–8 Oct to remove other contributors; pre-rewrite
  backups are in `.nested_git_backup/` (not pushed).
- **License:** `torsional-diffusion/LICENSE` (MIT, Jing & Corso) must stay.
- **Line endings:** LF (`.gitattributes`); bash/sbatch scripts break with CRLF.
- **Not in git:** paper PDFs and extracted fulltext (copyright), blog snapshots, `*.pkl`, `*.pt`, large CSVs.
- **Primary endpoint (pre-registered):** AMR-R mean and COV-R at δ = 0.5 Å for QM9 (0.75 Å for GEOM-DRUGS), paired
  molecule bootstrap, Holm correction. Key secondary: COV-R@0.1 on the success intersection. Do not change without
  asking.
- **GPU spending:** ask the user before GPU-heavy runs. Current Stage A cap: **45 GPU-h** (user, 17 Oct), with
  Experiment 1's gate (≤ 16 GPU-h) reserved first.
- **Downloads** need the user's approval (GEOM-DRUGS was approved on 16 Oct).
- **Honesty:** report corrections openly. Several earlier conclusions were overturned (§5).

## 2. Cluster

- SSH alias `ada` → `ada-gw1.iiit.ac.in`, user `nishanth.r`, key `~/.ssh/ada_iiit_new` (the old host
  `ada.iiit.ac.in` is gone). **Never log in with a password.**
- SLURM account and partition `plafnet2`, node **gnode118** (4× RTX 3090, 48 CPUs shared with other users, ~250 GB
  RAM). **Memory cap 5000 MB per CPU**: larger `--mem` silently raises the CPU count. Max job time 4 days.
- `/scratch/nishanth.r/tordiff` (3.5 TB) exists **only on gnode118**: reach it via
  `srun -n 1 -c 1 --mem=2G -t 00:10:00 -A plafnet2 -p plafnet2 -w gnode118 bash -c "..."` or sbatch. Always `-n 1`.
  Home (`/home2`) has a 30 GB quota.
- Layout on scratch: `torsional-diffusion/` (baseline), `flexitors/`, `venv/` (py3.9, torch 1.13.1, RDKit 2022.09.5),
  `data/QM9/`, `workdir/<run>/` (checkpoints), `results/<run>/<tag>/`, `logs/`. Paths are in `slurm/common.sh`.
- xtb 6.7.1 is available as the environment module `u22/xtb/6.7.1` on gnode118 (~0.14 s per QM9 molecule).
- Known cluster hazards: GPU out-of-memory from co-located jobs (twice; rerun the task); long-time-limit jobs wait
  hours behind other users' jobs (shorter limits backfill faster); round-1 TD code is not bit-reproducible run to run
  (COV-R 89.0286 vs 89.0293), and GPU sampling is non-deterministic (0.265 Å per-conformer deviation between reruns).

## 3. Repository map

| Path | What it is |
|---|---|
| `torsional-diffusion/` | TD baseline (upstream 5f713b4) + round-1/2 hooks (all default-off flags). **Frozen**: SHA-pinned by tests |
| `flexitors/` | New clean package (round 5): `td.py` (single import gateway to TD), `geometry/`, `models/`, `sampling/`, `training/`, `eval/` (incl. `validity.py`), `data/`, `tests/` |
| `experiments/` | Round-5 job scripts and READMEs: `common/`, `exp3/`, `exp5/`, `r5_cpu/` |
| `tools/` | Analysis tools (paired_compare, analyze_round2, analyze_leakfree, make_split_half, local_structure_analysis, lgeom, …) + `tools/tests/` |
| `slurm/` | Round-1/2/4 SLURM scripts, `common.sh` |
| `papers/` | Literature store: `INDEX.md` (index), `core/` (TD paper + `tordiff_analysis.md`), `related/` (PDFs + fulltext, local only), `gap_synthesis.md` |
| `orchestration/` | Multi-agent protocol v2: `README.md`, `roles/`, `templates/` |
| `round2/` … `round5/` | Each round's brief, research, verification, plans, votes, decisions, implementation records |
| `cluster_sync/` | Small result files copied from the cluster (summaries, paired tables, analyses) |
| `EXPERIMENTS.md` | **Every experiment: why, source, code change, outcome, interpretation, worth** |
| `RESULTS_QM9.md`, `PROGRESS_REPORT.md` | Round-1 reports (early 3 Oct; partly superseded) |
| `reports/round2_report.html` | Visual report (published artifact; round 1–2 era) |

Shared docs: experiments doc with charts `https://claude.ai/code/artifact/f6f439b5-684e-4087-9d40-9e6387f87950`;
visual report `https://claude.ai/artifact/D8fPvHneWhYtMKUMgEJZxn` (both private to the user).

## 4. How the work is organised (multi-agent protocol, `orchestration/README.md`)

The user wants research and coding done by **many agents that cross-check each other**:

1. **Research agents (scouts)** find evidence in papers and blogs. Every claim is an evidence entry with PDF path, PDF
   page and exact quote (≤ 40 words), or `path:line`, or result file + value. Papers are downloaded to `papers/related/`
   with fulltext and an `INDEX.md` row. Ideas become **proposal cards** or **experiment specs**.
2. **Verifiers** independently open every source and recompute every number; disputes go back to the author, who fixes
   or withdraws; the verifier closes each thread.
3. **Impact analyst / code grounder** (round 3) predict effects from our data and map ideas to code.
4. **Novelty red-team** rules NOVEL / INCREMENTAL / PRIOR ART.
5. **Three coding planners** plan independently with different emphases (minimal-clean, validity, cluster
   throughput), ask the research agents questions through `round<N>/consult/Q-*.md`, cross-review and **vote**.
6. **One implementer** builds the winning plan in milestones; **two reviewers** review each milestone; the research
   author of each spec does a **conformance check** (spec quote beside `path:line`); only then the orchestrator pushes.

Practical lessons: agent usage limits hit often (session and weekly). Run **≤ 2–3 agents at a time**, tell agents to
write files incrementally, and resume interrupted agents rather than restarting. Fresh agents are cheaper than resuming
agents with very large contexts.

## 5. What we know (verified) and what was overturned

**Round 1 (30 Sep–3 Oct), QM9, released TD checkpoint and retrains (3 seeds), paired statistics**
- TD reproduced: AMR-R **0.176 Å** (paper 0.178). The paper's 92.8 % COV-R cannot be reproduced even from its own
  released conformers (88.8 %): ~65–69 molecules fail RDKit embedding and count as 0.
- MMFF relaxation of RDKit geometry before diffusion: **0.150 Å** (−15 %), no retraining. Relaxing after: worse recall.
- Stereo from the true structure (G1): no effect. Conformer-matching variants (B2, B5): slightly worse.
- B1 (TD trained on true geometry): 0.234 Å on RDKit geometry (collapse).
- Run-matched floor: up to ~72 % of per-molecule error needs a geometry change (upper bound). Ring geometry dominates
  the floor.

**Round 2 (6–12 Oct): training-robustness arms, final pre-declared paired analysis**
(`cluster_sync/round2/results/analysis_r2/`)
- S3 (50/50 RDKit/true-L training): **+0.0032 Å worse** than control on RDKit L (significant); worse on MMFF L
  (+0.039). S2 (noise): +0.025 worse. S4 (λ-conditioned): +0.003 n.s. S5 (MMFF-matched training, MMFF at test):
  −0.004 vs MMFF at test alone. **No learned training recipe improved real use; only better geometry (MMFF) did.**

**Round 4 (11–15 Oct): leak-free split-half oracle** (`cluster_sync/round4/analysis_r4_leakfree/`)
- Giving a model its **own** target conformer's true geometry leaks conformer identity. Own-conformer oracle gains
  (B1 0.030, S3 0.053 Å) were leakage. **A1's −54 % and the ring-only 0.119 Å are upper bounds that include leakage.**
- The split-half test (one half of a molecule's true conformers supplies geometry, the other half is scored) is a
  candidate methodological contribution; no published paper defines one.

**Round 5 corrections (16–17 Oct)**
- Round 4's "another conformer's true geometry is worse than RDKit for every model" was a **sampling-budget artefact**.
  Budget-matched re-score (711 molecules): standard model **−0.038 Å better** (−0.015 on distinct targets), B1 −0.017
  (n.s.), S3 +0.017 (worse), S4 −0.007 (n.s.). Geometry–torsion **coupling is a hypothesis, not a finding**.
- Bug fixed: `tools/local_structure_analysis.py:302` counted the target among "other" conformers for molecules with
  ≤ 3 conformers (biased floor 0.032 → corrected 0.036 Å vs 0.114 Å for RDKit geometry).
- Stereo census: TD randomises every double-bond torsion and has no E/Z input (`featurization.py:101`,
  `sampling.py:107`); 11.1 % of QM9 conformers (3-seed mean) have an E/Z state absent from the true ensemble; double
  bonds twisted 16.6° on average (true 1.8°). GEOM-DRUGS: 140 test molecules with stereo double bonds (126 specified);
  QM9 only 37 (3 specified) → **DRUGS is the primary set for Experiment 3**.

**Position vs literature:** our non-oracle results (0.150–0.182 Å) sit at OMEGA level; published learned models report
0.073–0.103 Å on the same protocol. Evaluators differ: TD counts failures as 0 coverage, ET-Flow drops them (~3.8 COV-R
points on QM9). TD's sample count 2K uses the test-set conformer count K.

## 6. Round 5: the paper programme (current)

User's six-phase plan (`round5/PLAN.md`): 1 reproduce ✓, 2 error source ✓ (with corrections), 3 smallest viable
extension, 4 controlled ablations, 5 generalisation (strata, OOD, DRUGS), 6 evidence-based novelty.

Experiment specs (all verified by verifier V5; `round5/specs/`, ledgers in `round5/ledger/`):

| Exp | Spec | Status | Novelty verdict |
|---|---|---|---|
| 1 Joint local-geometry + torsion diffusion (V1 fixed L / V2 τ-conditioned L correction / V3 joint) | `EXP-1.md` rev 2 (§10 binding decisions) | Supported; first gate M1 ≈ 10.5 GPU-h (≤ 16); full table ≈ 300–350 GPU-h | Incremental (RINGER, GO-Flow, FoldingDiff prior art); new: coupling evidence, TD-protocol benchmark, projection construction |
| 2 Matching-free training | `EXP-2.md` | Supported; waits for Exp 1 gate | Incremental at best |
| 3 Stereo-constrained torsional diffusion | `EXP-3.md` (§9 binding) | Supported; DRUGS primary | Incremental method; the E/Z error measurement looks unpublished |
| 4 Full Boltzmann generator p(L,τ) | `EXP-4.md` | Feasibility section only (≈ 5–14 GPU-h) | Novel only as transferable p(L,τ) for drug-like molecules |
| 5 Adaptive sampling | `EXP-5.md` (§12 binding) | **Weakened**: oracle allocator beats best static rule by only 0.004 Å; adaptive allocation dropped; **selection** kept (0.178 → 0.170 Å at 2× proposals, needs matched-compute test) | Incremental |
| Shared evaluation | `EVAL.md` | Supported | — |

Recommended paper framing (not yet decided by the user): *"What limits torsional diffusion, and how to fix it"* —
analysis half (leak-free oracle, stereo census, protocol artefacts, corrected floor) is largely evidenced; method half
is Exp 1 (+ Exp 3 constraint). Main venue only if Exp 1 beats xTB-relaxed seeds and the A6b control in real use.

**Stage A (approved 16 Oct, cap raised to 45 GPU-h on 17 Oct):** build `flexitors/`; Exp 1 smoke + Gate 0 (M1);
Exp 3 (QM9 arms on **CPU** with their own baseline O for exact pairing; DRUGS on GPU); Exp 5 selection at matched
compute; CPU: budget-matched leak-free re-score (done), floor-bug fix (done).

Decisions: `round5/DECISION.md` (Stage A scope, plan vote — base `plans/plan_P2.md` + 10 grafts, budget revision).

### Implementation status
- **Milestone 1 — done, reviewed, conformance-checked, pushed:** `flexitors/` skeleton, `td.py` gateway + SHA pins,
  dead-code test, `eval/validity.py` (exact dihedrals, geometric E/Z and chirality labels, connectivity calibrated on
  GT: 99.91 % valid), floor-bug fix with regression test, QM9 validation set (1000 molecules, reproduces 99.1 % of
  test-set construction), budget-matched re-score. Record: `round5/IMPLEMENTATION.md`.
- **Milestone 2 — code committed, reviewed, review fixes NOT yet done:** Exp 3 sampler (one TD reverse loop with
  hooks, golden-equal to TD `sample()`; arms O, K1, K1-wrong, K1-canonical, K2-inf, F-filter, F-flip), Exp 3 analyses,
  Exp 5 arms P0–P6 and selectors, job machinery (OOM retry, resume, GPU-h ledger, STATUS.md), pre-declared analysis
  scripts. 20-molecule GPU smoke passed (0.23 GPU-h). Tests: 89 local / 99 cluster.
  **Open review findings** (`round5/review_P1_m2.md`, `round5/review_P2_m2.md`) — must be fixed before any full run:
  - BLOCKERS: DRUGS scored at 0.5 Å instead of 0.75 Å (`exp3/analysis.py:46`, `infer.sbatch:149` passes
    `--dataset qm9`); GPU-h estimate omits in-job evaluation (+3–8 GPU-h) — redo DRY_RUN; ledger must reserve Exp 1
    M1 first and check cumulative totals; `generate.py:83-98` drops a whole molecule when one ETKDG embedding fails
    (biases arm comparisons).
  - SHOULD-FIX: generic loop hooks so the Exp 1 geometry step attaches without editing the loop; per-arm affected
    sets for splicing; failure reporting and correct N; GT-matching cap; hierarchical bootstrap for K3
    non-inferiority; check flexitors O vs TD released runs; finite validity-adjusted AMR-R; heavy-atom selector
    torsions; ETKDG terminal-alkyne C–H artefact (find root cause; validity-adjusted mask should count only defects the
    method introduced; R-SAMPLE signs off); concrete GPU-h cut-order levers.
- **Milestone 3 — not started:** Exp 1 M1 per `EXP-1.md` §10 (geometry head predicting endpoint x₁, heavy-atom-first
  torsion definitions with exact torsion restoration after each geometry step, chirality guard, A1+ and A6b
  compute-matched controls, Gate 0 decided on validation).

### Cluster state (24 Oct)
- Gate-0 baseline job 24587, floor re-runs 24588/24589: completed 17 Oct (outputs under `results/round5/` on scratch;
  floor summary copied to `cluster_sync/round5/floor/`; **Gate-0 outputs not yet read or synced**).
- GEOM-DRUGS download (job 24342): model, parameters and generated conformers OK; **`drugs.tar.gz` is truncated
  (18.35 of 25.31 GB)**. Resume with the same script (`slurm/download_drugs.sbatch` on scratch; curl `-C -` resumes),
  then check the size and sha256 before extracting to `data/DRUGS/`.
- No jobs running.

## 7. Next steps (in order)

1. Fix the milestone-2 review findings (§6 list), re-run tests, redo DRY_RUN with evaluation included; R-CHEM and
   R-SAMPLE conformance checks; push.
2. Resume the GEOM-DRUGS download, verify, extract; read and sync the Gate-0 baseline outputs.
3. Build milestone 3 (Exp 1 M1) with the same review + conformance gates.
4. Run Stage A within 45 GPU-h (Exp 1 gate reserved first), then the pre-declared analyses.
5. Update `EXPERIMENTS.md`, the docs, and decide with the user on the paper framing and Stage B (Exp 1 full table,
   Exp 2, Exp 4 feasibility, DRUGS generalisation).

## 8. Suggested first message for a new chat

> Read `README.md` and `EXPERIMENTS.md` in `C:\Users\HP\Desktop\tor_diff`, then `round5/DECISION.md`,
> `round5/IMPLEMENTATION.md`, `round5/review_P1_m2.md` and `round5/review_P2_m2.md`. Continue from README §7 step 1
> using the multi-agent protocol in `orchestration/README.md` (≤ 2–3 agents at a time). Commit as the repo author with
> no co-author lines.
