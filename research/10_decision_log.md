# Decision log (UTC)

- 2026-10-26 08:18 — Phase 0: TorDiff-NAX copied with fresh single-author history (user asked to remove all contributors; no co-author trailers).
- 2026-10-26 08:20 — scholar-scout skill unavailable on the Linux login node (path given is a Windows path). Substituted OpenAlex (reference-lookup MCP), web search, WebFetch. Google Scholar recorded as not searched for every family.
- 2026-10-26 08:25 — Baseline frozen at /scratch/nishanth.r/tordiff; new work goes in /scratch/nishanth.r/jmcg only.
- 2026-10-26 — Phase I finding: TD's torsion set is state-dependent (utils/torsion.py:9-38 rotatable iff the cut fragment has >1 atom) → p(z,τ|G,h) is trans-dimensional; joint model must use per-state torsion sets or a heavy-atom torsion parameterisation.
- 2026-10-26 ~08:50 — Literature round 1 complete (157 queries, 78 registry papers: 27 full text, 48 abstract-only, 3 inaccessible). No PDF text extractor on login node → full text read via arXiv/PMC HTML. Citation expansion thin (4 citation_neighbours calls) — must be extended in Phase III/VIII.
- 2026-10-26 — AcepKa (arXiv 2604.00841) identified as the closest system: same interface (SMILES+pH → weighted microstates + conformers), sequential. Becomes mandatory Baseline A.
- 2026-10-26 — Task B (p(z|G,h)) judged largely solved by 2D tools (Epik 7, Uni-pKa); task C has no experimental supervision. Initial direction shifted from "joint generator" (I1) toward geometry-aware microstate populations on a conformation-sensitive stratum (I2/I5) + benchmark (I4). Provisional, pending Phase VIII.
- 2026-10-26 — LoQI audit based on code + dataset docs; paper PDF inaccessible (ChemRxiv Cloudflare 403). Must re-attempt full text before citing methodological details from the paper.
- 2026-10-26 ~11:30 — Phase VIII round 1: independent judge verdicts. I2 Reject (Rowan pKa = same mechanism with ETKDG); I5 Revise (Uni-pKa code already 3D-conformer-conditioned + exp fine-tuned; Starling LSE); I4 Revise (partly done: AcepKa builds pairs but no stratified/state-level evaluation). Merged into narrowed claim N1 (architecture_candidates/N1_intent.md): 2D-vs-3D-input and conformer-source ablation of Uni-pKa-class model, stratified by conformational sensitivity. I1 (full joint generator) remains parked — no supervision; only revisit if N1 shows geometry matters on the sensitive stratum.
- 2026-10-26 — Pre-registered strata definitions fixed in jmcg/data/build_sensitive_set.py before any model run.
- 2026-10-26 10:49 UTC (16:19 IST) — Step 2 built (job 40452): jmcg venv py3.10.12 (RDKit 2026.09.1) on gnode118 from offline wheels in /home2/nishanth.r/wheels_jmcg.
- 2026-10-26 10:52 UTC — Stratum v1 put ~75 % of IUPAC in 'sensitive' → replaced by stricter through-space v2 BEFORE any model output (job 40455); PREREG amended. Clean (non-Dwar-iBonD) entries: sensitive 2,229 / control 678.
- 2026-10-26 — Uni-pKa weights: bindresearch/pkasso_unipka @ e28a591 (Apache-2.0; model card says Uni-pKa by Luo et al.), 5 folds SHA-256 verified against unipkainfer 0.2.0 registry. unipkainfer gives per-conformer free energies (10 ETKDG+MMFF + 1 2D copy) → all N1 arms from one pass. Pilot driver + analysis written (jmcg/experiments/n1_{pilot,analyze}.py). Submission of the pilot job was blocked by the session's permission classifier (installs/runs third-party PyPI package + HF weights) — awaiting user approval.
- 2026-10-26 19:52 UTC — User approved installing unipkainfer 0.2.0 + running SHA-verified Uni-pKa weights (CPU only) for the N1 pilot.
- 2026-10-26 22:24 UTC (03:54 IST) — N1 pilot run 1 (job 41159) FAILED the pre-registered sanity gate (MAE ≈ 6.4). No conclusion drawn. Causes: missing pKa normalisation mean (+6.5049) in my conversion; unipkainfer's free-energy path uses EmbedMultipleConfs(randomSeed=0) → 10 identical conformers (RDKit 2024.03 & 2026.09; jobs 41187/41189/41192). Fixed in the driver by porting Uni-pKa's original seed-loop protocol; added a conformer-diversity guard. Note: in that protocol seeds 0 and 1 yield identical embeddings (faithfully reproduced). Run 2 = job 41194.
- 2026-10-26 — Possible upstream bug report to bindresearch/unipkainfer (identical conformers) — not filed; needs user decision (outward-facing).
- 2026-10-26 23:58 UTC (2026-10-27 05:28 IST) — Run 2 (job 41194) complete: 356/356 transitions, conformer guard passed. Gate failed only on "Novartis" (0.654 vs 0.810).
- 2026-10-27 00:05 UTC (05:35 IST) — Diagnosed: gate references came from the paper's Conclusions sentence, which mislabels Table 1 (0.810 = Novartis ACID; 0.493 = Novartis base; 0.554 = SAMPL6). Against Table 1 every column reproduces within 0.06. Post-hoc amendment 3 logged; re-analysis job 41255 (no model calls).
- 2026-10-27 00:05 UTC — **N1 FALSIFIED (pre-registered rule).** Sensitive Δ(2D−3D) = +0.005 [−0.014, 0.024], n=185; DiD +0.010 [−0.030, 0.051]. LSE/min/single-conformer aggregation identical to mean within 0.003. Released Uni-pKa is conformer-invariant (median per-conformer pKa SD ≈ 0.03). Consequence: I1 stays parked — no evidence from this model class that geometry carries pKa signal; stage 2 (IUPAC) not triggered (not underpowered). Next direction needs a user decision (see 00_project_state.md).
- 2026-10-27 ~10:15 UTC (15:45 IST) — User chose the pivot to tautomer populations (novelty loop + pre-registered pilot).
- 2026-10-27 10:28 UTC — Literature family F_tautomer (62 ops, subagent; Google Scholar not searched; ACS/ChemRxiv 403). Closest: Zankov 2019 JCIM conjugated pKa/logK_T QSPR (abstract only); Gunner 2020 + SAMPL6 overview (macro-pKa does not reveal tautomer ΔG). No per-microstate LSE model trained on pKa + tautomer ratios found.
- 2026-10-27 10:36 UTC — Independent judge: T1 B (diagnostic only; Uni-pKa authors already note "similar free energies to different tautomers"), T2 B (primary comparison must be joint vs tautomer-only, not vs pKa-only; sPhysNet-Taut on same folds needed before any paper claim). Both narrowed claims adopted verbatim.
- 2026-10-27 ~10:45 UTC — Design fix before any output: tautomer-only arm given the same step count as the joint arm (my draft gave it ~12× fewer steps — would have confounded the primary comparison).
- 2026-10-27 ~10:50 UTC — T pilot pre-registered (experiments/T_pilot/PREREG.md). Code: jmcg/experiments/taut_{zeroshot,features,headfit,analyze}.py, slurm_jmcg/t_pilot.sbatch. Not submitted: needs user approval (re-runs unipkainfer + released weights in a new context).
- 2026-10-27 10:46 UTC — User approved running the tautomer pilot (unipkainfer 0.2.0 + SHA-verified Uni-pKa weights, CPU only, gnode118).
- 2026-10-27 ~17:10 UTC (22:40 IST) — User asked to start explicit-state conformer generation on the existing TD implementation. New package jmcg/condgen (states.py, td_backend.py, validate.py, run.py, inventory.py, test_states.py), separate from N1/T pilot code and outputs. TD checkout, venv and checkpoints used read-only via sys.path (no file in /scratch/nishanth.r/tordiff touched). Deviation from TD generate_confs.py recorded: ETKDG seed embedding gets randomSeed (upstream unseeded). Smoke test 1 = job 43186 (CPU, 24 states incl. 2 negative controls, 8 confs, seed 0, + ETKDGv3 baseline + stereo-double-bond-locked TD arm).
- 2026-10-27 18:30–19:30 UTC (2026-10-28 00:00–01:00 IST) — Code audit of jmcg/condgen before any output existed (job 43186 never started; gnode118 fully allocated, est. start 2026-10-28 09:28 UTC). Fixes: (1) stereo policy "first" silently picked a stereoisomer -> removed; now "require" (reject) or "enumerate" (each stereoisomer an explicit state); (2) parent/state matching accepted redox changes (e.g. CCO -> CC=O) and charge-separated resonance forms, and labelled zwitterions "protonation" -> proton-transfer conservation (total dH == dq; per site dH == dq or dq == 0; charge changes on heteroatoms only) + classes tautomer / tautomer_zwitterionic / protonation / tautomer+protonation; (3) parent stereocentres now checked for configuration (mapped neighbour parity, not CIP); (4) SDF check now per atom (Z, charge, H) and per bond (type) in original order + isomeric SMILES from 3D; (5) bond check replaced (covalent radii +-25 % -> RDKit DG bounds +-0.10 A for 1-2 and +-0.20 A for 1-3); (6) H-defined imine E/Z would have produced false stereo failures -> heavy-atom view with removeDefiningBondStereo, limitation documented; (7) unchecked connectivity no longer counts as valid; (8) states tagged A/B/C/D (neutral in-dist / charged in vocab+size / size-or-charge OOD / unsupported). Added test_backend.py (featuriser width, charge one-hot, compatibility, torsion counts, seed reproducibility, E/Z lock, SDF round trip, NaN/stretch/break/clash/mirror/E-Z-flip detection). 43186 cancelled unstarted; replaced by 43422 (same resources; runs tests first). No result is claimed until the log is inspected.
- CSV counts (verified by script): 24 rows = 22 states of 10 distinct parents (10 states identical to their parent) + 2 negative controls (2 more parents).
- 2026-10-27 21:06 UTC (2026-10-28 02:36 IST) — T pilot job 43158 COMPLETED (stages 4–6, 25,681 microstates, 282,491 conformer rows). Sanity gate passed on all six reference columns. **T1 SUPPORTED** (zero-shot compression: M pairs n=221, MAE 3.21 log K, slope 0.10 [0.08, 0.12], predict-zero MAE 3.53). **T2 FALSIFIED** under the pre-registered rule. Out-of-fold tautomer MAE: R-taut 1.905 vs R-joint 1.912, Δ(taut − joint) = −0.007 [−0.118, 0.112]. External macro-pKa got worse with joint fitting: joint − pKa-only = +0.22 [+0.16, +0.28]. Joint MAE 1.91 > 1.0. R-joint+E 1.78 (Δ joint − joint+E = +0.13 [+0.04, +0.22]) is exploratory only, not pre-registered as primary. Full JSON: logs/jmcg_tpilot_s4_43158.log and experiments/T_pilot/run1/analysis.json.
- 2026-10-27 21:06–21:22 UTC (2026-10-28 02:36–02:52 IST) — Conditional generation, smoke 1 (job 43422) and follow-ups 43714 (diagnostic), 43717 (smoke 2), 43724 (tests). Two bugs found in outputs:
  - (a) `EmbedMultipleConfs(randomSeed=0)` gave 8 identical conformers for the ETKDG baseline and for the TD seeds. This repeats the N1 pitfall from 2026-10-26. Run 1's ETKDG comparison and diversity numbers are invalid. Fixed with per-conformer seeds from `numpy.default_rng(seed)`, plus a regression test.
  - (b) The clash rule flagged intramolecular H-bonds. A polar H···N/O pair is now flagged only below 1.5 Å, with a test. 5 verdicts changed on re-validation (3 in run 1, 2 in the seed-1 repro run), all polar-H···O.
  - Final tests: 20/20 + 10/10.
  - Smoke 2, valid/generated conformers: ETKDG 184/184; TD by category A 94/96, B 63/64, C 14/16; lock arm 7/8. Every failure is a clash or an E/Z flip, explained per conformer in experiments/CG_smoke/RESULTS.md.
  - Cross-process TD drift ≤ 1.1e-2 Å (≤ 0.11 Å with the old seeding) with identical verdicts. Not caused by threads or hash seed; cause not verified.
  - QM9 train charge inventory: 111/106,041 charged, all nitro, 0 net-charged → B states are representable but not in distribution.
  - The run used `pretrained/qm9_1order`, not `qm9_default`; the checkpoint must be pinned by path.
  - Joint architecture not started. Proposed next step: a wrapper-fidelity check against TD's generate_confs.py on about 50 GEOM-QM9 test molecules.

## 2026-10-28 — CG wrapper-fidelity check: FAITHFUL
- Job 44897 (10:41–11:00 UTC / 16:11–16:30 IST): TD generate_confs.py vs our TDGenerator on 50 GEOM-QM9 test molecules, qm9_default pinned (sha256 6a9a6ee9…), seeds 0/1, scored with TD evaluate_confs.py.
- Pre-registered rule (jmcg/condgen/fidelity.py) passed: W0−R0 ΔMAT-R −0.010 Å, W1−R1 −0.002 Å vs noise R1−R0 +0.011 Å (tolerance 0.031 Å); ΔCOV-R@0.5 +2.7 / +1.2 pp (tolerance 8.2 pp). Same 4 failures in all arms (3 cage embeddings, 1 corrupted 2-fragment entry).
- Decision: wrapper accepted as faithful for neutral QM9 inputs. Smoke 1/2 TD outputs used qm9_1order and should be re-run with qm9_default before citation. Report: research/experiments/CG_fidelity/RESULTS.md.

## 2026-10-28 — CG smoke 3 (qm9_default pinned)
- Job 45117 (15:56–15:59 UTC / 21:26–21:29 IST). Checkpoint pinned by path and hash. Tests 20/20 + 10/10.
- TD valid: A 94/96, B 63/64, C 16/16, lock 7/8; ETKDG 184/184 (bit-identical to smoke 2). acac Z-enol still flips E/Z under TD (2 of 8), so this is not checkpoint-specific.
- Cross-process TD drift reached 0.49 Å (acac) with identical verdicts, so the 0.11 Å tolerance rule is withdrawn; reproducibility is claimed only as identical verdicts. Smoke 1/2 TD outputs superseded.

## 2026-10-28 — Cross-process TD drift: root cause fixed
- Cause: TD diffusion/torus.py computes score_norm_ by unseeded numpy Monte Carlo at import, so each process scales model outputs by a slightly different table (relative ~1e-5); the SDE amplifies it to 0.03–1.7 Å. Not threads, MKL, ATen ISA, hash seed or ASLR (all tested). Located by a stage-by-stage trace (jmcg/condgen/drift_trace.py) and per-module hooks.
- Fix: td_backend imports TD under np.random.seed(20221024), restores the caller's RNG, refuses a prior unseeded torus import, and records the table hash. Verified bit-identical across processes and thread counts; new regression test (test_backend 11/11, job 45222).
- Reproducibility claims are now "bit-identical for a fixed seed on this node". TD's own generate_confs.py has the same issue (does not change the fidelity verdict). Details: CG_smoke/RESULTS.md section 9.

## 2026-10-28 — CG smoke 4 (fixed implementation, df5210a): passed
- Job 45272 (16:55–17:00 UTC / 22:25–22:30 IST). Provenance verified: HEAD df5210a clean, qm9_default sha 6a9a6ee9…, torus table 8cf17268cc0b0dee in both processes, inputs sha db070079….
- All 24 inputs accounted for. Valid/generated: ETKDGv3 A 96/96, B 64/64, C 16/16, D (cysteine) 8/8; TD A 94/96, B 63/64, C 16/16, lock 7/8; wrong-skeleton control rejected; cysteine TD hard-incompatible (S).
- Whole run bit-identical in a second process (46 SDFs, max |dx| 0); SDF re-read changed 0 verdicts.
- Decision: the smoke pipeline is validated. No quality claim, TD vs ETKDG. Open: E/Z lock is opt-in (proposed default for stereo-specified double bonds); DRUGS checkpoint absent. Details: CG_smoke/RESULTS.md section 10.

## 2026-10-28 — E/Z lock made the run default (committed 2026-10-29)
- run.py --ez_policy lock (default): td locks only double bonds whose E/Z the input SMILES specifies; enumeration-assigned and unspecified bonds stay free; --ez_policy free is the ablation opt-out; per-state reason in <id>.json ez_lock. Backend default unchanged (fidelity check still applies).
- Smoke 5 (job 45286): tests 20/20 + 15/15 (4 new lock regressions); 46/46 outputs bit-identical to smoke 4 with arms mapped (td = old td_lockdb). TD A 95/96, B 63/64, C 16/16. The lock preserves stereo only; 2 clash failures remain (acac 1.774 Å, histamine dication 1.345 Å).

## 2026-10-28 — M2–M5: canonical data, sequential baseline, E1/E6, E2 joint prototype
- **M2 canonical_v1 (job 45295):**
  - 15,502 observations, 29,280 states, 11,355 parents; 3 records dropped, with reasons.
  - Leakage L1/L1b/L2 = 0. 18 val and 18 test near-duplicates are flagged, not removed.
- **M3/M4 E1 (jobs 45306, 45308):**
  - External macro-pKa MAE 0.633 [0.56, 0.71] on oracle states; the full template matches it with 0 failures; the
    simple template fails 75/356.
  - N1 reproduced 356/356 (max |Δ| 1e-9) once a duplicate idx in N1's file was keyed around.
  - Tautomers: zero-shot measured MAE 3.21, slope 0.10. The enumerator never proposes the partner tautomer (0/221).
  - Conformers on neutral GEOM-QM9: TD beats ETKDGv3, COV-R@0.5 +6.4 pp, MAT-R −0.06 Å.
- **E2/E3 (job 45572):**
  - Pre-registered verdict SUPPORTED: Δ_S 0.39 [0.14, 0.64], Δ_C 0.23 [0.025, 0.49], on 150/221 eligible measured
    pairs.
  - The gain is carried entirely by φ3 (Boltzmann-weighted polar SASA); φ1 and φ2 contribute nothing.
  - Post hoc (job 45574): Δ_C lower bound > 0 in only 15/20 alternative fold seeds; φ3 vs 2D TPSA head to head CI
    includes 0; canonical test split (n = 19) agrees in direction.
  - Decision: report "a 3D polar-exposure term helps calibrate Uni-pKa tautomer scores", not "conformer coupling
    helps". φ1 is defective (no conformer de-duplication) and must be fixed before any entropy claim.
- **Bugs fixed before the fit:**
  - MMFF typing rewrote aromatic bond types in place, so features v1 were invalid (PREREG amendment 2).
  - Conformer-comparison NaN.
  - N1 join key.
- **M5 sampler (job 45584):** runs end to end; at θ = 0 it equals softmax(fold-mean G_U) to 1e-16.
  - Masking unscorable states gave neutral glycine 100% at pH 7.4; replaced by a flagged θ = 0 level fallback.
  - Open blocker: the validator rejects MMFF minima of zwitterions and carboxylates/phenolates, so the coupled model
    cannot score the dominant physiological states. A design decision for the user.
- **E4 not run:** no pH-resolved data.

## 2026-10-29 — P1 unminimised feature path (exploratory); P2 audit of the modified conformer model and EXP-1 M0
- **P1 (jobs 46561, 46673; `experiments/P1_unmin/RESULTS.md`), exploratory, not pre-registered:**
  - The MMFF path is kept as the control and reproduces E2 features_v2 exactly. Validator thresholds are unchanged.
  - With unminimised ETKDGv3 conformers, the share of states with features rises:
    - carboxylate 0.17 → 0.96;
    - phenolate 0.43 → 0.95;
    - carboxylate zwitterions 0.03 → 0.73;
    - neutral 0.84 → 0.97.
  - Cause measured: MMFF94s C–O⁻ (and C–N⁺, nitro N⁺–O⁻) minima fall outside the validator's distance-geometry
    bounds; the same conformers before minimisation pass.
  - φ3 agrees across paths (r = 0.994); φ2 (IMHB) does not (MMFF creates salt bridges).
  - On the 206 pairs eligible under unmin, J1 shows no gain over S1 (Δ_S 0.03 [−0.25, 0.27]).
  - M5 demo with the unmin path: no θ = 0 fallback levels and no pairs without a conformer.
  - TD covers only 621/1,080 panel states (QM9 vocabulary, state-builder checks).
  - **Decision for the user:** adopt the unmin path (and which φ terms) as the joint-model definition. That needs a
    new pre-registration before any confirmatory use. The E2 result and reports are unchanged.
- **P2 audit (`experiments/P2_conformer_arch/AUDIT.md`):**
  - No modified conformer architecture exists in code: `flexitors/models` and `flexitors/training` are stubs. EXP-1
    V2 was specified, not built.
  - The trained variants (round 1/2) are TD-architecture recipe changes, plus S4's λ input. None improved real use.
  - The E1 reference stays frozen TD `qm9_default` on the 183 common E1 molecules.
- **EXP-1 M0 (job 46664; `experiments/EXP1_M0/RESULTS.md`): PASS by the pre-declared rule.**
  - Validation cache, 30,420 primary pairs: recovery 99.6 %, skip share 4.74 % (limit 5 %), 0 flipped centres.
  - W drift after restoration ≤ 3e-15 rad; ring-dihedral > 10° share 41.7 % → 0.6 %.
  - No model, training run or GPU use. Next: build the V2 geometry head and a 1-epoch smoke, pending the user's GPU
    decision.

## 2026-10-29 — EXP-1 V2 built; training smoke FAIL on S2 (exploratory)
- Built (staging, uncommitted): batched torch projector, paired joint transform, TD model + 1o geometry head,
  projected loss, training loop, τ-frozen read-out, tests (`flexitors/tests/test_v2.py`), smoke checker.
- GPU smoke 46725 never ran (plafnet2 GPU-minute quota exhausted); CPU fallback 46744 + checks 46848.
- Smoke (`experiments/EXP1_V2_smoke/RESULTS.md`): S1, S3, S4, S5 PASS; **S2 FAIL** (val geo 0.01302 vs Δ̂=0 base 0.01288).
- Overfit diagnostic 46855 (32 mols, 400 iters): w_τ 0 reaches 0.856 × base (pass, narrowly); w_τ 1 0.913 (fail).
  Leading hypothesis: loss-scale imbalance (L_geo ≈ L_τ/80). Not tested.
- **Decision:** per the brief, no controlled training or E1 evaluation (smoke did not pass); Gate 0 not started.
  Next step needs the user: which fix to pre-declare (L_geo normalisation / w_geo) and where to run (GPU quota).
- Re-smoke 46998 with `--geo_norm base` (pre-declared, criteria unchanged): **S2 FAIL again** (val geo/base 0.989).
  Head gradient now stays O(1), so loss-scale imbalance is not the main cause. Still no controlled run, Gate 0 or E1.
- Re-smoke 3 (47102), `--l_head_init zero`: **PASS** (val geo/base 0.703; read-out closer 61 %, ring > 10° 56 % → 25 %).
- **Decision: train.** EXP1-V2-T1 pre-declared in `experiments/EXP1_V2_train/DESIGN.md`: geometry-only (w_τ 0), arms
  projector on / off (targeted ablation), 10 epochs on CPU (GPU quota exhausted), validation-selected checkpoints;
  E1 post-processing of frozen-TD samples vs TD, MMFF control and ETKDGv3, primary = MAT-R (TD+L) − TD.
  First submission 47158/47159 hit the 1 h default time limit (no `--time`); resubmitted 47246/47247; E1 job 47249
  queued on afterok.
- EXP1-V2-T1 finished (2026-10-31; `experiments/EXP1_V2_train/RESULTS.md`). Primary: MAT-R (TD + L(proj)) − TD =
  −0.0547 Å [−0.0628, −0.0473] seed 0, −0.0557 seed 1 → **improves on TD** under the pre-declared rule. Beats the MMFF
  control by ≈ 0.045 Å. Projector ablation: proj − noproj −0.0009 [−0.0020, 0.0002] on seed 0 → **projector not shown
  to matter**. Torsions unchanged (drift ≤ 2e-13°), so the gain is local geometry only. Exploratory (one training seed,
  CPU, 10 epochs, not converged).
- Addendum A (pre-declared in `EXP1_V2_train/DESIGN.md`, job 49899, CPU; GPU still blocked): torsion-constrained
  GFN2-xTB of the same TD conformers gives MAT-R 0.097 / 0.089 vs TD + L(proj)fix 0.111 / 0.106; L − xTB = +0.014
  [0.008, 0.019] seed 0, +0.017 seed 1 → **head does not beat xTB**, and is ~4.7× slower on CPU (1.67 vs 0.36 s per
  conformer). Read-out fixed for rigid molecules (0 errors). **Decision (user agreed 2026-10-31): stop the V2 geometry head as a QM9
  contribution**; adopt TD + constrained xTB as the baseline any geometry method must beat.
- EXP1-V2-T2 (2026-10-31; user brief "post-training scientific validation"; design and results in
  `experiments/EXP1_V2_T2/`). Reopens the stopped V2 line for verification and a pre-declared longer, three-seed run.
  Audit (49955): all T1 E1 outputs complete and every reported mean and CI reproduced exactly. OOM (49957): worker
  start-up kills only; re-run bit-identical, nothing truncated. Rigid molecules: root cause is TD's
  `out_nodes or n_atoms` fallback; fix + 3 tests; corrected complete evaluation = T1 fix arms. Geometry decomposition
  (49960): constrained xTB beats the head on bond, angle, ring and acyclic-dihedral error and validity; xTB bond MAE vs
  references 0.0027 Å, consistent with the references being xTB geometries. Training T2 (20 epochs, proj/noproj × seeds
  0–2, CPU) running: 49962–49972, resume jobs 49974–49979.
- EXP1-V2-T2 final (2026-11-04). Six 20-epoch runs completed (last 2026-11-04 03:39 UTC); final evaluation run once
  (57155, 8 CPUs / 6 workers per addendum D, 2 CPUs left free for other users). Seed-averaged MAT-R, 184 paired
  molecules: head − TD −0.063 [−0.072, −0.055] (3 / 3 seeds, both TD seeds); proj − noproj −0.0002 [−0.0013, +0.0011];
  head − constrained xTB +0.009 [+0.004, +0.014]; 20 vs 10 epochs −0.004 [−0.006, −0.002]; head 4.1× xTB CPU cost.
  H1 supported, H2 failed, H3 supported, H4 failed, H5 failed. **Decision (pre-declared rule, addendum B): C — freeze
  FlexiTors V2 as a completed negative result and prioritise joint microstate–conformer generation.** Not a methods
  paper; a benchmark note would need other refiners and GEOM-DRUGS. `experiments/EXP1_V2_T2/RESEARCH_DECISION.md`.
