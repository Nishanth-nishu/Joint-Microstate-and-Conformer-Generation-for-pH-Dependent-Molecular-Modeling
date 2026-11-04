# FlexiTors V2: final report and scientific assessment (2026-11-04)

All numbers come from `FINAL_RESULTS.md`, the single final evaluation run as job 57155. Rules come from DESIGN.md
Part D and addenda B–D, fixed before any T2 test number. Literature is from `LITERATURE_COMPARISON.md`.

Notation:
- "mÅ" = 0.001 Å.
- Intervals are 95 % paired bootstrap CIs over the 184 paired molecules.
- "proj" = projector on, "noproj" = projector off, seed-averaged over training seeds 0–2.
- The primary TD sampling seed is 0; TD sampling seed 1 is the replication.

## 1. Research question
Can a learned, torsion-preserving local-geometry head make a frozen Torsional Diffusion (TD) model's conformers more
accurate on GEOM-QM9? TD fixes bond lengths, angles and ring geometry at RDKit's values. The question is whether the
head beats the physical alternative: a torsion-restrained GFN2-xTB optimisation of the same samples.

## 2. Original TD architecture (unchanged, frozen for sampling)
- **Model:** TD's `TensorProductScoreModel`, an e3nn tensor-product network over the molecular graph with
  4 convolution layers in the QM9 configuration. It has a noise-level (σ) embedding and outputs one torsion score
  per rotatable bond.
- **Sampling:** reverse diffusion on the torus of rotatable-bond torsions. Local structure (bonds, angles, rings)
  comes from an RDKit/ETKDG seed conformer and never changes.
- **E1 samples:** frozen TD samples (checkpoint `qm9_default`, TD sampling seeds 0 and 1) on 200 GEOM-QM9 test
  molecules, budget 2 × n_true.

## 3. Our modifications (FlexiTors V2)
- **Geometry head:** an e3nn `o3.Linear` from the last convolution layer's irreps to `1x1o`. This gives one Cartesian
  vector per atom, which transforms as a vector, so the head is E(3)-equivariant. The last layer's output is captured
  with a forward hook; TD's forward pass is not copied. The head's weights start at zero, so the model initially
  predicts "no change".
- **Geometry clock t ∈ [0, 1]:** fed through TD's existing λ-embedding path (32 dimensions).
- **Torsion constraint:** W is the set of TD's rotatable-bond torsions. Let C = [W; J_rig], where W stacks the
  Blondel–Karplus gradients of the W torsions and J_rig stacks the 6 rigid-motion directions. Then

  P_W = I − Cᵀ (C Cᵀ)⁺ C,  rank 3n − 6 − |W|.

  P_W projects any Cartesian displacement onto the directions that change neither the W torsions (to first order) nor
  the rigid pose.
- **Target (pure L):** x_G* = Kabsch( SetTorsions(x_G, φ_W(x_R)) ). The reference conformer x_G is given the seed's W
  torsions, so the target differs from the seed x_R only in local structure.
- **Bridge (training):**
  - x_t = (1 − t) x_R′ + t x_G*′ + γ(t) P_W z;
  - γ(t) = σ_b √(t(1 − t)), σ_b = 0.15, z ~ N(0, I);
  - a chirality redraw if a stereocentre flips, at most 5 tries.
- **Loss:**
  - L_geo = ‖P_W(x̂₁ − x_G*′)‖² / rank, where x̂₁ is the head's endpoint prediction;
  - objective J = L_geo / sg(L₀), normalised by the stop-gradient base loss L₀ of predicting no change;
  - torsion-loss weight w_τ = 0.
  - The noproj arm is identical except that P_W is replaced by I in the bridge noise, the loss and the read-out.
- **Read-out (inference):** 20 steps.
  - Update: δ = Δt / (1 − t) · P_W(x̂₁ − x).
  - Guard: bonds within ±15 %, no clash below 0.75 × the covalent-radius sum, chirality preserved. A failed step is
    halved up to 3 times, then skipped.
  - At the end, Restore_W resets the W torsions exactly to TD's values (drift ≤ 10⁻¹² °).

  The head therefore never changes TD's torsional modes. It cannot create new conformational modes, and none is
  claimed.

## 4. Training procedure
Settings:
- **Data:** the round-2 paired cache, 82,895 train / 10,376 validation pairs. E1 test molecules were never used.
- **Optimiser:** Adam, lr 1e-3, batch 32, 100 validation iterations with fixed noise.
- **Checkpoint:** lowest validation objective.
- **Compute:** CPU, 8 CPUs per run (GPU quota exhausted).

Runs and budget:
- **T1 (exploratory):** 10 epochs, seed 0.
- **T2 (final):** 20 epochs, proj / noproj × seeds 0, 1, 2. Seed 0 continues T1 exactly, with optimiser, scheduler
  and RNG states restored.
- **T2 cost:** about 1,420 CPU hours.
- **Selection:** every run selected its last epoch, so validation was still falling slowly (about −0.004 objective
  per epoch over the last 5 epochs).

## 5. Experimental design (pre-declared)
- **Arms on identical TD samples:**
  - TD;
  - ETKDGv3;
  - TD + torsion-restrained MMFF94s;
  - TD + torsion-restrained GFN2-xTB (k = 1.0 Eh, torsions then restored);
  - TD + L(proj) / L(noproj), T1 and T2.
- **Metrics:** COV / MAT recall and precision (TD's `evaluate_confs.py`), and geometry MAEs against the nearest
  reference (bond, angle, ring dihedral, acyclic dihedral). Validity: identity, clash, stereo.
- **Endpoints 1–7 and the decision rule:** DESIGN.md Part D and addendum B. The strata and near-duplicate checks are
  exploratory (addendum C).
- **Reference caveat:** GEOM-QM9 reference conformers are GFN2-xTB geometries (xTB's bond MAE to them is 0.0027 Å).
  A method that converges to xTB minima is therefore favoured by construction.

## 6. Results summary
| | MAT-R (TD seed 0 / 1) | vs TD | vs constrained xTB | CPU s / conformer (like for like) |
|---|---|---|---|---|
| TD | 0.170 / 0.165 | — | — | — |
| TD + MMFF (restrained) | 0.163 / 0.152 | | | 0.018 (16 workers, secondary) |
| TD + L(proj), seed-avg | **0.106 / 0.103** | −0.063 [−0.072, −0.055] | +0.009 [+0.004, +0.014] | 1.23 |
| TD + constrained xTB | **0.097 / 0.089** | −0.073 [−0.083, −0.063] | — | 0.30 |

## 7. Statistical analysis (Task D)
**D1. Does the head improve TD across all three training seeds?**
Yes.
- **Per seed:** proj − TD MAT-R is −0.063 [−0.071, −0.055], −0.064 [−0.073, −0.056] and −0.064 [−0.073, −0.056].
- **Size:** about 37 % of TD's MAT-R.
- **Other metrics:** MAT-P improves by −0.076 and COV-R@0.25 by +5.4 points [+3.1, +8.0]. This is practically
  meaningful.

**D2. Is it consistent across the two TD sampling seeds?**
Yes. On seed 1 the differences are −0.062 to −0.063 for each training seed, and −0.0625 [−0.0708, −0.0548]
seed-averaged.

**D3. Does projector on beat projector off consistently?**
No.
- **TD seed 0:** seed-averaged −0.0002 [−0.0013, +0.0011]. Per training seed: +0.0001, −0.0002, −0.0005, every CI
  including 0.
- **TD seed 1:** −0.0010 [−0.0022, −0.0001], and 2 of 3 seed CIs exclude 0.
- **MAT-P:** favours proj by 1.4–1.5 mÅ on both seeds.
- **Geometry:** proj has slightly *worse* bond (+0.0008 Å), angle (+0.13°) and ring-dihedral (+0.19°) error than
  noproj; those CIs exclude 0.

**D4. Is the projector effect practically meaningful?**
No.
- **Size:** the largest estimate (1.0 mÅ, TD seed 1) is 1.6 % of the head's gain over TD. It is 7–11 % of the remaining
  gap to xTB and 0.2 % of the 0.5 Å coverage threshold.
- **Context:** it is the same order as the training-seed SD (0.2–0.5 mÅ).

**D5. Does longer training help?**
Yes, slightly.
- **Effect:** T2 seed 0 − T1 seed 0 (proj) is −3.9 mÅ [−6.1, −2.1] on TD seed 0 and −3.2 mÅ [−5.0, −1.6] on seed 1;
  noproj −5.0 mÅ.
- **Size:** 10 more epochs closed nearly 30 % of the T1 gap to xTB (proj s0 − xTB: +13.8 → +9.8 mÅ).
- **Reproducibility:** the three seeds agree to 0.2–0.5 mÅ SD.
- **Caveat:** validation was still falling at epoch 20, so the head is not converged.

**D6. Does the head approach constrained xTB?**
It is closer than after T1 but still significantly worse.
- **MAT:** MAT-R +9.3 mÅ [+3.6, +14.1] (seed 1: +14.2 [+9.6, +19.9]); MAT-P +16.7 mÅ [+11.3, +22.9].
- **Coverage:** COV-R at 0.25 and 0.5 Å is indistinguishable from xTB (−0.2 [−0.8, +0.5] and +0.3 [−1.0, +1.9]
  points). COV-P@0.25 is 1.2 points lower [−2.3, −0.4].
- **Size:** the MAT-R gap is about 13 % of the head's own gain over TD.

**D7. Does constrained xTB remain superior on local geometry and MAT?**
Yes, on every component.
- **Local geometry** (proj − xTB): bond +0.0088 Å [+0.0083, +0.0093], angle +0.47° [+0.41, +0.53], ring dihedral
  +1.5° [+1.0, +2.0], acyclic dihedral +6.7° [+3.6, +10.2].
- **MAT:** MAT-R and MAT-P favour xTB on both TD seeds.
- **Size:** xTB's bond error (0.0027 Å) is a quarter of the head's (0.0115 Å).

**D8. Are gains concentrated in particular structures?** (exploratory, n ≥ 24 per stratum)
- **Over TD:** gains are largest for rigid molecules (0 rotatable bonds: −81 mÅ), polycyclic molecules (≥ 2 rings:
  −78 mÅ) and molecules with 3- or 4-membered rings (−76 mÅ). They are smallest for flexible (≥ 2 rotatable: −33 mÅ)
  and acyclic molecules (−30 mÅ). This is where TD's frozen RDKit local structure is most wrong.
- **Gap to xTB:** largest for polycyclic and strained molecules (+17 and +14 mÅ).
- **Projector:** the only projector effect excluding 0 on TD seed 0 is in flexible (−2.8 mÅ) and acyclic (−1.9 mÅ)
  molecules. On TD seed 1 the flexible-stratum CI includes 0. These are hypothesis-generating only.

**D9. Sensitive to the stereo-insensitive near-duplicate?**
No. Without `N#CC12C3C4CN1C2C43` (n = 183), on TD seed 0:
- proj − TD is −0.063 to −0.064 for every seed;
- proj − noproj −0.0002 [−0.0014, +0.0012];
- proj − xTB +0.0092 [+0.0036, +0.0140].

**D10. Does the head improve validity or stereochemical preservation?**
It improves validity over TD; it preserves stereochemistry as well as TD's torsions allow.
- **Identity failures:** 136 → 26–30 of 5,388 conformers. Clashes 22 → 6.
- **Stereo:** stereo changed vs the TD input in 2 conformers (0 on TD seed 1). MMFF changed 79.
- **Against xTB:** xTB is better again (identity 3, clash 1, stereo changes 4 / 0).
- **Inherited failures:** the 268–272 stereo-vs-SMILES failures are TD's own and pass through unchanged. The head
  cannot fix wrong stereochemistry in TD's samples, and does not try to.

**Variability** (reported separately, not combined; TD seed 0):

| source | size |
|---|---|
| training-seed SD of mean MAT-R | 0.5 mÅ |
| TD-sampling-seed shift of the mean | about 3 mÅ, the same for every arm |
| molecule-level bootstrap CI width of proj − TD | 17 mÅ |

## 8. Hypotheses (Task E), judged by the pre-declared endpoint rules
The brief's H1–H5 map onto the DESIGN.md Part D endpoints as follows. Each verdict uses that endpoint's rule
unchanged.

| hypothesis | pre-declared test | result | verdict |
|---|---|---|---|
| **H1** The head improves local geometry over unrefined TD | Endpoint 1: all 3 training-seed CIs of MAT-R(proj − TD) < 0 and TD seed 1 same sign. Endpoint 6 describes where. | 3 / 3 CIs < 0, seed 1 same sign. Bond −64 %, angle −65 %, ring dihedral −5.7°. | **Supported** |
| **H2** The projector contributes a meaningful accuracy improvement | Endpoint 3: seed-averaged MAT-R(proj − noproj) CI < 0 and ≥ 2 / 3 seeds same sign | −0.2 mÅ [−1.3, +1.1]; local geometry slightly worse with proj | **Failed** |
| **H3** Longer training and multiple seeds improve reproducibility | Endpoint 5: T2 s0 − T1 s0 paired CI. Endpoint 2: between-seed SD and seed-averaged effect. | −3.9 mÅ [−6.1, −2.1]; SD 0.5 mÅ vs effect 63 mÅ | **Supported** (see note) |
| **H4** The head is competitive with constrained xTB | Endpoint 4: "beats" if CI < 0, "indistinguishable" if it includes 0, else "worse" | +9.3 mÅ [+3.6, +14.1] (seed 1 +14.2) | **Failed ("worse")** |
| **H5** The head offers a computational advantage over constrained xTB | Endpoint 7: like-for-like CPU s / conformer | 1.23 vs 0.30, 4.1× more expensive | **Failed** (on CPU; GPU unmeasured) |

- **Note on H3:** "improves reproducibility" has no pre-declared comparator, since T1 had one seed.
  - What is shown: the result is reproducible across 3 training seeds (SD 0.5 mÅ) and 2 TD sampling seeds, and
    20 epochs beat 10 epochs by a small, significant margin.
  - What is not shown: that longer training reduced seed variance.
- **Note on H4 and the xTB references:**
  - The references are xTB minima, so H4 tests whether the head reproduces an xTB-style local relaxation as well
    as xTB itself does. It does not.
  - Both methods reduce the same error components (bond, angle, ring), but this does not show that the head
    performs the same optimisation. The head is a learned map trained toward those references.
- **Note on numbering:** the interim deck (slide 16) used the earlier T1 numbering (H3 "beats MMFF", H4 "beats xTB",
  H5 "longer training"). On that numbering the verdicts are:
  - beats MMFF: yes, −57 mÅ [−66, −47];
  - beats xTB: no;
  - longer training: yes.

## 9. Comparison with constrained xTB (the key control)
On GEOM-QM9, torsion-restrained GFN2-xTB of TD's own samples is better than the trained head on:
- MAT-R and MAT-P;
- every local geometry component;
- validity.

It is also 4.1× cheaper per conformer on CPU. The head matches xTB only on recall coverage (COV-R), where neither
differs significantly.

The head is therefore dominated on QM9-CPU in every pre-declared respect. This is not evidence about larger
molecules, GPU read-out or other datasets, where xTB's cost grows faster; none of these was testable here.

## 10. Literature (Task G, details in `LITERATURE_COMPARISON.md`)
- **Prior refinement of TD samples:** DiSCO (AAAI 2024) already refines TD samples on GEOM-QM9, among others, with a
  Cartesian bridge on all atoms and no torsion constraint. Its physics comparison is unconstrained MMFF94s. Its paper
  body was inaccessible (bot check); this is from the abstract and code.
- **Other refiners:** Xu & Gao (2025) refine ET-Flow / MCF / DMT samples in Cartesian space and use xTB only for
  properties. TD, CoarsenConf and Uni-Mol do not refine TD samples with a torsion constraint.
- **Not found in the works and searches checked:** torsion-exact refinement (P_W + restoration), seed-torsion
  targets, and a torsion-restrained physics control on the same samples. This is not proof of absence.

## 11. Novelty assessment
| element | new? | does the evidence reward it? |
|---|---|---|
| Learned refinement of TD samples | No (DiSCO) | — |
| Exactly torsion-preserving refinement (P_W + Restore_W) | Not found in checked works | No: the projector ablation fails (H2) |
| Seed-torsion targets (pure-L) | Not found in checked works | Indirectly: they make the head's gain purely local, which is what lets xTB be compared fairly |
| Torsion-restrained xTB control on the same samples | Not found in checked works | **Yes**: it changes the conclusion, from "the learned head improves TD by 37 %" to "a cheaper physical step improves it more" |

The genuine contribution is an evaluation finding, not a method:
- **Finding:** on GEOM-QM9, a learned local-geometry refiner looks strong against raw TD but is dominated by a
  torsion-restrained xTB step on the same samples. GEOM-QM9 references are xTB geometries, which explains the
  ranking.
- **Scope:** the finding is limited to QM9, one TD checkpoint and one refiner.
- **Projector:** it is a correct and well-tested constraint, but not a contribution supported by its ablation.

## 12. Limitations
- **Scope:** GEOM-QM9 only; one TD checkpoint; 184 paired molecules; one refiner of ours. DiSCO and other refiners
  were not run through the constrained-xTB control.
- **References:** they are xTB geometries, so the comparison favours xTB by construction. A DFT-referenced or
  DRUGS-scale evaluation could rank the methods differently. It was not available here.
- **Compute:** CPU only. Read-out GPU cost unmeasured; xTB force constant not varied.
- **Training:** 20 epochs, not converged. Extending it is not justified by the decision rule (`RESEARCH_DECISION.md`).
- **Strata:** exploratory, uncorrected for multiplicity.

## 13. Reproducibility
- **Code:** `jmcg_stage`, rsynced to `/scratch/nishanth.r/jmcg`. The code md5s are in `e1/code.md5`; 116 tests passed
  in the evaluation job.
- **Inputs:** `e1_manifest.sha256` (24 files). Checkpoints: `train/T2_*/checkpoints.sha256`.
- **Commands:** `slurm_jmcg/exp1_v2_t2_e1.sbatch` (evaluation), the training commands in `train_jobs.txt`, and
  `research/experiments/COMMANDS_M1-M7.md`.
- **Outputs:** predictions and RMSD matrices in `/scratch/nishanth.r/jmcg/results/exp1/T2_e1/`; JSON and logs in
  `e1/`. All bootstraps are seeded: `default_rng(0)` for the endpoints and `default_rng(1)` for addendum C.

## 14. Final research decision
**Option C: freeze FlexiTors V2 as a completed investigation and prioritise joint microstate–conformer generation.**

This is the outcome of the addendum-B rule:
- endpoint 1 is "yes";
- endpoint 4 is "worse";
- the head is not cheaper than xTB.

Option D does not qualify, because more training could not reach rule A while the head costs 4× more than xTB. The
comparison of options A / B / C on the brief's criteria is in `RESEARCH_DECISION.md`.
