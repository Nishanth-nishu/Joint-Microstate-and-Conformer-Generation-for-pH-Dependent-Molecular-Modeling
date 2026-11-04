# P1 (EXPLORATORY, not pre-registered): unminimised conformer features for the M5 extraction blocker

**Question.** Does a feature path without MMFF minimisation remove the M5 blocker? The blocker: the validator
rejects MMFF94s minima of zwitterions, carboxylates and phenolates, so those states have no φ.

**Constraints observed:**
- the MMFF path (E2 definition) is kept unchanged as the control;
- validator thresholds unchanged;
- the M3, E2 and M5 reports are not modified;
- nothing here re-tests the E2 pre-registered claim.

**Code:**
- `jmcg/joint/features.py`: new `path='unmin'` and `features_from_conformers`; the default is still `'mmff'`;
- `jmcg/joint/p1_unmin.py`: experiment driver;
- `jmcg/joint/sample.py`: `--feature_path`, default `mmff`;
- `jmcg/joint/test_joint.py`: 4 tests, all pass.

**Jobs:**
- 46561 (2026-10-29 10:27–11:55 UTC / 15:57–17:25 IST): panel, both RDKit paths, TD generation. Two failures:
  - one of 16 parallel TD processes was OOM-killed;
  - `compare` raised KeyError `psi1`, because the extract step did not write the 2D ψ columns.
- 46662: the resume failed on a `grep` exit code in its own diagnostics.
- 46673 (11:34–12:02 UTC / 17:04–17:32 IST): re-ran the killed TD chunk alone (peak RSS 3.5 GB), then
  td_features, compare and the sampler.
- Sbatch scripts: `slurm_jmcg/p1_unmin.sbatch`, `p1_unmin_resume.sbatch`.

## Arms (same panel, same 30 ETKDGv3 per-conformer seeds, same validator)
| arm | conformers | weights |
|---|---|---|
| **mmff** (control) | ETKDGv3 → MMFF94s minimisation | Boltzmann on MMFF94s, 298 K |
| **unmin** | ETKDGv3, no minimisation | uniform (φ1 = −ln n_valid carries no information on this arm) |
| **td** | TD `qm9_default` (sha `6a9a6ee9…`), E/Z lock default, n = 30, seed 0, re-validated in the jmcg venv | uniform |

**Panel** (`panel.csv`, 1,080 states):
- the 402 Tautobase E2 states;
- the 16 M3 smoke states;
- up to 150 states per category sampled (seed 0) from the states Uni-pKa scored in E1.

Categories come from formal charges:

| category | n |
|---|---|
| neutral | 491 |
| zwitterion | 161 |
| cation | 156 |
| other anion | 151 |
| carboxylate | 84 |
| phenolate | 37 |

**Caveat on "zwitterion":** the category means "both charge signs present". It includes 40 nitro/N-oxide states, so
the carboxylate-containing subset (66) is reported separately.

## 1. Control check
On the 402 E2 states, the mmff arm reproduces E2 `features_v2` exactly: n_valid equal 402/402, max |Δφ| = 0.0.

## 2. Feature-extraction success (state has ≥ 1 valid conformer), Wilson 95 % CI
| category | n | mmff | unmin | td (all states) | td where TD ran |
|---|---|---|---|---|---|
| neutral | 491 | 0.84 [0.80, 0.87] | **0.97** [0.95, 0.98] | 0.75 | 368/368 |
| zwitterion (any) | 161 | 0.17 [0.12, 0.23] | **0.77** [0.70, 0.83] | 0.45 | 73/74 |
| zwitterion with carboxylate | 66 | **0.03** [0.01, 0.10] | **0.73** [0.61, 0.82] | 0.14 | 9 of the 66 have TD output |
| carboxylate | 84 | 0.17 [0.10, 0.26] | **0.96** [0.90, 0.99] | 0.45 | 38/39 |
| phenolate | 37 | 0.43 [0.29, 0.59] | **0.95** [0.82, 0.99] | 0.54 | 20/20 |
| cation | 156 | 0.59 [0.51, 0.66] | 0.84 [0.77, 0.89] | 0.31 | 49/50 |
| other anion | 151 | 0.61 [0.53, 0.68] | 0.94 [0.89, 0.97] | 0.44 | 67/70 |

**Paired:** no state had φ under mmff but not under unmin. Unminimised gains, by category:
- carboxylate 67;
- zwitterion 97;
- phenolate 19;
- neutral 63;
- cation 39;
- other anion 50.

**Remaining unmin failures.** Mostly `rt_isomeric,stereo`: states with unassigned stereocentres, whose ETKDG
embedding picks a configuration the validator compares against a stereo-free target. This fails on every arm.

**TD's low state-level rates are coverage, not validity.** TD did not run on 459 of 1,080 states:
- atoms outside the QM9 vocabulary (S, Cl, Br, …);
- the state builder's stereo, skeleton and radical checks.

Among conformers TD did produce, most rejections are `clash`.

## 3. Why the MMFF path fails (diagnostic, first minimised conformer of each failing state)
The E2 RESULTS hypothesis (MMFF anion geometry against UFF-based bounds) is now measured. For example, MMFF's
C–O⁻ length falls outside the ±0.10 Å window around the distance-geometry bound.

| category | failing states | worst 1–2 pair after MMFF (top) | median worst 1–2 dev after MMFF (Å) | median worst 1–3 dev: after MMFF / same conformer unminimised (Å) |
|---|---|---|---|---|
| carboxylate | 70 | C–O⁻ 62/70 | 0.111 | 0.097 / 0.017 |
| phenolate | 21 | C–O⁻ 14/21 | 0.104 | 0.131 / 0.014 |
| zwitterion | 134 | C–O⁻ 52, N⁺–O⁻ 39 (nitro) | 0.122 | 0.134 / 0.021 |
| cation | 64 | C–N⁺ / N⁺–N | 0.144 | 0.189 / 0.026 |
| neutral | 80 | C–C 43 | 0.103 | 0.139 / 0.013 |

Per-state detail is in `compare/mmff_failure_diagnostics.csv`. The failure is a mismatch between MMFF94s minima and
the validator's distance-geometry bounds. It is not a sign that the structures are chemically wrong: the same
conformers before minimisation sit well inside the bounds.

## 4. Feature distributions on matched states (φ under both arms)
| | n | φ3 (polar SASA/100 Å²) Pearson | φ3 mean diff (arm − mmff) | φ2 (IMHB) Pearson | φ2 mean: mmff → arm |
|---|---|---|---|---|---|
| unmin vs mmff, all | 652 | 0.994 | +0.023 | 0.62 | 0.32 → 0.09 |
| unmin vs mmff, charged | 241 | 0.988 | +0.040 | 0.55 | 0.70 → 0.19 |
| td vs mmff, all | 459 | 0.995 | +0.010 | 0.42 | 0.24 → 0.07 |

- **φ3, the descriptor that carried the E2 result, is nearly identical across paths.** Unminimised conformers
  expose slightly more polar surface.
- **φ2 is not:** MMFF minimisation creates intramolecular H-bonds (salt bridges in charged states) that the
  unminimised ensembles rarely contain. φ2 is therefore path-dependent and should not be compared across paths.

## 5. E2 predictions (exploratory; same pre-registered fitting procedure, out of fold, scaffold folds, seed 0)
**128 measured pairs eligible under all three paths:**

| arm | MAE |
|---|---|
| S0 | 3.55 |
| S1 | 2.73 |
| C1 | 2.51 |
| J1-mmff | 2.31 |
| J1-unmin | 2.25 |
| J1-td | 2.18 |

- J1-mmff − J1-unmin = +0.06 [−0.08, 0.21]; J1-mmff − J1-td = +0.12 [−0.02, 0.28]. No path is distinguishable.
- Out-of-fold predictions agree (Pearson 0.97, mean |Δ| 0.6 log units).

**All pairs eligible under unmin (206 of 221; 150 under mmff):**

| arm | MAE |
|---|---|
| S1 | 2.68 |
| C1 | 2.64 |
| J1-unmin | 2.66 |
| S1 + φ3 only | 2.47 |

- Δ_S = 0.03 [−0.25, 0.27] and Δ_C = −0.01 [−0.34, 0.31]: on this larger set the full J1 shows **no** gain.
- On the same set, the one-term φ3 model is better than S1 by 0.21 (no CI computed for it here).
- The J1 refit gives φ2 a negative coefficient (−3.3), against +0.66 under mmff. This is consistent with φ2 being
  path-dependent.

**All pairs eligible under td (167):**

| arm | MAE |
|---|---|
| S1 | 2.79 |
| C1 | 2.65 |
| J1-td | 2.50 |

Δ_S = 0.30 [0.09, 0.50]; Δ_C = 0.16 [−0.05, 0.41].

These are post hoc fits on the same Tautobase data E2 used. They neither confirm nor overturn the pre-registered E2
verdict. The 56 pairs newly eligible under unmin are the ones the E2 selection caution warned about, and adding them
removes J1's margin over S1.

## 6. M5 sampler with the unmin path (smoke parents, pH 7.4, coefficients `compare/coef_unmin.json`)
- `fallback_levels_theta0` is empty and `n_pairs_without_conformer` = 0 for S1 and J1. Before (MMFF path, job
  45584) there were 6 fallback levels and 40 of 120 pairs without a conformer.
- **Every level now uses the coupled score.** Glycine's zwitterion (p = 0.996) and ibuprofen's carboxylate
  (p = 0.999) get conformers.
- J1 and S1 state probabilities differ by ≤ 3e-4 on these parents. The coupling still changes nothing visible here.
- The unmin coefficients come from the exploratory refit above, fitted on neutral tautomer pairs only. They are
  applied unchanged to charged levels: a modelling assumption, not a validated one.

## Conclusion (what the measurements support)
1. **The blocker is the MMFF-minima-vs-bounds mismatch, and the unminimised path removes most of it.** State-level
   extraction rises to:
   - 0.96 for carboxylates (0.17 under mmff);
   - 0.95 for phenolates (0.43);
   - 0.73 for carboxylate zwitterions (0.03).
   In the M5 demo, no level falls back to θ = 0. Validator thresholds were not changed.
2. **It does not solve everything.**
   - 23 % of zwitterions and 16 % of cations still fail, mostly on unassigned-stereo states (a separate,
     path-independent issue).
   - The unmin path changes the meaning of φ2 (and φ1 becomes constant per n_valid). So J1 coefficients are not
     transferable between paths.
   - On the enlarged unmin-eligible set, the full J1 model shows no gain over S1.
3. **TD is not a practical generator for this panel.** It covers only 621 of 1,080 states (QM9 vocabulary and
   state-builder limits), and its rejections are mostly clashes. Where it runs, its φ3 matches mmff (r = 0.995).
4. **Open for decision:**
   - whether to adopt the unmin path (with φ = {φ3} or a re-defined φ1/φ2) as the joint model's feature definition;
     that change needs a new pre-registration before any confirmatory claim;
   - whether to keep φ2 at all.
