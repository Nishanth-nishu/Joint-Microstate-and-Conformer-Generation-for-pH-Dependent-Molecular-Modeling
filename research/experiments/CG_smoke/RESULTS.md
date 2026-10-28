# Explicit-state conformer generation: smoke test decision report (2026-10-28)

All numbers below come from logs and CSVs under `run1/`, `run1_diag/` and `run2/`, and from `logs/jmcg_cg*_*.log`.
Times are given in UTC and IST. No performance claim is made. With 8 conformers per state and no reference ensembles,
these runs check plumbing, chemistry handling and validity only. The joint architecture was not started.

## 1. Job status and logs

| job | what | state | UTC | IST |
|---|---|---|---|---|
| 43186 | smoke 1 (first version) | cancelled before it started (queue: Priority); replaced by 43422 | – | – |
| 43158 | tautomer pilot, stages 4–6 | COMPLETED, exit 0 | 2026-10-27 17:28:51 → 21:06:26 | 22:58:51 → 02:36:26 (28 Oct) |
| 43422 | smoke 1 (`run1/`) | COMPLETED, exit 0 | 21:06:39 → 21:10:13 | 02:36:39 → 02:40:13 |
| 43714 | diagnostic: re-validation and cross-process drift (`run1_diag/`) | COMPLETED | 21:14:10 → 21:16:23 | 02:44:10 → 02:46:23 |
| 43717 | smoke 2 = smoke 1 rerun after fixes (`run2/`) | COMPLETED | 21:17:40 → 21:19:29 | 02:47:40 → 02:49:29 |
| 43719, 43721 | test-only job, my `--wrap` used `sh` and failed at line 4–5 | FAILED, did no work | – | – |
| 43724 | test-only job (`slurm_jmcg/cg_tests.sbatch`) | COMPLETED | 21:21:41 → 21:21:55 | 02:51:41 → 02:51:55 |

**Checkpoint used:** `/scratch/nishanth.r/tordiff/workdir/pretrained/qm9_1order/best_model.pt` (SHA-256 4199663e…).
It has 44 node features and types H/C/N/O/F. This is the **first-order** variant, not `pretrained/qm9_default` (second-order, TD's
default). The selector took the first `pretrained` QM9 directory it found. No GEOM-DRUGS checkpoint exists on disk.

## 2. Unit tests (job 43724, final code)

- `test_states`: **20 passed, 0 failed.** Covers accepted transformations, rejection of redox, charge moved without a
  proton, a carbanion, a different skeleton, a salt, an unparsable SMILES, a valence error, a radical, unassigned
  stereo and an inverted parent stereocentre, plus explicit stereo enumeration.
- `test_backend`: **10 passed, 0 failed.** Covers:
  - featuriser width = 44 for all 22 TD-compatible states
  - charge one-hot on the glycine zwitterion
  - S gives a hard incompatibility; ibuprofen is out of distribution (OOD) by size
  - state-dependent torsion counts: phenol 1, phenolate 0, 2-hydroxypyridine 1, 2-pyridone 0, pyridinium 0
  - same-process reproducibility: identical for the same seed; 3.07 Å apart for seed 0 vs seed 1
  - **new:** conformers within one call are distinct for seed 0 and seed 1
  - E/Z lock: 8/8 locked vs 7/8 unlocked correct
  - SDF round trip on 23 accepted states: atom order, charges, H, bond types, isomeric SMILES, SD tags
  - geometry controls on ibuprofen and acac Z-enol: NaN, +0.5 Å stretch, +1.5 Å break, clash, mirror image and
    E/Z flip all detected, unperturbed conformers pass
  - **new:** clash screen accepts an intramolecular O–H···O=C H-bond at 1.64 Å and 1.70 Å and rejects one at 1.20 Å

Earlier runs of the same suites: in job 43422, 20/20 and 8/8 (before the new tests). In job 43714, 8/9 and in job 43717, 9/10. Both of those failures were
the new H-bond test. Each time the fault was in the **test code**, not the code under test:
- the O–H atom lookup used `GetTotalNumHs()` without `includeNeighbors`;
- the ETKDG geometry had the methyl, not the C=O, facing the OH.

## 3. Per-state results: smoke 2 (`run2/qm9/summary.csv`, seed 0, 8 conformers requested per method)

CSV: 24 rows = 22 states of 10 distinct parents (10 of the 22 are the parent itself) + 2 negative controls (12 parents
in total). Categories, from the code:
- **A** = neutral, inside vocabulary and size: 12 states.
- **B** = charged, inside vocabulary and size: 8 states.
- **C** = more than 9 heavy atoms: 2 states.
- **D** = unsupported or invalid: 2 controls.

Every accepted state prepared as specified, and the SDF round trip passed for all 23 accepted states. For every
generated conformer:
- connectivity rebuilt from 3D matched the state;
- bonds and angles were inside the bounds;
- no connectivity check was skipped.

The failures in the table below come only from clash or stereo.

Column notes: *valid* counts conformers that pass every check. *cl/uq* is the number of heavy-atom Butina clusters at
0.5 Å and of unique conformers at 0.1 Å. *tors* is TD's torsion count for that state.

| state | cat | class | q | heavy | tors | ETKDG valid | TD valid | TD+lock valid | ETKDG cl/uq | TD cl/uq | TD s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| propylamine_neutral | A | parent | 0 | 4 | 3 | 8/8 | 8/8 | – | 1/2 | 1/3 | 3.49* |
| propylamine_cation | B | protonation | +1 | 4 | 3 | 8/8 | 8/8 | – | 2/3 | 3/3 | 0.78 |
| glycine_neutral | A | parent | 0 | 5 | 3 | 8/8 | 8/8 | – | 1/1 | 1/2 | 0.70 |
| glycine_zwitterion | B | tautomer_zwitterionic | 0 | 5 | 2 | 8/8 | 8/8 | – | 1/1 | 1/3 | 0.47 |
| gaba_neutral | A | parent | 0 | 7 | 5 | 8/8 | 8/8 | – | 6/8 | 6/8 | 0.94 |
| gaba_zwitterion | B | tautomer_zwitterionic | 0 | 7 | 4 | 8/8 | **7/8** | – | 7/8 | 5/7 | 0.88 |
| hydroxypyridine | A | parent | 0 | 7 | 1 | 8/8 | 8/8 | – | 1/1 | 1/1 | 0.55 |
| pyridone | A | tautomer | 0 | 7 | 0 | 8/8 | 8/8 | – | 1/1 | 1/1 | 0.01 |
| acac_diketo | A | parent | 0 | 7 | 4 | 8/8 | 8/8 | – | 4/6 | 4/8 | 0.77 |
| acac_Z_enol | A | tautomer | 0 | 7 | 5 (4 locked) | 8/8 | **6/8** | **7/8** | 2/2 | 4/6 | 0.76 |
| histamine_Ntau | A | parent | 0 | 8 | 3 | 8/8 | 8/8 | – | 6/8 | 5/8 | 0.88 |
| histamine_Npi | A | tautomer | 0 | 8 | 3 | 8/8 | 8/8 | – | 4/6 | 4/6 | 0.91 |
| histamine_cation | B | protonation | +1 | 8 | 3 | 8/8 | 8/8 | – | 3/4 | 5/8 | 0.97 |
| histamine_dication | B | protonation | +2 | 8 | 3 | 8/8 | 8/8 | – | 2/5 | 5/7 | 1.07 |
| ibuprofen_S_neutral | C | parent | 0 | 15 | 8 | 8/8 | **6/8** | – | 3/6 | 7/8 | 2.58 |
| ibuprofen_S_anion | C | protonation | −1 | 15 | 7 | 8/8 | 8/8 | – | 5/6 | 6/8 | 2.41 |
| pyridine | A | parent | 0 | 6 | 0 | 8/8 | 8/8 | – | 1/1 | 1/1 | 0.01 |
| pyridinium | B | protonation | +1 | 6 | 0 | 8/8 | 8/8 | – | 1/1 | 1/1 | 0.01 |
| phenol | A | parent | 0 | 7 | 1 | 8/8 | 8/8 | – | 1/1 | 1/1 | 0.60 |
| phenolate | B | protonation | −1 | 7 | 0 | 8/8 | 8/8 | – | 1/1 | 1/1 | 0.01 |
| acetamidine | A | parent | 0 | 4 | 3 | 8/8 | 8/8 | – | 1/1 | 1/1 | 0.46 |
| acetamidinium | B | protonation | +1 | 4 | 3 | 8/8 | 8/8 | – | 1/1 | 1/1 | 0.52 |
| NEG_cysteine_S | D | (accepted state) | 0 | 7 | – | 8/8 | not run: S not in qm9 vocabulary | – | 6/6 | – | – |
| NEG_wrong_skeleton | D | – | – | – | – | rejected: heavy-atom skeleton differs | – | – | – | – | – |

\* The first TD call includes model warm-up. ETKDG took 0.01–0.13 s per state for 8 conformers. TD ran on 4 CPU
threads.

**Totals by category (not merged).** All counts are valid / generated conformers.

| category | states | ETKDG | TD | TD + E/Z lock |
|---|---|---|---|---|
| A | 12 | 96/96 | 94/96 | 7/8 (acac Z-enol only) |
| B | 8 | 64/64 | 63/64 | – |
| C | 2 | 16/16 | 14/16 | – |
| D | 2 | 8/8 (cysteine only) | none, by design | – |

### Every failure, explained

All failures are TD conformers; ETKDG had none.

| state, method, conformer | failed check | cause |
|---|---|---|
| acac_Z_enol, td, 0 | stereo, round-trip SMILES, clash | The C=C flipped Z→E: TD treats a non-ring C=C as a torsion. The flip also left a carbonyl O···methyl H contact at 1.63 Å. The lock arm does not show this failure. |
| acac_Z_enol, td, 2 | clash | Enol O···methyl C–H at 1.85 Å. The C3–C4 rotation placed the methyl, not the C=O, next to the enol O. |
| acac_Z_enol, td_lockdb, 6 | clash | Same motif: O···H at 1.903 Å against a 1.904 Å threshold. **Borderline.** |
| gaba_zwitterion, td, 4 | clash | Carboxylate C···ammonium H at 1.99 Å against a 2.03 Å threshold: a collapsed gas-phase salt bridge. **Borderline.** |
| ibuprofen_S_neutral, td, 0 | clash | H···H at 1.35 Å against a 1.44 Å threshold. |
| ibuprofen_S_neutral, td, 7 | clash | H···H at 1.39 Å against a 1.44 Å threshold. |

### Reproducibility (run 2)

- **ETKDG:** bit-identical across separate processes for the same seed.
- **TD within one process:** bit-identical for the same seed.
- **TD across processes, same seed:** max atom deviation of 1.5e-3 Å (glycine zwitterion), 3.1e-3 Å (acac Z-enol) and 1.1e-2 Å
  (histamine dication). Validity and cluster counts did not change.
  - The diagnostic job 43714 ran this with 4 or 1 threads and with a random or fixed `PYTHONHASHSEED`. The drift
    stayed (4e-3 to 1.1e-1 Å with the old seeding), so neither threads nor hash randomisation is the cause.
  - My working explanation is floating-point summation order in CPU kernels. **That is not verified.**
- **Seed 0 vs seed 1:** conformers differ by 3.1–7.8 Å.

### Comparison with ETKDG (descriptive only)

Both methods produced valid geometry for every A and B state, except the TD failures listed above. ETKDG had no
failures at all.

TD's bond lengths and angles are copied from its ETKDG seed, so for TD the bond/angle check only tests the seed, not
the model.

TD does not apply the model at all for states with 0 torsions: pyridone, pyridine, pyridinium and phenolate. For
those, the TD output is the ETKDG seed geometry.

Hydroxypyridine, phenol and acetamidine have torsions, but those torsions move only H atoms: 1 unique heavy-atom
conformer.

For the 14 remaining flexible states, counting unique heavy-atom conformers at 0.1 Å out of 8:
- TD has more in 9 states, the same number in 4 and fewer in 1 (gaba zwitterion).
- With 8 conformers and no reference ensemble this is descriptive only. More unique conformers does not mean better
  conformers.

MMFF94s single-point energies (unrelaxed) are in the per-state JSONs. They are only comparable within one state.
For the charged and zwitterionic states they are gas-phase values with no physical meaning: for example, 116–120
kcal/mol for the histamine dication. They are not used for any conclusion.

## 4. Bugs found and fixes

1. **Identical conformers with seed 0.** This was the most serious bug.
   - With `randomSeed=0`, `EmbedMultipleConfs` returned 8 identical conformers. In run 1 every ETKDG baseline had
     1 cluster and median RMSD 0.00. TD seed geometries were identical too, so TD diversity came from torsions
     alone.
   - This is the same RDKit pitfall already recorded for the N1 pilot (2026-10-26), and I repeated it.
   - **Fix:** each conformer is embedded separately, with its own seed from `numpy.default_rng(seed)`
     (`td_backend.conformer_seeds` / `embed_n`). ETKDG and the TD seeds use the same per-conformer seeds.
   - New regression test added. **Run 1's ETKDG comparison and diversity numbers are invalid.**
2. **H-bonds counted as clashes.**
   - A polar H···N/O hydrogen bond (about 1.6–2.0 Å) falls below 0.7 × the vdW sum, which is 1.90 Å for H···O. So the clash
     rule rejected correct H-bonded geometries.
   - Re-validating the run 1 SDFs changed 3 verdicts; the seed-1 repro run changed 2 more. Every changed verdict
     was a polar-H···O contact.
   - **Fix:** a polar H···N/O pair is now flagged only below 1.5 Å. Such contacts are counted (`n_hbond_contacts`),
     and the clashing pairs are now listed per conformer. New test added.
3. **Ambiguous checkpoint selection.** The run picked `qm9_1order`, not `qm9_default`. **Fixed in smoke 3**
   (`slurm_jmcg/cg_smoke3.sbatch`): the checkpoint is pinned by path and the job stops if its SHA-256 does not start
   with 6a9a6ee9b0e51bda. See section 8.
4. **Cross-process TD drift.** Measured above and not fixed. The rule from now on: reproducibility claims are
   "same verdicts and clusters, coordinates within 0.11 Å", not bit-identical.
5. **Fixed in the earlier code audit** (before any output; listed in `10_decision_log.md`):
   - silent stereoisomer choice
   - redox and resonance forms accepted
   - zwitterion mislabelled
   - no parent stereo check
   - weak SDF check
   - weak bond check
   - imine E/Z false failures
   - unchecked connectivity counted as valid
   - lock arm run with no E/Z element

## 5. Which outputs can be interpreted

- **Pipeline correctness: yes, for all 23 accepted states.** State preparation, rejection of wrong transformations,
  SDF fidelity, stereo handling, compatibility gating and failure reporting all behaved as specified.
- **TD on A states (neutral, at most 9 heavy atoms): only as a sanity check.** The model runs, and its output is valid
  geometry at a rate similar to ETKDG. There is no reference ensemble, so nothing about accuracy or coverage can be
  said.
- **TD on B states (charged): not interpretable as learned behaviour.**
  - QM9 train has 111 of 106,041 molecules with any formal charge (0.10%), and **all 111 are nitro groups**: N⁺ and O⁻
    with no H. **None is net-charged.**
  - So ammonium, carboxylate, pyridinium, amidinium, imidazolium and phenolate environments never appeared in
    training.
  - The featuriser can represent these states, and the outputs are valid geometries, but the torsion distribution is
    an extrapolation.
- **TD on C states (ibuprofen, 15 heavy atoms, QM9 maximum 9): not interpretable.**
- **The E/Z lock arm** removed the stereo flip: 8/8 vs 7/8 correct E/Z in the unit test, and no stereo failure in
  run 2. It did not prevent the C3–C4 orientation clash.

## 6. States that need a different checkpoint or retraining

- **Different checkpoint:**
  - **C (size):** both ibuprofen states. Needs a GEOM-DRUGS checkpoint, and none is on disk.
  - **D:** cysteine (S is not in the qm9 vocabulary). Also needs a DRUGS checkpoint.
- **Retraining or a different training set:**
  - **B (charge):** all 8 states. They need a model trained on data with protonated or deprotonated environments.
  - Whether GEOM-DRUGS contains enough of them is unknown; it should be inventoried in the same way as QM9 before
    anyone assumes so.
- **Model changes, not just a checkpoint:** states where E/Z must be kept fixed. TD's torsion set includes non-ring
  stereo double bonds. The lock arm handles this at sampling time; training would need the same mask.

## 7. Smallest defensible next experiment (not started)

**Wrapper-fidelity check.**
- Run `jmcg.condgen` against TD's own `generate_confs.py` on about 50 neutral GEOM-QM9 test molecules.
- Use the same pinned checkpoint, `qm9_default`, and the same number of conformers.
- Report COV-R and MAT-R for both against the GEOM reference ensembles.

This checks that state handling, the seeding fix and the E/Z lock do not degrade TD on its own domain. It runs on CPU
in minutes. Nothing about charged states can be claimed until this passes. Separately, and without training, a
read-only count of charged environments in GEOM-DRUGS would show whether a B-category model can be trained from
existing data.

**Status:** done. The wrapper-fidelity check passed (FAITHFUL; `../CG_fidelity/RESULTS.md`).

## 8. Smoke 3: `qm9_default` pinned (job 45117, 2026-10-28 15:56–15:59 UTC / 21:26–21:29 IST)

Identical to smoke 2 except for the checkpoint: `pretrained/qm9_default/best_model.pt`, SHA-256 6a9a6ee9b0e51bda…,
pinned by path and checked by hash. Same 24 states, seed 0, 8 conformers, same code. Outputs: `run3/`.

- **Tests:** test_states 20/20, test_backend 10/10.
- **Valid/generated, by category (smoke 2 → smoke 3):**

| category | ETKDG | TD | TD + E/Z lock |
|---|---|---|---|
| A (12 neutral in-distribution states) | 96/96 → 96/96 | 94/96 → 94/96 | 7/8 → 7/8 |
| B (8 charged states) | 64/64 → 64/64 | 63/64 → 63/64 | – |
| C (2 ibuprofen states, OOD size) | 16/16 → 16/16 | 14/16 → 16/16 | – |
| D (2 controls) | cysteine ETKDG 8/8, TD hard-rejected (S); wrong skeleton rejected: unchanged | | |

- **The four smoke-3 failures:**

| state (category) | method | conf | reason |
|---|---|---|---|
| acac_Z_enol (A) | td | 0 | Z→E flip of the enol C=C **and** O···H 1.78 Å (O to a non-polar H, so not H-bond exempt) |
| acac_Z_enol (A) | td | 4 | Z→E flip only (clash-free, bonds and angles OK) |
| acac_Z_enol (A) | td_lockdb | 2 | O···H 1.77 Å (methyl H to the carbonyl O; limit 1.90 Å); E/Z kept |
| histamine_dication (B) | td | 5 | H···H 1.35 Å (limit 1.44 Å) |

  - Every bond, angle, connectivity, SDF and stereocentre check passed for every conformer.
  - The acac Z-enol flips again on the same seed embedding (conformer 0), as it did with `qm9_1order`. The E/Z
    problem is not checkpoint-specific. The lock arm removes the flips but not the methyl orientation clash.
- **Per-state changes vs smoke 2 (TD only; ETKDG is bit-identical because its code and seeds are unchanged):**
  - gaba_zwitterion 7/8 → 8/8.
  - ibuprofen_S_neutral 6/8 → 8/8.
  - histamine_dication 8/8 → 7/8.
  - With 8 conformers per state these are single-conformer changes, not evidence that either checkpoint is better.
- **Clusters (TD vs ETKDG, 12 states with more than one cluster in either):** TD more in 6, the same in 3, fewer in 3.
  Smoke 2 was 9 more / 4 same / 1 fewer over 14 states. Descriptive only, at n = 8.
- **Reproducibility (seed 0, two processes):**
  - ETKDG is bit-identical.
  - TD max |dx|: glycine zwitterion 7e-4 Å, histamine dication 6.8e-2 Å, **acac Z-enol 0.49 Å** (smoke 2: at most
    1.1e-2 Å).
  - Verdicts were identical across the two processes, with the same failing conformers and the same valid counts.
  - The 0.11 Å rule in item 4 of section 4 is therefore **broken** for this checkpoint. Small floating-point
    differences between processes are amplified by the SDE; the root cause is still unverified. Reproducibility claims
    must be stated as "same verdicts", not as a coordinate tolerance, until the drift is explained.

**Interpretation is unchanged** (sections 5–6). The TD run on A states is now on the checkpoint that passed the
fidelity check, so it can be cited as a sanity check. B and C remain uninterpretable. Smoke 1 and 2 TD outputs
(`qm9_1order`) are superseded by smoke 3.
