# Milestone 2: canonical data layer v1 (2026-10-28)

Built by `jmcg/data/build_canonical.py` (job 45295, 17:39 UTC / 23:09 IST) into `data_manifests/canonical_v1/`. Audited
by `jmcg/data/audit.py` (`audit/audit.json`, `audit/audit.md`). Identity policy tests: `jmcg/data/test_canonical.py`,
7/7 passing. Command: `sbatch slurm_jmcg/m2_data.sbatch`. RDKit 2026.09.1, Python 3.10.12.

## Tables
- **`states.csv`**: one row per microstate.
  - Key: `state_id = 'S' + sha1(canonical isomeric SMILES)[:16]`.
  - Columns: formal charge, heavy-atom count, elements, number of unspecified stereo elements, parent, and origin
    datasets.
- **`parents.csv`**: one row per parent.
  - Key: `parent_id = 'P' + InChIKey` of the parent form (uncharged, then canonical tautomer, stereo kept).
  - `family_key` is the first InChIKey block (connectivity only). Also: Murcko scaffold, split group, merged split
    component, and the method used to build the parent form.
- **`observations.csv`**: one row per experimental label.
  - Columns: task, source and source record, role, split, value, units, sign convention, temperature, solvent, label
    status, the acid/base or tautomer state ids, the reported state, the linked parents, and quality flags.
- **`dropped.csv`**: every record not loaded, with the reason.
- **`manifest.json`**: SHA-256 of inputs and outputs, split policy, counts.

## Identity policy
Each choice below is tested in `test_canonical.py`.
- **Distinct microstate ids:** distinct protonation states, tautomers, zwitterions and stereoisomers always get
  distinct `state_id`s. Unspecified stereo stays unspecified.
- **Rejected, never stripped:** salts and mixtures, radicals, and SMILES that don't parse are rejected with a reason.
- **Protomers and tautomers** share a parent. The shared parent comes from RDKit's canonical tautomer of the uncharged
  form.
  - RDKit's default canonicaliser **removes stereo at tautomer-active atoms**, so (R)- and (S)-alanine collapsed into
    one parent (found by the test, job 45290). The canonicaliser now keeps sp3 and double-bond stereo.
  - 40 parents fall back to the microstate itself, because their canonical-tautomer or uncharged form does not
    round-trip through SMILES. These are recorded in `parent_method`.
- **Stereoisomers** are distinct parents but share a `family_key`. Splits group by scaffold, or by family when acyclic,
  so stereoisomers never straddle splits.
- **Split assignment** is a deterministic hash of the split component (`sha256('jmcg-split-v1|component') % 100`:
  under 70 train, under 85 val, else test).
  - Components are unions of split groups linked by any observation, for example a tautomer pair whose tautomers
    standardise to different families (16 measured-or-other pairs).
  - Curated observations whose linked families occur in an external test set are `excluded_external_overlap`.

## Counts (`manifest.json`)
- 15,502 observations, 29,280 microstates, 11,355 parents.
- 3 records dropped: Dwar-iBonD transitions whose acid and base sets are not one proton apart.
- 137 split groups merged through links.

| task / role | train | val | test | other |
|---|---|---|---|---|
| macro-pKa, IUPAC curated (acid/base microstates unknown) | 4,758 | 666 | 968 | 3 excluded (external overlap) |
| macro-pKa, external test (SAMPL6/7/8, Novartis) | – | – | – | 356 external_test |
| macro-pKa, Dwar-iBonD (released Uni-pKa's training set) | – | – | – | 8,229 unipka_train |
| tautomer log K, all (221 measured, 166 estimated, 135 calculated) | 331 | 79 | 111 | 1 excluded |
| tautomer log K, measured only | 171 | 21 | 28 | 1 excluded |

## Audit (`audit/audit.json`)
- **No conformer references.** None of these datasets has a 3D reference ensemble for these states. The only
  conformer references available anywhere in the project are GEOM-QM9 test molecules, which are neutral.
- **Chemistry:**
  - State charges range from −5 to +4: 13,233 neutral, 8,482 at +1 and 5,317 at −1.
  - Elements: C, H, N, O, S, Cl, F, Br, P, I, Se, Si, As.
  - Heavy atoms: median 13, 95th percentile 26, maximum 59.
  - 5,597 states have unspecified stereo.
  - Only 4,070 states are inside TD-qm9's vocabulary and size limit (H, C, N, O, F and at most 9 heavy atoms).
- **Missing temperature:** Uni-pKa sets don't record it, and Tautobase leaves it blank for 516 of 522 records. It is
  stored as "not recorded", never imputed.
- **Duplicates and conflicts:**
  - 85 IUPAC family × pKa-index duplicates come from different input forms of one compound; 6 of them differ by more
    than 1 unit.
  - 24 IUPAC entries have a spread greater than 1 among their own repeat measurements.
  - 4 Tautobase pairs are duplicated, with no conflicting values and no sign conflicts.
- **Leakage:**
  - L1, L1b and L2 = 0: no family, linked family or scaffold component crosses train/val/test.
  - L3: 18 val and 18 test parents have a training parent at Morgan-2 Tanimoto ≥ 0.9, mostly homologous series. These
    are reported, not removed.
  - L4: 267/857 test and 145/529 val families are in Dwar-iBonD. This matters only when evaluating the released
    Uni-pKa on these splits.
  - L5: the old `sensitive_set_v2` split had 26 families in more than one split (protonation-specific keys). The new
    split agrees with it for 3,077/6,916 records, so v2 results are not comparable split-for-split.
- **Task-specific splits:** the measured tautomer set is small (28 test pairs). Tautomer experiments keep the T pilot's
  5-fold scaffold-grouped CV over the measured pairs, which uses no validation-based selection. The canonical split
  above is used for pKa-task fitting and selection.
