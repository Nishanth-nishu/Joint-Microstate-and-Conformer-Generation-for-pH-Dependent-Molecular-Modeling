# Code review P2 (validity focus): Stage A milestone 2, commits 0f70dc5..53029b6

Reviewer: planner P2, 2026-10-17. Checked against:
- EXP-3 §9, EXP-5 §12, EVAL §1.1/§3.1/§7;
- `round5/DECISION.md`;
- the `## Answer` sections in `round5/consult/`;
- `round5/IMPLEMENTATION.md` (milestone 2).

No code was edited and nothing was committed.

**Local re-run** (Windows, Python 3.11, RDKit 2026.03.6):

| Suite | Result | Note |
|---|---|---|
| `flexitors/tests` | 89 passed, 13 skipped | matches IMPLEMENTATION §2 |
| `tools/tests/test_lsa_floor.py`, `test_leakfree.py`, `test_leakfree_budget.py` | 15 passed | |
| `tools/tests/test_round2_model.py` | does not collect locally (`No module named 'utils'`) | pre-existing round-2 test that needs the TD venv; not a milestone-2 regression |

**Milestone-1 findings are fixed** (`a5b5d32`, `76785ae`, `54e3173`, `53290e9`, `b6e4aa8`, `02e39f6`, `d68ba5d`):
- B1 (validity vs strained cages), B2 (Gate-0 threshold recipe), S1–S8, N1, N3, N4, N8 are all fixed.
- The extra-bond rule is now relative to the reference.

**Good in milestone 2:**
- Per-conformer streams keyed by SMILES.
- The reverse loop is bit-identical to TD's `sample()` on CPU.
- K1 is excluded from the seed randomisation **and** masked from updates (DECISION fix 3).
- `UPDATE_SIGN` is measured against TD.
- The late planarity fold follows §9.2.
- F_filter fills the output to 2K and flags the fills.
- PG goes through TD's own `sample()` at batch 128.
- Analyses were committed before any result.

## BLOCKER

**B1. Failure sets depend on how many conformers a run draws, so contrasts compare different universes**
(`flexitors/sampling/generate.py:83-98`, `:138-140`; used by `experiments/common/stageA_tasks.tsv` and both analyses).

- **Mechanism.** `embed()` runs one `EmbedMolecule` per conformer index and returns `None` for the whole molecule if
  **any** index in the requested range fails. A run that draws more conformers is therefore more likely to lose the
  molecule:
  - the seed-0 TD pool draws max(128, 4K);
  - Exp 5 seeds 1–2 and PG4K draw 4K;
  - PG2K and the Exp 3 arms draw 2K.
  Stochastic ETKDG failures exist: some seeds of mildly strained molecules fail while others succeed.
- **Failure scenarios:**
  - a molecule whose conformer 150 fails is a failure in seed 0's P0, P1 and P3 (all read the TD pool), and in Exp 3's
    `O_s0` (the pool prefix, `stageA_tasks.tsv` row `exp3/qm9/O_s0`). The 2K arms (K1, K2_inf, …) and PG2K/PG4K
    succeed on it;
  - under the TD convention (`experiments/exp5/analysis.py:98-99`, COV = 0), P3 − P2b and every Exp 3 arm − O then
    gain coverage from an embedding accident, not from the method;
  - the Exp 3 analysis instead silently drops the unpaired molecule (`experiments/exp3/analysis.py:76-77`, `:85-86`),
    so its failure asymmetry is invisible.
- **Fix:** make a molecule's failure in a run that reads a prefix [0, n) depend only on indices < n:
  - record per-index embedding failures;
  - keep the pool up to the first failing index;
  - an arm needing n columns fails iff a failed index is < n.
  That is exactly what a standalone n-conformer run does. Then assert in both analyses that the two arms of every
  contrast have identical failure sets on their universe, and print failures per arm and seed (EVAL §3).

## SHOULD-FIX

**S1. Splice on the molecules where the arm's operator acts, not on the whole double-bond subset**
(`flexitors/sampling/splice.py:16-24`; `stageA_tasks.tsv`: K1, K1_wrong and F_filter use `keys_double.txt`).

- K1 and K1_wrong freeze **specified** bonds only (`stereo_guard.py:99-100`), and F_filter acts only on molecules with
  a target (`generate.py:180-181`). On QM9 these operators are the identity on 99 of the 102 molecules, which the arm
  nevertheless regenerates.
- On CPU those 99 molecules equal O bit for bit. On GPU they are independent re-draws (0.265 Å deviation), so the
  K1 − O contrast is 3 molecules of signal plus 99 molecules of pure noise.
- Define the affected set per arm, as EXP-3 §9.1's own logic implies ("identity on molecules without a constrained
  bond"):
  - K1 / K1_wrong / F_filter / F_flip: molecules with a specified structural bond;
  - K2_inf: every molecule with a TD-rotatable double bond (the late fold acts there);
  - K1_canonical: double or triple bond.
- The pre-declared non-inferiority subset stays the 102; the extra molecules are spliced from O. See ruling 1.

**S2. Exp 3 analysis: failures and the two conventions are not reported** (`experiments/exp3/analysis.py:73-87`,
`:105-113`).
- `paired` and `hierarchical` keep only finite pairs, i.e. the success intersection, with no failure count. EVAL §3
  requires both conventions, with failures per arm and seed.
- `derived_all_molecule_diff` divides by `n_all = len(first recall frame)` (`:106`), the number of generated
  molecules, not N = the csv rows (1000; §9.1 "(n_affected / N)").
- `ez_mols` and `plan_mols` come from O seed 0 alone (`:96-99`). A molecule that O seed 0 failed on disappears from
  every seed's analysis. Take the subsets from the key files.

**S3. GT-ensemble matching can truncate and miss GT states** (`flexitors/eval/stereo_eval.py:57-66`).
- `GetSubstructMatches(..., uniquify=False, maxMatches=2000)` runs on the **H-explicit** skeleton, so hydrogen
  permutations multiply the match count: each methyl contributes 3!, so 4 methyls give 1,296 per heavy-atom
  automorphism.
- Once the 2000 cap is hit, heavy-atom automorphisms that map a different cis/trans state can be cut. The GT state set
  is then incomplete, and correct generated conformers count as "inconsistent".
- The structural quadruples are all heavy atoms (each end has a heavy substituent, and heavy bonds precede H bonds in
  bond order). Fix: match the heavy-atom skeletons (assert that the quadruples are heavy), or flag any molecule whose
  match count reached the cap.

**S4. Exp 5 K3 non-inferiority uses the molecule bootstrap** (`experiments/exp5/analysis.py:199-200`, `:209-211`).
EVAL §7 states that non-inferiority claims use the hierarchical bootstrap (seeds, then molecules). Use
`alf.hier_ci` on the seed × molecule array, as Exp 3 already does for its NI test.

**S5. The baseline embedding differs from TD's protocol, and nothing checks it** (`generate.py:83-98` vs
`torsional-diffusion/generate_confs.py:64-69`).
- O uses per-conformer `EmbedMolecule`; TD uses `EmbedMultipleConfs(numThreads=5)`. Both are ETKDG, so this is
  acceptable in principle, but O is the reference of every Exp 3 and Exp 5 contrast and the absolute level quoted
  beside them.
- The ETKDG alkyne defect seen on the cluster (ruling 2) may depend on this call. My local check under RDKit
  2026.03.6, 150 `EmbedMolecule` embeddings of three terminal alkynes, gave sp C–H lengths of 1.04–1.08 Å and no
  missing bond.
- Before wave 2, run a cheap distributional golden: flexitors O (seeds 0–2) against the released R0–R2 on the 1000
  QM9 test molecules. Compare AMR-R and COV-R within the seed spread, failure counts, and the rate of the alkyne
  defect in R0's own seeds.

**S6. The validity-adjusted AMR-R becomes +∞** (`flexitors/eval/pool_metrics.py:135-136`).
- If every column of a molecule is invalid, `AMR-R_valid` = +∞ for that molecule, and every dataset mean in the Exp 5
  `levels` table (`analysis.py:212`) and the descriptive contrasts becomes ∞/NaN.
- With ≈ 15 % of conformers flagged by the alkyne defect, molecules where every seed carries it are likely.
- Pre-declare: a molecule with no valid slot is a validity failure (AMR excluded, COV 0), counted. See ruling 2.

**S7. The selector's torsion features are not heavy-atom dihedrals** (`flexitors/sampling/select.py:70-72`).
- `td_reference_quadruple` uses TD's first-neighbour rule, so c/d can be H on a heavy rotatable bond.
- The validated prototype used heavy-atom torsions (`round5/ledger/sample_data/signals.py:14-28`). H-terminated
  dihedrals add per-seed L noise to the FPS distance.
- Use heavy-first quadruples (the EXP-1 §10.1.2 rule), so the ported selector is the one that was validated.

## NIT

- **N1.** `K1_canonical` freezes unspecified double bonds at the seed's RDKit isomer (`stereo_guard.py:100`). EXP-3
  §9.2 rules that this isomer is worse than TD's own choice (E-CHEM-061). K1_canonical − K1 therefore confounds
  "freeze triple bonds" with "pin RDKit's random isomer". Say so beside that secondary row.
- **N2.** The smoke's K1_canonical row spliced into an O generated on `keys_double.txt` only ("arm molecules 20; O
  molecules 18"). The production table uses the full-pool O, so no defect in production; keep the smoke consistent.
- **N3.** `F_filter` on a molecule without targets regenerates O (`generate.py:179-181`). It is redundant, and adds
  GPU noise; S1 removes it.
- **N4.** The `tools` test count "15 passed" covers the three named files only; state the scope (see the
  `test_round2_model.py` note above).

## Rulings on the three raised issues

### Ruling 1: GPU nondeterminism (arm vs O on untouched molecules deviates 0.265 Å on GPU, 0.0 on CPU)

**The splice design of EXP-3 §9.1 is not invalidated.**
- **Unbiasedness:** arm and O run the same code with the same per-conformer streams. GPU float non-associativity is
  symmetric between them, so differences on affected molecules stay unbiased (the smoke's K1 − O on untouched
  molecules is −0.0004 Å AMR-R).
- **The untouched contribution is exactly right:** untouched molecules are *copied* from O, so their paired difference
  is exactly 0, which is the true effect of an identity operator.
- **What is lost is only the common-random-number variance reduction** on the molecules that are regenerated.

Required changes:
1. **Per-arm acting sets (S1), for QM9 and DRUGS.** This removes regenerated-but-untouched molecules from every arm,
   which is where most of the GPU noise enters (on QM9, K1 acts on 3 molecules, not 102).
2. **QM9 Exp 3 on CPU (option b), including its own O.**
   - QM9 Exp 3 is cheap on CPU (≈ 20 min per arm and seed on 8 cores; no GPU-h), and CPU gives exact pairing on the
     acting molecules too, so K2_inf's 102 molecules get the full CRN benefit.
   - The O of these contrasts must come from the **same CPU run configuration**, not from the GPU pool prefix: even on
     CPU the pool prefix equals a 2K run only to 1e-4 Å (batch-composition float32 effects), and 20 steps can amplify
     that.
   - All-molecule values stay derived as (n_affected / N) × subset difference, with absolute levels from one source
     (§9.1 condition 4).
3. **DRUGS on GPU, accepted** (CPU is too slow for 160 drug-size molecules × 2K × arms). Restate §9.1 condition 3 for
   GPU as two checks:
   - (i) the CPU bitwise identity test of the code path, which already passes;
   - (ii) a GPU distributional check on ≥ 50 unaffected molecules: arm − O AMR-R with a 95 % CI including 0 and
     |mean| < 0.002 Å.
   Report that DRUGS contrasts are unpaired draws (no CRN), so their CIs are wider. R-CHEM must sign this restatement
   off before wave 3.

Rejected options:
- **Deterministic CUDA.** `torch_scatter` and e3nn use atomic CUDA reductions with no deterministic implementation in
  torch 1.13. `torch.use_deterministic_algorithms(True)` would either raise or leave them nondeterministic; it is not
  attainable without rewriting kernels.
- **No splicing.** It costs more GPU, fixes no pairing, and adds noise on untouched molecules.

### Ruling 2: ETKDG terminal-alkyne C–H stretched in seeds (≈ 15 % of conformers invalid, O included)

1. **Root cause first** (cheap, before wave 1):
   - measure the defect rate in TD's own R0 seeds/outputs (`EmbedMultipleConfs`, cluster RDKit 2022.09.5) and in
     flexitors' per-conformer `EmbedMolecule`;
   - my 2026.03.6 check shows none.
   - If only flexitors' call produces it, fix the embedding call (it is then an S5 protocol difference, not a property
     of TD).
   - If R0 has it too, it is part of TD's L source and stays.
2. **Primary metrics: unchanged.** The evaluator is heavy-atom RMSD and does not see the H, and invalid conformers are
   never dropped (EVAL §3.1 (a), (c)).
3. **Any-cause invalid count (EVAL §3.1 definition): unchanged and reported.** A 1.37 Å C–H is a real geometric
   defect.
4. **Validity-adjusted secondary (EVAL §3.1 (b)): the mask should be method-induced invalidity only.** That means
   defects relative to the conformer's own seed: `n_missing_bonds_unexplained`, unexplained extra bonds, clashes the
   seed did not have, chirality/E/Z changes.
   - Reason: the secondary exists to penalise what a method breaks. A defect inherited unchanged by rigid TD moves is
     identical in every arm that shares seeds, so counting it just deletes ≈ 15 % of the slots in every arm and
     contaminates the secondary with an ETKDG artefact.
   - Report seed-inherited defects per arm beside it.
   - For Exp 1, where V2 can *repair* L, report "repaired" defects separately. The credit is legitimate, but it must
     not be presented as the validity-adjusted gain.
5. **No valid slot left** (S6): treat the molecule as a validity failure (AMR excluded, COV 0), counted.
6. This amends EVAL §3.1 (b), so **R-SAMPLE (EVAL owner) must sign it off.** Until then compute and store both masks
   (any-cause and method-induced), so no GPU work is redone either way.

### Ruling 3: K1 acts on only 3 QM9 molecules

1. **No change to the pre-declared QM9 family** (Holm over 8). The count was known before any result (E-CHEM-066;
   EXP-3 §9.3 calls QM9 confirmatory), so editing the family now would be a design change without information gain,
   and removing rows would weaken the pre-registration.
2. **Splice per acting set (S1).** On QM9, K1, K1_wrong, F_filter and F_flip then differ from O only on the 3
   specified molecules, and their contrasts are not drowned in regeneration noise.
3. **Report `n_acting` beside every row.** Label the QM9 K1, K1_wrong, F_filter and F_flip rows "uninformative by
   construction (operator acts on 3 molecules)"; no claim is drawn from them, including K1_wrong − K1 on QM9.
   QM9 K2_inf, through the late planarity fold on all 102, is the only informative QM9 arm (GT-consistency on 37,
   planarity on 102).
4. **EXP-3 §9.5 inheritance rule.** "DRUGS decides, or QM9 if DRUGS is missing" cannot work for K1 vs K2_inf, because
   QM9 cannot distinguish them. If the DRUGS arms are missing, defer the inheritance decision instead of falling back
   to QM9. R-CHEM to acknowledge this before the Exp 3 analysis runs.
