# EXP-1 M0 (revised, EXP-1 §10.3): oracle-velocity integration test of the V2 geometry sampler

**What it tests.** Can the planned sampler's L update follow an oracle field to the pure-L target? The update is
dihedral-preserving projected (P_W), torsion-restored and constraint-guarded. No model is trained; this is the CPU
gate the binding spec puts before any FlexiTors training.

**Data.** Validation paired cache only (`/scratch/nishanth.r/tordiff/data/QM9/cache_paired.val`, sha256
`1671196c0e1ea394…`). Read only. No test-set pair is used.

**Code.**
- `flexitors/geometry/projector.py`, `flexitors/geometry/l_step.py`, `experiments/exp1/m0_oracle.py`.
- Tests in `flexitors/tests/test_projector.py`. With the existing internal, torsion_types and dead-code suites:
  37 passed, 1 skipped.

**Job.** 46664, 2026-10-29 11:28 UTC / 16:58 IST, 16 CPUs. Outputs: `m0_summary.json` and `m0_pairs.csv.gz` (one row
per pair). An earlier attempt, 46574, stopped at the tests:
- a test helper named `setup` clashed with pytest's module-setup hook;
- the TD gateway needed `FLEXITORS_TD_ROOT` pointing at the frozen TD tree, which holds the torus tables.

Neither error touched the M0 code.

## Procedure (as specified)
For each pair (x_R = matched RDKit seed, x_G = its own GT) with `pair_ok`:
1. **Pure-L target x_G\*:** the W dihedrals of x_G (heavy-first rule, §10.1 item 2) are set to x_R's by side
   rotations, then x_G\* is Kabsch-aligned to x_R on heavy atoms.
2. **Integration:** N = 20 steps of x ← x + P_W(x)(x_G\* − x)/(N − n). After each step:
   - exact torsion restoration;
   - the guard: bonds within ±15 % of the seed, no non-bonded heavy pair below 0.75 × the covalent radii sum, and
     every potential tetrahedral stereocentre keeping its signed-volume sign;
   - halving of the step up to 3 times, then a skip.
3. **Kill rule (pre-declared):** recovery < 50 %, or > 5 % of primary pairs ending with a skipped step or a flipped
   centre.

## Results (78,158 pairs from 10,017 validation molecules)
| stratum | pairs | median RMSD R→G\* (Å) | median RMSD end→G\* (Å) | recovery | pairs with a skip | flipped at end | ring-dihedral RMSD > 10° (seed → end) |
|---|---|---|---|---|---|---|---|
| **primary** (≥ 1 heavy rotatable bond + ring ≥ 4) | 30,420 | 0.131 | 0.0005 | **99.6 %** | **4.74 %** | 0 | 41.7 % → 0.64 % |
| rigid (n_rot_heavy = 0), sanity | 5,785 | 0.124 | 0.0005 | 99.6 % | 3.27 % | 0 | 42.6 % → 0.90 % |
| flexible, no ring ≥ 4 (not in the gate) | 41,953 | 0.099 | 0.0002 | 99.8 % | 10.8 % | 0 | – |

- **W-dihedral drift:**
  - after restoration, max 2.7e-15 rad;
  - before restoration, per step: median 2.1e-5 rad, p95 2.7e-4 rad (primary).
- **W rows:** 316,915 rows used, 44,014 dropped as near-linear (mostly alkynes and nitriles), 48 % H-terminated.
- **Restoration off (ablation, 257 primary pairs, first conformer):**
  - final W drift up to 0.010 rad with restoration off, against ≈ 0 with it on;
  - recovery 99.5 % either way;
  - 8.2 % of these pairs had a skip (higher than the full primary rate, on a different, smaller set).

**Verdict by the pre-declared rule: PASS.** The margin on the skip criterion is small: 4.74 % against a 5 % limit.

## What it does and does not show
- **Confirmed:** with an oracle field, the projected, restored, guarded update reaches the pure-L target almost
  exactly. Ring-dihedral errors > 10° fall from 42 % of pairs to 0.6 %. W dihedrals are held to machine precision,
  and no stereocentre is flipped (the guard rejected 5,082 candidate steps on primary pairs). So the subspace
  P_W carries the RDKit → GT L error, and the sampler mechanics work.
- **Not shown:** that a network can learn this field. M0 has no learning; that is Gate 0 (M1).
- **Where the skips come from:**
  - 2.1 % of primary pairs have a GT target whose configuration at an RDKit potential stereocentre differs from the
    seed's (e.g. unspecified centres or N inversion). The guard then blocks the crossing by design.
  - The remaining 2.7 % (the skip share among pairs whose target stereo matches the seed) end with a skip for
    reasons not recorded per step: bond window, clash, or chirality on an intermediate geometry. Logging the reason
    per step is a needed fix before M1.
  - The flexible, no-ring stratum has more skips (10.8 %); it is outside the gate's population.
