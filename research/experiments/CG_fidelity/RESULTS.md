# CG wrapper fidelity: our TDGenerator vs TD's own generate_confs.py

2026-10-28. Job 44897 on gnode118 (CPU only, 8 cores), 10:41–11:00 UTC / 16:11–16:30 IST. Smoke run (N=3): job 44894. Inventory: job 44891.

**Verdict: FAITHFUL** under the rule pre-registered in `jmcg/condgen/fidelity.py` (rule written before the 50-molecule run). On 50 GEOM-QM9 test molecules, our wrapper's coverage and matching RMSD differ from TD's own script by less than TD's own seed-to-seed difference.

## Setup
- **Molecules:** 50 of the 1000 GEOM-QM9 test rows (`data/QM9/test_smiles.csv`), drawn with `numpy.default_rng(0)`. Row indices are in the job log and `n50/subset.csv`.
- **Checkpoint:** `pretrained/qm9_default/best_model.pt` (sha256 6a9a6ee9b0e51bda…), pinned by path. Smoke runs 1 and 2 used `qm9_1order`; this run does not.
- **Conformer budget:** 2 × n_true per molecule, TD's default. 20 SDE steps, no MMFF.
- **Arms:**
  - R0 and R1: TD `generate_confs.py --seed 0/1 --no_energy`.
  - W0 and W1: `TDGenerator.generate(Chem.AddHs(MolFromSmiles(corrected_smiles)), 2n, seed)`, i.e. the same Mol that TD builds.
  - The state layer was not used for generation. Its verdict on each molecule was only recorded.
- **Scoring:** all four arms were scored with TD's `evaluate_confs.py --dataset qm9` against `test_mols.pkl`.

## Results (n_evaluated = 46 in every arm; the 4 failures count as 0 coverage in the summary means)

| arm | COV-R@0.5 % | MAT-R Å | COV-P@0.5 % | MAT-P Å | gen time (s, 50 mols) |
|---|---|---|---|---|---|
| R0 (TD, seed 0) | 84.79 | 0.1919 | 80.07 | 0.2296 | 299 |
| R1 (TD, seed 1) | 81.81 | 0.2033 | 79.88 | 0.2379 | 299 |
| W0 (wrapper, seed 0) | 87.25 | 0.1820 | 82.50 | 0.2214 | 273 |
| W1 (wrapper, seed 1) | 82.95 | 0.2015 | 80.90 | 0.2359 | 273 |

Paired per-molecule differences, n = 46, with bootstrap 95% CIs:

| pair | ΔMAT-R Å | ΔCOV-R@0.25 pp | ΔCOV-R@0.5 pp |
|---|---|---|---|
| noise: R1 − R0 | +0.0114 [−0.0103, +0.0322] | +3.80 [−0.35, +9.37] | −3.23 [−7.23, −0.08] |
| W0 − R0 | −0.0100 [−0.0431, +0.0173] | +5.38 [−1.80, +13.93] | +2.68 [−1.41, +8.20] |
| W1 − R1 | −0.0019 [−0.0240, +0.0178] | −2.00 [−6.92, +3.21] | +1.24 [−2.44, +5.07] |
| W1 − W0 | +0.0195 [−0.0157, +0.0622] | −3.58 [−10.65, +3.43] | −4.67 [−11.20, +0.41] |

**How the rule applies:**
- The MAT-R tolerance is 0.0114 + 0.02 = 0.0314 Å. Both wrapper pairs are within it.
- The COV-R@0.5 tolerance is 3.23 + 5 = 8.23 pp. Both wrapper pairs are within it.
- The wrapper has 4 failures and TD has 4, so the failure condition holds.

All wrapper-vs-TD CIs include 0. The noise pair's COV-R@0.5 CI just excludes 0, so with 46 molecules, seed-to-seed differences of about 3 pp are normal.

## Failures (identical molecules in all four arms)
| molecule | TD script | wrapper |
|---|---|---|
| `C1OC[C@H]2[C@H]3CN=C1N32` | Failed to embed | RuntimeError: ETKDG seed embedding failed |
| `C[C@H]1[C@@H]2[C@H]3C[C@](C)(C3)[C@H]12` | Failed to embed | same |
| `C[C@]12C[C@](CCO)(C1)C2` | Failed to embed | same |
| `N#N.[NH]C1=N[CH]OC1=O` (two fragments, radicals; a corrupted GEOM entry) | Failed to get seed (rejected cleanly) | **bare AssertionError** inside TD's featuriser |

The three embedding failures are strained polycyclic cages, where RDKit's default ETKDG fails. TD and the wrapper fail on them in the same way. In the wrapper the fourth molecule raises a bare AssertionError instead of a clear message. In normal use `build_state` rejects it first ("several fragments"), so this only happens when the generator is called directly. It is a small robustness gap, not a fidelity problem.

## State layer on the same 50 molecules (recorded only)
- 44 accepted.
- 5 rejected for unassigned stereo elements: the GEOM SMILES leave some stereocentres or double bonds unspecified, and our "require" policy refuses to pick one silently.
- 1 rejected as multi-fragment (the corrupted entry above).
- 1 accepted molecule has 0 rotatable bonds, so the model is not applied to it (the same as TD).

## What this shows and what it does not
- It shows: on neutral in-distribution QM9 molecules with the default checkpoint, the wrapper's changes do not measurably degrade TD at n = 46. Those changes are per-conformer seeded embedding, `EmbedMolecule` per conformer, batch size equal to n, and keeping atom order.
- It does not show:
  - anything about charged (B) or larger (C) states;
  - differences smaller than about 0.03 Å MAT-R or about 8 pp COV-R, which this n cannot resolve;
  - anything about `qm9_1order`, the checkpoint used in smoke 1 and smoke 2. Those runs should be repeated with `qm9_default` pinned before their TD outputs are cited.
- No claim is made about TD vs ETKDG.

Files: `n50/compare.json` (rule, verdict, CIs), `n50/compare_per_mol.csv`, `n50/W*.pkl.log.json` (per-molecule wrapper records), and the generation and evaluation logs. Smoke run: `n3/`.
