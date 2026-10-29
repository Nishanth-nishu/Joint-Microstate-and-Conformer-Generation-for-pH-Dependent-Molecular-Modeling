# M3/M4: sequential baseline p(z | G, h) · p(x | z, G, h). E1 (state and conformer components) and E6 strata

The sequential baseline:
- parent + pH;
- Uni-pKa enumeration (simple or full SMARTS template, inclusive charge window −2…+2);
- released Uni-pKa free energies (5-fold mean);
- model populations;
- selected states;
- conformers (TD qm9_default with the E/Z lock default, or ETKDGv3), then validation.

All metrics are against external experimental data with scaffold-cluster bootstrap CIs. Nothing here was fitted.

## Runs
| job | when (UTC / IST) | what |
|---|---|---|
| 45304 | 2026-10-28 | end-to-end smoke (`smoke/`): 6 parents, glycine pKa 2.31/9.83, histamine 6.49/9.96 |
| 45306 | 17:37–21:44 / 23:07–03:14 | full E1 states (`e1/full/`): 845 enumeration seeds, 9,945 enumerated states, 4,264 scored |
| 45308 | 20:58–21:37 / 02:28–03:07 | E1 conformers (`e1_confs/`): 200 GEOM-QM9 test molecules, TD seeds 0/1, ETKDGv3 seeds 0/1, ETKDGv3+MMFF94s seed 0 |
| 45572 | 21:48–21:49 / 03:18–03:19 | E6 strata (`e6/`); N1 reproduction recomputed; conformer comparison v2 (NaN fix) |

## E1 state component: external macro-pKa (356 transitions, 270 scaffold groups; SAMPL6/7/8 + Novartis)
| arm | n | MAE [95% CI] | RMSE | slope |
|---|---|---|---|---|
| oracle states (dataset microstates) | 356 | 0.633 [0.563, 0.710] | 0.873 | 0.88 |
| enumerated, full template | 356 (0 failed) | 0.633 [0.562, 0.710] | 0.873 | 0.88 |
| enumerated, simple template | 281 (75 failed) | 0.650 [0.567, 0.741] | 0.902 | 0.86 |
| oracle on the same 281 | 281 | 0.642 [0.558, 0.733] | 0.897 | 0.86 |

- **Simple template:** all 75 failures are a missing charge level (e.g. "levels present [0, 1, 2], need 0,−1"). These
  are acid transitions the simple template has no SMARTS for. Recall of the dataset's microstates on successes: acid
  0.88, base 0.90.
- **Full template:** acid recall 0.99, base 1.00, and the same MAE as the oracle. Use the full template.
- **By source:** Novartis acids 0.82, Novartis bases 0.54, SAMPL6 0.52, SAMPL7 0.52, SAMPL8 0.64.
- **Reproduction against N1:** 356/356 transitions, max |Δ| 9.7e-10.
  - The first summary reported 4.34. That came from a key collision in N1's transitions.csv, not from scoring: two
    sampl8 records are both numbered idx 17.
  - Fixed by matching on (record, target) (`e1/full/reproduction_vs_N1_v2.json`).

## E1 tautomer component (Tautobase; released Uni-pKa, zero shot)
- Oracle log K MAE:
  - measured: 3.21 [2.81, 3.62], slope 0.10 (predict-0 MAE 3.53);
  - estimated: 3.21;
  - calculated: 3.19.
- **Enumerator recall of the partner tautomer: 0/221 measured (0/522 overall)** with either template. Uni-pKa's
  enumerator only moves protons between charge levels; it has no tautomer step. An enumerated sequential pipeline
  therefore cannot propose the second tautomer at all. This is a confirmed structural gap.

## E1 conformer component (neutral GEOM-QM9 test only: the only reference ensembles available)
183 molecules common to all arms. Each arm had 14–15 molecules with model failures (cage systems where ETKDG embeds
0 conformers) and 1 RMSD NaN.

| arm | MAT-R (Å) | MAT-P (Å) | COV-R@0.5 | COV-P@0.5 | COV-R@0.25 |
|---|---|---|---|---|---|
| TD seed 0 / 1 | 0.170 / 0.164 | 0.215 / 0.205 | 94.8 / 96.2 | 90.9 / 92.9 | 82.7 / 81.3 |
| ETKDGv3 seed 0 / 1 | 0.230 / 0.216 | 0.228 / 0.220 | 89.1 / 91.0 | 89.9 / 91.7 | 69.3 / 68.6 |
| ETKDGv3 + MMFF94s (ablation) | 0.209 | 0.159 | 88.2 | 94.7 | 73.4 |

Paired differences, arm minus mean of the TD seeds (bootstrap 95% CI):
- **ETKDGv3 seed 0:** MAT-R +0.063 [0.046, 0.080] Å; COV-R@0.5 −6.4 [−9.0, −3.9] pp; COV-R@0.25 −12.7 [−16.6, −8.7] pp.
- **ETKDGv3 seed 1:** similar.
- **ETKDGv3 + MMFF94s:** MAT-P −0.051 [−0.067, −0.034] Å (better precision); COV-R@0.5 −7.3 pp (worse recall).

Confirmed: TD beats unrefined ETKDGv3 on recall-type metrics for neutral in-distribution molecules. Nothing here
measures p(x | z) for charged or tautomeric states, because no reference ensembles exist for them.

## E6 strata (external pKa, oracle states; exploratory; groups with n < 20 flagged)
| stratum | n | MAE [95% CI] |
|---|---|---|
| acid-level charge +1 (base protonation) | 204 | 0.55 [0.48, 0.64] |
| acid-level charge 0 (neutral acid) | 151 | 0.74 [0.61, 0.88] |
| acid-level charge −1 | 1 | too few |
| non-zwitterionic | 302 | 0.62 [0.54, 0.70] |
| zwitterion-capable | 54 | 0.72 [0.52, 0.93] (mean signed +0.31) |
| 10–20 heavy atoms | 129 | 0.58 [0.48, 0.69] |
| > 20 heavy atoms | 226 | 0.66 [0.57, 0.77] |
| unspecified stereo | 158 | 0.59 [0.50, 0.70] |
| in Uni-pKa's Dwar-iBonD training families | 4 | 0.09 (small sample) |
| inside TD-qm9's domain (H,C,N,O,F, ≤ 9 heavy atoms) | 1 | too few |

Notes on the strata:
- **TD domain:** only 1 of 356 external transitions is inside it, so TD-qm9 cannot serve as p(x | z) for the external
  pKa set.
- **Tautomers by size:** measured Tautobase MAE is worse for ≤ 9 heavy atoms (4.01, n = 73) than for 10–20 (2.82,
  n = 146). It is likewise worse for pairs inside TD's domain (4.16, n = 66) than outside it (2.81).
- **E4 (pH generalisation): not run.** The data support only macro-pKa transitions (covered above), with no
  pH-resolved microstate populations or conformer labels.
