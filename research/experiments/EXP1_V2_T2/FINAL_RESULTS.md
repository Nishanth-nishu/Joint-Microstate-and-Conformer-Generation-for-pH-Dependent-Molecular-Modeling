# EXP1-V2-T2: final results (three training seeds, 20 epochs)

The final evaluation was run once, as job 57155. It started 2026-11-04 05:05 UTC / 10:35 IST and ended 08:27 UTC /
13:57 IST (elapsed 3 h 29 min, exit 0). Settings: 8 CPUs, 22 G, 6 workers (DESIGN.md addendum D).

Every number below is read from the job's outputs in `e1/`:
- `final_analysis.json` / `.log`;
- `geom_validate_seed{0,1}.json`;
- `*.time`, `timing_*.pkl.log.json`;
- `v2_diag_val.json`, `train_sacct.txt`.

Rules and analyses are those of DESIGN.md Part D and addenda B–D, all fixed before any T2 test number. Nothing was
re-run or re-selected after the results were seen.

MAT values are in Å; 1 mÅ = 0.001 Å. CIs are 95 % paired percentile bootstraps over molecules (5,000 draws).
"Seed-averaged" means the per-molecule metric averaged over training seeds 0, 1 and 2, then bootstrapped.

## 0. Checks before any number
| check | result |
|---|---|
| six runs complete (epochs 0–19, `best_model.pt`) | yes |
| checkpoints vs hashes written by the training jobs | 12 / 12 OK |
| E1 inputs vs `e1_manifest.sha256` | 24 / 24 OK |
| tests | 116 passed, 3 skipped |
| run-once guard | first and only run (`STARTED`, `DONE` present) |
| read-out errors | 0 of 5,390 conformers in each of the 12 read-outs |
| torsion drift after read-out (max over conformers) | ≤ 1.0 × 10⁻¹² ° in every read-out |

`best_model.pt` sha256 prefixes:

| | seed 0 | seed 1 | seed 2 |
|---|---|---|---|
| proj | 0d188b8e | 708c3b11 | f93244bf |
| noproj | 798f7bbe | 20a1fc48 | c33c6348 |

## 1. Denominators
| level | TD sampling seed 0 | TD sampling seed 1 |
|---|---|---|
| molecules attempted | 200 | 200 |
| molecules with TD output | 185 (TD fails on 15) | 185 |
| conformers read out per arm | 5,390 | 5,390 |
| paired molecules, all arms (RMSD defined everywhere) | **184** (1 dropped: NaN RMSD) | **184** |
| conformers in geometry checks | 5,388 | 5,388 |

The paired set is the same 184 molecules as in T1. ETKDG fails on 14 molecules, not 15, but the extra molecule is
outside the paired set.

## 2. Training (validation only)
| run | job(s) | epochs | best epoch | best val objective | train wall (sacct) | CPU time (sacct) | MaxRSS | non-finite loss |
|---|---|---|---|---|---|---|---|---|
| proj s0 | 49962 | 20 (0–9 = T1) | 19 | 0.2006 | 40.5 h (epochs 10–19) | 144.5 h | 16.5 G | 0 |
| noproj s0 | 49964 | 20 (0–9 = T1) | 19 | 0.2483 | 39.0 h | 139.8 h | 16.4 G | 0 |
| proj s1 | 49966 | 20 | 19 | 0.2022 | 75.5 h | 270.7 h | 16.4 G | 0 |
| noproj s1 | 49968 | 20 | 19 | 0.2517 | 75.0 h | 268.8 h | 16.2 G | 0 |
| proj s2 | 49970 | 20 | 19 | 0.2055 | 78.9 h | 279.6 h | 16.5 G | 0 |
| noproj s2 | 49972 + 52355 | 20 (resumed at 8) | 19 | 0.2522 | 31.9 + 58.2 h | 112.9 + 202.7 h | 16.4 G | 0 |

- **Selection:** every run's best epoch is its last, so validation was still improving at epoch 20.
  - Last six validation objectives, proj s0: 0.2121, 0.2116, 0.2040, 0.2047, 0.2058, 0.2006.
  - The same for noproj s0: 0.2640, 0.2593, 0.2553, 0.2585, 0.2538, 0.2483.
  - The slope is small and noisy: about −0.004 per epoch over the last five epochs.
- **Comparing the arms:** validation objectives cannot be compared between arms. Each arm divides its own loss by its
  own base: the projected loss for proj, the unprojected loss for noproj. The arms are compared on E1 only.
- **Validation diagnostic** (`v2_diag.py`, 1,000 validation samples): model / base 0.186–0.191 (proj s0–s2) and
  0.193–0.196 (noproj s0–s2). The ratio is similar across the five t quintiles.
- **Totals:** about 1,420 CPU hours for T2 training, epochs 10–19 of seed 0 included, plus T1's epochs 0–9 for seed 0.
- **Times:** sacct's wall-clock columns are in IST (cluster local time); the end times in FINAL_EXPERIMENT_STATUS.md
  are in UTC.

## 3. E1 accuracy (184 paired molecules)
### 3a. Means, TD sampling seed 0 (seed 1 in brackets)
| arm | MAT-R | MAT-P | COV-R@0.25 % | COV-R@0.5 % | COV-P@0.25 % | COV-P@0.5 % |
|---|---|---|---|---|---|---|
| ETKDGv3 | 0.229 (0.217) | 0.227 (0.221) | 69.5 (68.2) | 89.2 (91.1) | 67.9 (67.6) | 90.0 (91.8) |
| TD (frozen) | 0.170 (0.165) | 0.215 (0.206) | 82.8 (80.9) | 94.9 (96.2) | 71.4 (71.4) | 91.0 (92.9) |
| TD + MMFF94s (torsions restrained) | 0.163 (0.152) | 0.182 (0.170) | 82.5 (83.0) | 94.2 (96.3) | 78.8 (79.7) | 92.9 (93.8) |
| TD + L(proj), T1 10 epochs, seed 0 | 0.111 (0.106) | 0.145 (0.134) | 87.8 (88.0) | 95.7 (96.8) | 81.5 (83.0) | 93.8 (94.8) |
| TD + L(proj), T2 seed 0 / 1 / 2 | 0.107 / 0.106 / 0.106 | 0.139 / 0.139 / 0.139 | 88.2 / 88.3 / 88.0 | 95.9 / 95.9 / 95.9 | 82.2 / 82.1 / 81.9 | 93.9 / 94.1 / 94.1 |
| TD + L(noproj), T2 seed 0 / 1 / 2 | 0.107 / 0.106 / 0.107 | 0.141 / 0.141 / 0.140 | 87.8 / 88.4 / 88.5 | 95.9 / 95.9 / 95.8 | 81.7 / 81.9 / 82.3 | 93.8 / 94.0 / 94.2 |
| **TD + L(proj), seed-averaged** | **0.106 (0.103)** | **0.139 (0.129)** | 88.2 (88.5) | 95.9 (96.8) | 82.1 (83.5) | 94.0 (94.8) |
| TD + L(noproj), seed-averaged | 0.107 (0.104) | 0.141 (0.131) | 88.2 (88.3) | 95.9 (96.7) | 81.9 (83.1) | 94.0 (94.7) |
| **TD + constrained GFN2-xTB** | **0.097 (0.089)** | **0.123 (0.113)** | 88.3 (88.9) | 95.6 (96.6) | 83.3 (84.5) | 93.8 (94.8) |

### 3b. Paired differences, seed-averaged (negative MAT = first arm better)
| comparison | TD seed | MAT-R | MAT-P | COV-R@0.25 (points) | COV-P@0.25 (points) |
|---|---|---|---|---|---|
| proj − TD | 0 | **−0.0634 [−0.0724, −0.0552]** | −0.0755 [−0.0862, −0.0656] | +5.4 [+3.1, +8.0] | +10.7 [+7.7, +13.9] |
| | 1 | −0.0625 [−0.0708, −0.0548] | −0.0767 [−0.0870, −0.0667] | +7.6 [+4.5, +11.1] | +12.1 [+8.7, +15.6] |
| proj − noproj | 0 | **−0.0002 [−0.0013, +0.0011]** | −0.0015 [−0.0024, −0.0005] | −0.1 [−1.0, +0.6] | +0.1 [−0.4, +0.6] |
| | 1 | −0.0010 [−0.0022, −0.0001] | −0.0014 [−0.0024, −0.0006] | +0.1 [−0.6, +0.8] | +0.4 [−0.1, +1.0] |
| proj − xTB | 0 | **+0.0093 [+0.0036, +0.0141]** | +0.0167 [+0.0113, +0.0229] | −0.2 [−0.8, +0.5] | −1.2 [−2.3, −0.4] |
| | 1 | +0.0142 [+0.0096, +0.0199] | +0.0164 [+0.0116, +0.0219] | −0.4 [−2.0, +0.9] | −1.0 [−2.3, +0.1] |
| proj − MMFF | 0 | −0.0566 [−0.0662, −0.0472] | −0.0423 [−0.0512, −0.0339] | +5.7 [+3.2, +8.5] | +3.2 [+1.2, +5.7] |
| | 1 | −0.0496 [−0.0581, −0.0417] | −0.0407 [−0.0486, −0.0325] | +5.5 [+3.2, +8.1] | +3.8 [+1.9, +6.2] |
| proj − ETKDG | 0 | −0.1226 [−0.1411, −0.1046] | −0.0881 [−0.1056, −0.0702] | +18.7 [+14.3, +23.3] | +14.2 [+9.5, +19.0] |
| xTB − TD (reference) | 0 | −0.0726 [−0.0827, −0.0632] | −0.0922 [−0.1053, −0.0794] | +5.6 [+3.2, +8.2] | +11.9 [+8.8, +15.3] |
| T2 proj s0 − T1 proj s0 (20 vs 10 epochs) | 0 | **−0.0039 [−0.0061, −0.0021]** | −0.0053 [−0.0077, −0.0033] | +0.4 [−0.2, +1.2] | +0.8 [+0.1, +1.5] |
| | 1 | −0.0032 [−0.0050, −0.0016] | −0.0042 [−0.0062, −0.0023] | +0.7 [+0.1, +1.5] | +0.8 [+0.2, +1.6] |
| T2 noproj s0 − T1 noproj s0 | 0 | −0.0050 [−0.0072, −0.0030] | −0.0055 [−0.0078, −0.0033] | +0.0 [−0.7, +0.8] | +0.4 [−0.0, +0.9] |

### 3c. Per training seed
| | proj − TD (TD seed 0) | TD seed 1 sign | proj − noproj (TD seed 0) | proj − noproj (TD seed 1) | proj − xTB (TD seed 0) |
|---|---|---|---|---|---|
| seed 0 | −0.0628 [−0.0714, −0.0546] | −0.0625 | +0.0001 [−0.0012, +0.0015] | −0.0007 [−0.0030, +0.0009] | +0.0098 [+0.0046, +0.0144] |
| seed 1 | −0.0637 [−0.0729, −0.0556] | −0.0627 | −0.0002 [−0.0016, +0.0012] | −0.0014 [−0.0030, −0.0001] | +0.0089 [+0.0028, +0.0139] |
| seed 2 | −0.0636 [−0.0729, −0.0557] | −0.0624 | −0.0005 [−0.0020, +0.0012] | −0.0009 [−0.0017, −0.0001] | +0.0090 [+0.0030, +0.0140] |

## 4. Variability: three sources, reported separately (MAT-R, proj)
| source | TD seed 0 | TD seed 1 |
|---|---|---|
| training seed: SD of mean MAT-R over 3 seeds | 0.0005 Å (0.1069, 0.1060, 0.1061) | 0.0002 Å (0.1028, 0.1026, 0.1030) |
| TD sampling seed: seed-averaged mean MAT-R | 0.106 | 0.103 (shift ≈ 3 mÅ, same for every arm) |
| molecules: width of the bootstrap CI of proj − TD | 0.0172 Å | 0.0160 Å |

noproj training-seed SD: 0.0003 / 0.0002 Å. Molecule-to-molecule variation dominates. The training-seed SD is about
1 % of the effect over TD and about 3 % of the molecule-level CI width. These are not combined into one interval.

## 5. Where the gain is (geometry vs nearest reference, TD seed 0; seed 1 in brackets)
| arm | bond MAE Å | angle MAE ° | ring-dihedral MAE ° | acyclic-dihedral MAE ° |
|---|---|---|---|---|
| TD | 0.0319 (0.0319) | 3.30 (3.33) | 14.31 (14.36) | 40.1 (37.5) |
| MMFF (restrained) | 0.0194 (0.0194) | 2.19 (2.18) | 10.90 (11.05) | 40.2 (36.2) |
| L(proj) seed-avg | 0.0115 (0.0115) | 1.15 (1.13) | 8.58 (8.62) | 32.5 (28.0) |
| L(noproj) seed-avg | 0.0107 (0.0107) | 1.02 (1.00) | 8.39 (8.42) | 31.5 (28.3) |
| constrained xTB | **0.0027 (0.0025)** | **0.68 (0.63)** | **7.11 (7.50)** | **25.8 (22.4)** |

Paired differences, TD seed 0:

| comparison | bond (Å) | angle (°) | ring dihedral (°, n = 160) | acyclic dihedral (°, n = 139) |
|---|---|---|---|---|
| proj − TD | −0.0204 [−0.0216, −0.0193] | −2.15 [−2.27, −2.03] | −5.7 [−6.8, −4.8] | −7.7 [−11.5, −4.1] |
| proj − xTB | +0.0088 [+0.0083, +0.0093] | +0.47 [+0.41, +0.53] | +1.5 [+1.0, +2.0] | +6.7 [+3.6, +10.2] |
| proj − noproj | +0.0008 [+0.0005, +0.0011] | +0.13 [+0.10, +0.16] | +0.19 [+0.11, +0.26] | +0.9 [−0.7, +2.7] |

In words:
- **Against TD:** the head removes about 64 % of TD's bond-length error and 65 % of its angle error. xTB removes
  92 % and 79 %.
- **Projector:** noproj has slightly lower bond, angle and ring errors than proj; the CIs exclude 0.
- **Coverage:** proj − TD gains +5.4 points at 0.25 Å but only +1.0 [+0.0, +2.4] at 0.5 Å on TD seed 0
  (+0.6 on seed 1). The gain sits at tight thresholds.

## 6. Validity (conformers failing, of 5,388; TD seed 0, seed 1 in brackets)
| arm | identity fail | clash | stereo changed vs TD input |
|---|---|---|---|
| TD | 136 (141) | 22 (20) | 0 (0) |
| ETKDGv3 | 108 (102) | 1 (0) | 0 (0) |
| MMFF | 6 (9) | 8 (9) | 79 (79) |
| constrained xTB | **3 (2)** | **1 (1)** | 4 (0) |
| L(proj) s0 / s1 / s2 | 30 / 26 / 26 (27 / 33 / 31) | 6 / 6 / 6 (8 / 7 / 7) | 2 / 2 / 2 (0 / 0 / 0) |
| L(noproj) s0 / s1 / s2 | 23 / 22 / 34 (23 / 25 / 35) | 5 / 5 / 5 (6 / 6 / 6) | 2 / 2 / 2 (0 / 0 / 0) |

- **Stereo vs SMILES:** 270–272 conformers fail this check for every TD-derived arm (TD seed 1: 268), so the
  failures are inherited from TD's samples. MMFF raises it to 347.
- **Torsions:** TD's torsions are restored exactly by every arm (drift ≤ 10⁻¹² °). Constrained xTB before
  restoration drifts by a median 0.23° (p95 0.52°) on the timing subset.
- **Read-out guard, T2 proj:** 4.8–5.8 % of conformers had at least one skipped step. Mean step halvings
  1.14–1.40 per conformer. 0–1 conformers failed the final geometry check; 0 failed chirality.
- **xTB convergence:** xTB did not converge on 5 of 5,390 conformers.

## 7. Cost (CPU)
**Primary, like for like** (step 4b: the first 20 molecules of TD seed 0, 498 conformers, 6 workers, same job):

| method | CPU s / conformer | wall s |
|---|---|---|
| L(proj) read-out, T2 proj s0 | **1.23** | 116 |
| constrained GFN2-xTB | **0.30** | 26 |

The head is **4.1× more expensive** than constrained xTB per conformer.

**Secondary** (full read-outs; not like for like because worker counts differ):
- T2 heads: 1.06–1.11 CPU s / conformer at 6 workers. Wall time per read-out 966–1,012 s. Peak cgroup memory
  11.9–18.3 GB. 15–17 of 48 node cores busy.
- T1 heads: 1.64–1.67 at 32 workers.
- xTB: 0.36 at 32 workers.
- MMFF: 0.018 at 16 workers.

GPU cost of the read-out remains unmeasured (GPU quota exhausted throughout).

## 8. Exploratory (addendum C): strata, near-duplicate
These cannot change an endpoint verdict. There is no multiplicity correction, and every stratum here has n ≥ 24.

Seed-averaged MAT-R, TD seed 0 (seed 1 means in brackets):

| stratum | n | proj − TD | proj − xTB | proj − noproj |
|---|---|---|---|---|
| rotatable bonds 0 | 78 | −0.081 [−0.098, −0.066] (−0.078) | +0.013 [+0.000, +0.024] (+0.022) | +0.0015 [−0.0005, +0.0043] (+0.0000) |
| rotatable bonds 1 | 53 | −0.068 [−0.082, −0.055] (−0.068) | +0.008 [+0.003, +0.012] (+0.012) | −0.0001 [−0.0020, +0.0017] (−0.0026) |
| rotatable bonds ≥ 2 | 53 | −0.033 [−0.039, −0.027] (−0.034) | +0.005 [+0.002, +0.008] (+0.004) | −0.0028 [−0.0045, −0.0014] (−0.0009) |
| rings 0 | 24 | −0.030 [−0.038, −0.023] (−0.027) | +0.003 [+0.001, +0.007] (+0.004) | −0.0019 [−0.0034, −0.0001] (−0.0016) |
| rings 1 | 76 | −0.057 [−0.071, −0.046] (−0.054) | +0.003 [−0.010, +0.013] (+0.012) | +0.0001 [−0.0022, +0.0031] (−0.0008) |
| rings ≥ 2 | 84 | −0.078 [−0.093, −0.065] (−0.080) | +0.017 [+0.013, +0.021] (+0.019) | −0.0000 [−0.0013, +0.0012] (−0.0010) |
| 3- or 4-membered ring | 111 | −0.076 [−0.089, −0.064] (−0.074) | +0.014 [+0.011, +0.018] (+0.016) | −0.0005 [−0.0017, +0.0007] (−0.0012) |
| no 3- or 4-membered ring | 73 | −0.044 [−0.052, −0.036] (−0.046) | +0.002 [−0.011, +0.012] (+0.012) | +0.0002 [−0.0019, +0.0031] (−0.0007) |

What the strata show:
- **Gain over TD:** largest for rigid and polycyclic, strained molecules (where TD's RDKit local structure is
  worst); smallest for flexible and acyclic molecules.
- **Gap to xTB:** largest for polycyclic and small-ring molecules. Its CI includes 0 in two TD-seed-0 strata, with
  wide CIs.
- **Projector:** its only CI that excludes 0 in favour of proj is in the flexible (≥ 2 rotatable) and acyclic strata
  on TD seed 0: −2.8 mÅ and −1.9 mÅ. On TD seed 1, the ≥ 2 rotatable CI includes 0; acyclic still excludes it
  (−1.6 mÅ). This is exploratory and small.

**Near-duplicate removed** (`N#CC12C3C4CN1C2C43`, n = 183), TD seed 0:
- proj s0 / s1 / s2 − TD: −0.0629 / −0.0639 / −0.0637;
- proj − noproj: −0.0002 [−0.0014, +0.0012];
- proj − xTB: +0.0092 [+0.0036, +0.0140].

No conclusion changes.

## 9. Endpoint verdicts (computed mechanically by `t2_final_analysis.py`, rules of DESIGN.md Part D)
| endpoint | rule | result | verdict |
|---|---|---|---|
| 1 improves on TD | all 3 training-seed CIs < 0 and TD seed 1 same sign | 3 / 3 CIs < 0; seed 1 −0.062 to −0.063 | **yes** |
| 2 persists across seeds | report SD and seed-averaged effect | SD 0.5 mÅ; seed-averaged −0.0634 [−0.0724, −0.0552] | persists (SD ≪ effect) |
| 3 projector | seed-averaged CI < 0 and ≥ 2 / 3 seeds same sign | −0.0002 [−0.0013, +0.0011]; 2 / 3 seeds favour proj | **not shown** |
| 4 vs constrained xTB | CI < 0 beats; includes 0 indistinguishable; else worse | +0.0093 [+0.0036, +0.0141] | **worse** |
| 5 longer training | T2 s0 − T1 s0, paired CI | −0.0039 [−0.0061, −0.0021] | **helped** (small) |
| 6 where | MAE components, COV 0.25 vs 0.5 | local bonds / angles, ring dihedrals; tight-threshold coverage | descriptive |
| 7 cost | like-for-like CPU s / conformer | 1.23 vs 0.30 | head **4.1× more expensive** |

The interpretation, H1–H5 and the decision are in `FINAL_SCIENTIFIC_ASSESSMENT.md` and `RESEARCH_DECISION.md`.

## 10. Limitations
- **Scope:** one TD model (`qm9_default`), one E1 subset (200 molecules, 184 paired), GEOM-QM9 only. Nothing here
  speaks to GEOM-DRUGS or larger molecules.
- **References:** the GEOM-QM9 reference conformers are GFN2-xTB geometries, so constrained xTB has a built-in
  advantage on these metrics. That advantage is itself the finding.
- **xTB control:** the restraint force constant (k = 1.0 Eh) was not varied.
- **Cost:** the read-out was timed on CPU only, and GPU cost is unknown. CPU timings depend on node load: 15–18 of
  48 cores were busy during the read-outs.
- **Training:** 20 epochs, with validation still falling. The budget was fixed in advance and not extended.
- **Seeds:** three training seeds and two TD sampling seeds, for one TD checkpoint.
