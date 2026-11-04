# Commands: master-prompt milestones M1–M7 (all on gnode118, partition/account plafnet2)

Staging is `/home2/nishanth.r/jmcg_stage`; every job rsyncs `jmcg/` into `/scratch/nishanth.r/jmcg`.
- **Python environments:** the jmcg venv (Python 3.10, rdkit 2026.09.1) for data, seq and joint; the TD venv for
  condgen and TD.
- **NFS caching:** wait about 60 s after editing before submitting. `joint_sample.sbatch` additionally rsyncs with
  `--checksum` and logs the md5 of the source it ran.

| step | command | job |
|---|---|---|
| M2 canonical data v1 | `sbatch slurm_jmcg/m2_data.sbatch` | 45295 |
| M3 smoke (parents + pH → states → conformers) | `sbatch slurm_jmcg/m3_seq_smoke.sbatch` | 45304 |
| M4 E1 states, full | `sbatch slurm_jmcg/m3_e1.sbatch` | 45306 |
| M4 E1 conformers, n = 200 | `sbatch slurm_jmcg/e1_confs.sbatch` | 45308 |
| E2 features v1 (superseded) | `sbatch slurm_jmcg/e2_features.sbatch` | 45307 |
| E2 bug diagnosis | `sbatch slurm_jmcg/e2_debug_arom.sbatch`; `sbatch slurm_jmcg/e2_debug_arom2.sbatch` | 45568, 45570 |
| E2 features v2 + pre-registered fit + E6 + N1 reproduction + conformer comparison v2 | `sbatch slurm_jmcg/e2_fit.sbatch` | 45572 |
| E2 exploratory robustness (post hoc) + refactor check | `sbatch slurm_jmcg/e2_robust.sbatch` | 45574 |
| M5 joint sampler demo (S0/S1/J1, pH 7.4) | `sbatch slurm_jmcg/joint_sample.sbatch` | 45584 |

Module entry points (jmcg venv unless noted):
```
python -m jmcg.data.build_canonical ...                   # M2
python -m jmcg.data.audit --data data_manifests/canonical_v1
python -m jmcg.seq.state_stage --parents parents.csv --ph 2 7.4 12 --out <dir>
python -m jmcg.seq.e1_states --data data_manifests/canonical_v1 --out <dir> --n1_transitions <N1 transitions.csv>
python -m jmcg.seq.e6_strata --data data_manifests/canonical_v1 --e1 <E1 dir> --out <dir>
python -m jmcg.condgen.e1_confs etkdg|compare ...          # TD venv
python -m jmcg.joint.e2_tautomer features|fit|robust --data data_manifests/canonical_v1 --features <csv> --energies <E1 energies.csv> --out <dir>
python -m jmcg.joint.sample --state_dir <state stage dir> --coef <e2_result.json> --arm S0|S1|J1 --ph 7.4 --out <dir>
```

## 2026-10-29: P1 unminimised features (exploratory), EXP-1 M0
| step | command | job |
|---|---|---|
| P1 panel + mmff/unmin features + TD (16 chunks) + compare + sampler | `sbatch slurm_jmcg/p1_unmin.sbatch` | 46561 (TD chunk 1 OOM-killed; compare KeyError psi1) |
| P1 resume: TD chunk 1 alone, td_features, compare, unmin sampler | `sbatch slurm_jmcg/p1_unmin_resume.sbatch` | 46662 (failed on a grep exit code), 46673 |
| EXP-1 M0 oracle-velocity test (TD venv, validation cache) | `sbatch slurm_jmcg/exp1_m0.sbatch` | 46574 (test-name clash), 46664 |
| read-only inspections (flexitors, specs, paired cache) | `slurm_jmcg/inspect_*.sbatch` | 46545, 46552, 46565, 46569 |
```
python -m jmcg.joint.p1_unmin panel|extract|td_features|compare ...      # jmcg venv; see module docstring
python -m jmcg.joint.sample ... --feature_path unmin --coef <compare>/coef_unmin.json
FLEXITORS_TD_ROOT=$TD/torsional-diffusion python experiments/exp1/m0_oracle.py --cache $TD/data/QM9/cache_paired.val --out <dir>
```

### EXP-1 V2 smoke (2026-10-29)
| step | command | jobs |
|---|---|---|
| GPU smoke (held, cancelled: GPU quota) | `sbatch slurm_jmcg/exp1_v2_smoke.sbatch` | 46725 |
| CPU smoke: tests, resume check, 400-iteration training | `sbatch slurm_jmcg/exp1_v2_smoke_cpu.sbatch` | 46738, 46739 (failed), 46744 |
| smoke checks only, on the 46744 run | `sbatch slurm_jmcg/exp1_v2_smoke_checks.sbatch` | 46848 |
| overfit diagnostic (32 mols, w_tau 0/1) | `sbatch slurm_jmcg/exp1_v2_overfit.sbatch` | 46855 |
| re-smoke, `--geo_norm base` + secondary overfit | `sbatch slurm_jmcg/exp1_v2_resmoke_norm.sbatch` | 46998 |
| re-smoke 3, `--l_head_init zero` (+ cancelled 2,000-iteration secondary) | `sbatch slurm_jmcg/exp1_v2_resmoke_zero.sbatch` | 47095 (test fail), 47102 |
| EXP1-V2-T1 training, one arm | `sbatch --export=ALL,ARM=proj\|noproj --job-name=jmcg_v2_t1_<arm> slurm_jmcg/exp1_v2_t1_train.sbatch` | 47158/47159 (1 h default limit), 47246/47247 |
| E1 read-out pipeline dry run (3 molecules, smoke checkpoint, no metrics) | `sbatch slurm_jmcg/exp1_t1_e1_dryrun.sbatch` | 47163 |
| EXP1-V2-T1 E1 evaluation | `sbatch --dependency=afterok:47246:47247 slurm_jmcg/exp1_v2_t1_e1.sbatch` | 47249 |
| leakage check, E1 test molecules vs paired caches | `sbatch slurm_jmcg/exp1_t1_leak_check.sbatch` | 49875 |
| GPU availability probe (held: AssocGrpGRESMinutes, cancelled) | `sbatch --gres=gpu:1 --wrap "nvidia-smi -L"` | 49881 |
| addendum A: constrained xTB control + fixed L read-out, E1 | `sbatch slurm_jmcg/exp1_v2_t1_xtb.sbatch` | 49899 |
| review follow-up: raw W drift of xTB / MMFF before restore, seed 0 | `sbatch slurm_jmcg/exp1_t1_xtb_signcheck.sbatch` | 49940 |
| read-only diff of a would-be commit vs HEAD | `sbatch slurm_jmcg/review_diff.sbatch` | 49935 (staging-only, misleading), 49936 |
### EXP1-V2-T2 (2026-10-31)
| step | command | jobs |
|---|---|---|
| audit of saved T1 E1 outputs | `sbatch slurm_jmcg/exp1_t2_audit.sbatch` | 49949 (print bug), 49954 (SMILES-key bug), 49955 |
| OOM probe, seed-0 re-run with oom_kill monitor | `sbatch slurm_jmcg/exp1_t2_oom_probe.sbatch` | 49957 |
| geometry validation / decomposition, all arms | `sbatch slurm_jmcg/exp1_t2_geom_validate.sbatch` | 49960 |
| T2 training, one run (+ resume if not ok) | `sbatch --export=ALL,ARM=<proj\|noproj>,SEED=<0\|1\|2> --job-name=jmcg_v2_t2_<arm>_s<seed> slurm_jmcg/exp1_v2_t2_train.sbatch`; resume: add `RESUME=1`, `--dependency=afternotok:<job>` | 49962–49972 (even), 49974–49979 |
| E1 input manifest + analysis dry run (T1 stand-ins) | `sbatch slurm_jmcg/exp1_t2_manifest.sbatch` | 51353 |
| addendum-C analysis dry run (geometry stand-ins) | `sbatch slurm_jmcg/exp1_t2_analysis_dryrun.sbatch` | 53494 |
| noproj s2 moved to 6 CPUs (cancel at epoch boundary, resume) | `sbatch --cpus-per-task=6 --export=ALL,ARM=noproj,SEED=2,RESUME=1 --job-name=jmcg_v2_t2_noproj_s2 slurm_jmcg/exp1_v2_t2_train.sbatch` | 52355 (+ resume 52356) |
| final evaluation, run once (8 CPUs, 22 G, 6 workers; addendum D) | `sbatch slurm_jmcg/exp1_v2_t2_e1.sbatch` | 51354 → 56777 (self-resubmitted, cancelled unstarted), 57155 |

Reproduce the final analysis from saved outputs (inside a job; `/scratch` is not readable on the login node):
`python experiments/exp1/t2_final_analysis.py --n200 $P/results/e1_confs/n200 --t1 $P/results/exp1/T1_e1 --t2 $P/results/exp1/T2_e1 --geom_dir $P/results/exp1/T2_e1 --train_dir $P/workdir/exp1 --out <new dir>`
(`P=/scratch/nishanth.r/jmcg`; the run-once guard is only in the sbatch script, the analysis itself is read-only).
