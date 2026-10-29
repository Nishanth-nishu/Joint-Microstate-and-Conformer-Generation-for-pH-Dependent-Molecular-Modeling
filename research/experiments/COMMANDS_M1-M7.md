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
