# Project state

Updated 2026-10-27 10:55 UTC (16:25 IST). Repo: Nishanth-nishu/Joint-Microstate-and-Conformer-Generation-for-pH-Dependent-Molecular-Modeling.
Work dir on gnode118: /scratch/nishanth.r/jmcg (staging on /home2/nishanth.r/jmcg_stage, rsynced into scratch by jobs).

| Phase | Status |
|---|---|
| 0 Fresh repo from TorDiff-NAX | done, pushed (44e5ed1), job 40060 |
| I Baseline audit | done (01_baseline_audit.md); no code changed |
| II Literature families A–E | round 1 done (157 queries, 02_literature_search_log.csv) |
| III Citation expansion | thin — to extend |
| IV LoQI/RAPA audit | done round 1 (04_loqi_rapa_audit.md; LoQI PDF inaccessible) |
| V Feature matrix | conformer methods done (05_feature_matrix_conformer_methods.csv); needs AcepKa/Uni-pKa/LoQI/RAPA columns |
| VI Dataset audit | round 1 done (06_dataset_audit.md) |
| VII Candidate ideas | initial set (07_candidate_ideas.md) |
| VIII Novelty loop | round 1 done: I2 Reject, I4/I5 Revise → N1 (08_novelty_claims.csv) |
| IX Minimum experiment | N1 pilot done: **N1 falsified** (run2/RESULTS.md). Released Uni-pKa is conformer-invariant; sanity gate passes vs Table 1 (post-hoc reference fix, amendment 3) |
| VIII-b Novelty loop round 2 (tautomers) | done: F_tautomer search + judge → T1 B (diagnostic), T2 B (joint vs tautomer-only) |
| IX-b Tautomer pilot | pre-registered (experiments/T_pilot/PREREG.md); code ready; awaiting user approval to run (CPU, ~7 h) |
| X–XIV | not started |

Tooling limitation: scholar-scout not installed on this machine; Google Scholar not searched (see CONVENTIONS.md).
