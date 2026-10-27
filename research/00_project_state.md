# Project state

Updated 2026-10-27 00:10 UTC (05:40 IST). Repo: Nishanth-nishu/Joint-Microstate-and-Conformer-Generation-for-pH-Dependent-Molecular-Modeling.
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
| X–XIV | not started — direction pending: (a) write up the negative/invariance result (outcome B), (b) pivot to tautomer populations (Uni-pKa Tautobase log K MAE 3.4), (c) test geometry with a physics-based ΔG source instead of Uni-pKa |

Tooling limitation: scholar-scout not installed on this machine; Google Scholar not searched (see CONVENTIONS.md).
