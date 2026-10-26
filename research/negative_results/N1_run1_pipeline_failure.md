# N1 pilot run 1 — pipeline failure (not a scientific result)
Job 41159, 2026-10-26 21:42–22:24 UTC. Sanity gate failed: A3D-mean MAE 6.39/6.85/6.70/6.43 (SAMPL6/7/8/Novartis) vs
reported 0.49/0.55/0.62/0.81. Mean signed error −6.47 ± 0.87 → missing normalisation constant. All 10 "3D" conformers
identical for all 1,322 microstates → unipkainfer free-energy path conformer bug (EmbedMultipleConfs seed 0).
Lesson: pre-registered sanity gates caught both; always verify conformer diversity before ensemble ablations.
Raw outputs: research/experiments/N1_pilot/run1_failed_sanity_gate/.
