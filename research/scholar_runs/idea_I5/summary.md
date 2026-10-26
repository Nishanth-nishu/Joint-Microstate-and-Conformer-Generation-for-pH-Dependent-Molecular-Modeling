# idea_I5 summary (Phase VIII iterations 2–4) — written by orchestrator from the search agent's hand-back, 2026-10-26
Uni-pKa released code (commit 7cfbf65): 10 ETKDG+MMFF conformers + one 2D copy per microstate; evaluation = arithmetic mean of free energies across conformers, then log-sum-exp across microstates. Uni-pKa is NOT a 2D baseline.
Closest prior art
1. Uni-pKa — SAME OBJECTIVE, largely SAME MECHANISM (3D-conformer-conditioned microstate free energy fine-tuned on experimental pKa).
2. Starling (Wagen 2025) — SAME MECHANISM (log-sum-exp over 10 conformers per microstate); worse than Uni-pKa on 4/5 sets.
3. Wieder 2021 — SAME MECHANISM, tautomer half (ML potential retrained on experimental ratios via ensemble free energies).
4. sPhysNet-Taut 2025 — 3D single-conformer, fine-tuned on Tautobase (paper_notes/sphysnet_taut_pan2025.md).
5. Liu/Luo 2025 Angew (3D non-covalent descriptors for pKa) — UNCLEAR, not accessible; read first.
Novel residue: only a residual-head design with geometry on/off ablation.
Cheapest decisive test: Uni-pKa with its 2D-coordinate copy vs its 3D conformers.
Data: IUPAC zwitterion-capable ~1,656 molecules (434 with ≥15 heavy atoms); SAMPL6/7/8 ~70; Tautobase ~476 aqueous pairs. Expected effect ~0.5–1 kcal/mol ≈ label noise → paired, repeated scaffold splits; null result likely.
