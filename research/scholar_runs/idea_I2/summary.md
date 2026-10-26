# idea_I2 summary (Phase VIII iterations 2–4) — written by orchestrator from the search agent's hand-back, 2026-10-26
Coverage caveats: Google Scholar not searched; OpenAlex hit HTTP 429 mid-run; ChemRxiv/RSC full texts 403.

Closest prior art
1. Rowan pKa (Wagen & Wagen 2024, ChemRxiv 10.26434/chemrxiv-2024-8489b) — SAME OBJECTIVE + SAME MECHANISM except conformer source (ETKDG): per-protonation-state conformer search, AIMNet2 + CPCM-X, ensemble free energies, fitted to experiment; micro-pKa only (SAMPL6 RMSE 1.21). Note: paper_notes/rowan_pka_wagen2024.md
2. Pracht et al. 2018 SAMPL6 (abstract only) — same mechanism with DFT, macro-pKa.
3. Wieder et al. 2021 — tautomer counterpart; own ablation: finding the global minimum matters more than weighting many conformers.
4. Uni-pKa — Table 2 "Jaguar++" baseline weights 10 conformers; Uni-pKa still beats Rowan and Starling.
5. Starling (Wagen 2025) — same objective, learned energies (paper_notes/starling_wagen2025.md).
Only untested piece: a *generative* conformer model as the conformer source.
Attack: the pipeline already exists; any gain would be "better minima", not a new method.
