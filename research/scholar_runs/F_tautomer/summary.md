# F_tautomer summary (2026-10-27 10:30 UTC / 16:00 IST)

Closest prior work
- **Zankov et al. 2019 JCIM (abstract only)**: "conjugated" ridge regression and "twin" NN predicting log K_T and pKa
  jointly with the Kabachnik relation log K_T = ΔpKa built in; 639 log K_T + 2371 pKa, mixed solvents, ISIDA/CGR fragment
  descriptors. Reported that the relation "was respected only for conjugated but not individual models". → the generic
  idea "train pKa and tautomer constants jointly and thermodynamically consistently" is DONE (for QSPR).
- **Gunner et al. 2020 JCAMD; SAMPL6 overview 2021**: physical-chemistry statement that macro-pKa does not reveal tautomer
  ΔG ("there is no pKa connecting them"). → the chemistry half of the identifiability argument is known.
- Tautomer-only predictors with Tautobase numbers: sPhysNet-Taut 2025 (MAE 1.16 kcal/mol on a 100-pair test; fine-tuned on
  306 Tautobase pairs), MolTaut (MAE 2.61), Wieder 2021 (RMSE 2.8 after fit), Pan 2026 (classification). Kazan QSPR
  (Gimadiev 2018 RMSE 0.7 log units CV; multi-solvent).
- pKa-trained microstate models (Uni-pKa, AcepKa, Starling, DTi-pKa, Epik 7): no tautomer data in training found; Uni-pKa
  reports Tautobase log K MAE 3.435 (subset/protocol unknown). No published analysis of why.

Not found (62 operations; Google Scholar not searched; ACS/ChemRxiv 403)
- A per-microstate free-energy network trained with BOTH an LSE macro-pKa loss and experimental tautomer-ratio loss.
- An ML-identifiability analysis: an LSE macro-pKa loss gives the minor tautomer's free energy a gradient proportional to
  its population, so within-macrostate ΔG is (near-)unidentified by macro-pKa data.

Open unknowns: Zankov 2019 full text (micro vs macro targets, data availability); Uni-pKa SI (Tautobase protocol);
Starling (tautomer eval?); Epik 7 tautomer training data; Tautobase license; G-pKa dataset contents.
