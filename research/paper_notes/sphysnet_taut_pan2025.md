# sPhysNet-Taut (Pan, Zhang, Xia, Zhang, JCTC 2025) — 3D NN tautomer ratios fine-tuned on experiment

## Citation
Pan X., Zhang X., Xia S., Zhang J. Z. H.(Y.). "Fast and Accurate Prediction of Tautomer Ratios in Aqueous Solution via a
Siamese Neural Network." J. Chem. Theory Comput. 21(6):3132–3141 (2025). DOI 10.1021/acs.jctc.5c00041 (PMC11948319).
Checked 2026-10-26 UTC.

## Fulltext status
full_text_read (PMC, via WebFetch summariser).

## Method
- Base sPhysNet-MT pretrained on Frag20-Taut DFT (B3LYP/6-31G*/SMD) gas/water/transfer energies (~0.9 M samples).
- Geometry: 300 ETKDG conformers per tautomer, MMFF94-optimised, "retaining the lowest-energy optimized conformation" — single
  conformer per tautomer; no Boltzmann averaging tested (not described).
- Siamese: shared base network on both tautomers, ΔG = difference of predicted aqueous energies.
- Fine-tuned on experimental Tautobase aqueous ratios (483 records; 306 train / 77 val, 5-fold; fixed 100-pair test).

## Results
100-tautomers test: MAE 1.16, RMSE 1.93 kcal/mol (MolTaut 3.64; B3LYP/SMD 6.27). SAMPL2 RMSE 1.0 kcal/mol. Fine-tuning on
75 pairs → RMSE 2.80, better than from scratch.

## Relevance to I2/I5
I5 (tautomer half): same objective, overlapping mechanism — a 3D-geometry-conditioned free-energy model whose
experimental-data fine-tuning is effectively a learned correction to a QM-pretrained energy. Differences: single MMFF
minimum, not an ensemble; no 2D-without-geometry ablation; no pKa. Also shows the experimental tautomer set is ~480 pairs.
