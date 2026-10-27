# Tautomer candidates T1/T2 — intent (2026-10-27 10:35 UTC / 16:05 IST)
Claim T1 (diagnostic): released Uni-pKa, trained only on macro-pKa with an LSE microstate loss, predicts aqueous tautomer
ratios poorly, and its errors are concentrated where the loss gives the minor tautomer least gradient.
Claim T2 (method): fitting the same per-microstate free-energy head jointly to macro-pKa (Dwar-iBonD) and measured
Tautobase log K reduces out-of-fold tautomer log K MAE versus a pKa-only refit of the same head, with macro-pKa MAE on
SAMPL6/7/8 + Novartis non-inferior (≤ +0.05).
Falsifiers: T1 — zero-shot MAE < 1.0 log K on measured pairs; T2 — Δ(joint − pKa-only) tautomer-MAE CI includes 0, or
pKa non-inferiority fails.
Evidence base: scholar_runs/F_tautomer/summary.md; paper_notes/{unipka_2024,sphysnet_taut_pan2025,wieder_2021_tautomer_ratios}.md.
