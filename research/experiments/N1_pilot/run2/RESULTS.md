# N1 pilot run 2 — results (job 41194, 2026-10-26 22:30–23:58 UTC / 2026-10-27 04:00–05:28 IST; re-analysis job 41255)

Pre-registered verdict (gate referenced to Uni-pKa Table 1, PREREG amendment 3 — post hoc, see PREREG):
**N1 FALSIFIED on this pilot.**

## Sanity gate (A3D-mean = default Uni-pKa, MAE; Table 1 reference in brackets) — all pass (±0.10)
Novartis acid 0.818 [0.810], base 0.544 [0.493], total 0.654 [0.620]; SAMPL6 0.516 [0.554]; SAMPL7 0.522 [0.570];
SAMPL8 0.638 [0.619/0.631]. RMSE also tracks (acid 1.106 [1.061], SAMPL8 0.885 [0.878]). 356/356 transitions, 0 failures,
0 microstates where a 3D slot fell back to the 2D copy.

## Primary endpoint (stratum_v2; molecule-paired bootstrap 10k, 95 % CI)
| stratum | n | MAE A3D-only-mean | MAE A2D | Δ = MAE(2D) − MAE(3D) |
|---|---|---|---|---|
| sensitive | 185 | 0.664 | 0.669 | +0.005 [−0.014, +0.024] |
| control | 62 | 0.557 | 0.551 | −0.006 [−0.042, +0.029] |
| other | 109 | 0.630 | 0.610 | −0.020 [−0.054, +0.012] |
DiD (sensitive − control) = +0.010 [−0.030, +0.051]. Both CIs include 0 → falsified by the pre-registered rule.
Power: sensitive n=185 ≥ 30; the CI upper bound (0.024) excludes the pre-stated expected effect range (0.05–0.2).

## Secondary
- Aggregation is irrelevant: sensitive Δ(LSE − mean) = +0.0001 [−0.001, +0.001]; Δ(min − mean) = +0.003 [−0.004, +0.010];
  Δ(single conformer − mean) = 0.0000 [−0.004, +0.004].
- The model's output barely moves with the conformer: median per-conformer pKa SD 0.029 (sensitive), 0.024 (control/other).
- Conformer spread vs |error| (Spearman): sensitive ρ = 0.11 (p = 0.13); control −0.22 (p = 0.08); other 0.08 (p = 0.38).

## Interpretation (scope-limited)
Released Uni-pKa is effectively conformer-invariant at the pKa level: 10 MMFF conformers, one conformer, or 2D
coordinates give the same MAE within ~0.02 units, on sensitive and control molecules alike. This is a statement about
the model, NOT evidence that conformation is physically irrelevant to pKa — the model may simply have learned to
ignore geometry (training samples one random conformer per epoch). The stratum is a topological heuristic, and the test
sets are small/drug-like. It removes the motivation the project had for I1 from this model class.
