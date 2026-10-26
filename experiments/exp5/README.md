# Exp 5: selection at matched compute, Stage A (EXP-5 §12)

Pools (`flexitors/sampling/generate.py`): test seed 0 TD max(128, 4K) with the 128-column K-free block; seeds 1–2 TD 4K;
PG 2K and 4K per seed (batch 128) with the setting chosen on validation; validation TD 4K and PG 2K for the four README
settings. Each pool is scored once by the unchanged evaluator, featurised (torsions, MMFF94s, `flexitors/sampling/select.py`)
and screened for validity (`flexitors/eval/validity.py screen`).

Order: wave w1 (validation + TD test pools) → `python experiments/exp5/analysis.py validation --root $RES/round5/exp5
--csv $VAL_DIR/val_smiles.csv --out $RES/round5/exp5/analysis` (selector and PG setting; hashed decision file) → wave
w2 (PG test pools) → `analysis.py test`. Arms P0, P1, P2, P2b, P3, P6 and the ORACLE row are column subsets; the
primary Holm family {P3 − P1, P3 − P2b} × {AMR-R, COV-R@0.5} and the K3 rule are fixed in `analysis.py`.
