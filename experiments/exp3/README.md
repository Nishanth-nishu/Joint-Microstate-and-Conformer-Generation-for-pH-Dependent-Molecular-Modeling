# Exp 3: stereo-constrained torsional diffusion, Stage A inference arms (EXP-3 §9)

Arms: O, K1, K1_wrong, K1_canonical, K2_inf (sampler, `flexitors/sampling/stereo_guard.py`), F_filter (generation arm),
F_flip (`flexitors/sampling/stereo_fix.py`, CPU). QM9 test with the released `qm9_default`, 3 sampling seeds; O is the
2K prefix of Exp 5's TD pool of the same seed (per-conformer streams); arms are generated on the affected subsets
(102 molecules with a TD-rotatable double bond; 211 with a double or triple bond for K1_canonical) and spliced into O
(EXP-3 §9.1). GEOM-DRUGS (wave w3) after the staged download, seed 0 first, ≤ 4.5–5 GPU-h.

Rows: `experiments/common/stageA_tasks.tsv` (waves w2, w3), run by `experiments/common/infer.sbatch`.
Endpoints: `flexitors/eval/stereo_eval.py`. Pre-declared analysis: `analysis.py` (primary Holm family of 8 per
dataset, subset non-inferiority, secondary K1_wrong/K1_canonical, mixed-GT sensitivity).
