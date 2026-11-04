# EXP1-V2-T1: controlled training of the V2 geometry head + E1 evaluation (pre-declared 2026-10-29, before training)

Status: **exploratory, reduced budget** (CPU; the plafnet2 GPU-minute quota is exhausted, job reason
`AssocGrpGRESMinutes`). It is not EXP-1 Gate 0 (spec: 3 σ_b × 30 epochs on GPU) and makes no GPU-performance claim.

## Gate that opened this
Smoke 3 (job 47102, `EXP1_V2_smoke/v2_resmoke_zero_s0/`): S1–S5 PASS. Validation geometry loss 0.00831 vs Δ̂ = 0 base
0.01183 (0.703 × base, S2 needs < 0.9). Smoke configuration history: raw loss (46744) fail, `--geo_norm base` (46998)
fail, + `--l_head_init zero` (47102) pass, one change at a time, criteria never changed.

## Question
Does a learned, dihedral-preserving local-geometry correction, applied to frozen-TD conformers, move them closer to
the GEOM-QM9 reference ensembles than TD alone, and is the projector (the principal architectural change) needed?

## Design
Training: geometry objective only (w_τ = 0), as EXP-1 M1. Torsions in the evaluation come from frozen TD
`qm9_default`, so the torsion model is held fixed and identical across arms.

| arm | training flags (else identical) | checkpoint dir (own path) |
|---|---|---|
| T1-proj (V2) | `--geo_projector on` | `$P/workdir/exp1/T1_proj_s0` |
| T1-noproj (ablation of the projector) | `--geo_projector off` (loss with P_W = I; read-out step unprojected; torsion restoration and guard unchanged) | `$P/workdir/exp1/T1_noproj_s0` |

Common: `--w_tau 0 --w_geo 1 --geo_norm base --l_head_init zero --sigma_b 0.15 --seed 0 --lr 1e-3 --batch_size 32
--n_epochs 10 --limit_val_iters 100 --num_workers 4`, TD qm9 hyper-parameters, round-2 paired cache (82,895 train,
10,376 val molecules). Checkpoint = `best_model.pt`, lowest validation geometry objective (validation split only).
Fixed budget of 10 epochs; nothing is chosen on the E1 (test) molecules.

The w_geo = 0 ablation (A3) is not trained: with w_τ = 0 it is an untrained head, i.e. the TD arm below.

## E1 evaluation (run once per finished checkpoint)
Molecules: the 200-molecule E1 subset (`M3_seq/e1_confs/subset.csv`, GEOM-QM9 **test**), paired on the molecules common
to all arms (183 in the reference run). Budget 2 × n_true; scored with TD `evaluate_confs.py`; paired bootstrap over
molecules (`jmcg.condgen.e1_confs compare`, reference = TD).

| arm | conformers |
|---|---|
| TD | existing TD seed 0 and seed 1 samples (job 45308), unchanged |
| TD + L(proj) | the same TD conformers, each refined by the τ-frozen read-out (20 steps, torsion restoration, guard) of T1-proj |
| TD + L(noproj) | the same, with T1-noproj and unprojected steps |
| TD + MMFF94s (non-learned control) | the same TD conformers, MMFF94s minimisation (≤ 200 iterations), W dihedrals then restored to TD's |
| ETKDGv3 | existing seed 0 and seed 1 samples, unchanged |

**Primary endpoint:** MAT-R difference (TD + L(proj)) − TD, seed 0, paired bootstrap 95 % CI.
**Secondary:** COV-R@0.25/0.5, MAT-P, COV-P@0.25/0.5; seed 1 replicate; projector ablation (proj − noproj); the MMFF control.
**Also recorded:** training curves (train/val geometry loss and base, gradient norms), guard halvings/skips/chirality
rejections per conformer, W-dihedral drift after read-out (torsions must be unchanged), read-out runtime per
conformer (CPU), evaluate_confs failures and RMSD-NaN drops.

## Interpretation rules (fixed now)
- "Improves on TD" only if the primary CI excludes 0 in the favourable direction and seed 1 has the same sign.
- The projector "matters" only if proj − noproj on MAT-R excludes 0 in favour of proj.
- If TD + L(proj) does not beat the MMFF control, the learned head adds nothing over a force field; say so.
- Single training seed and reduced budget: any result is exploratory.

## Addendum A (pre-declared 2026-10-31, after the T1 E1 results, before any xTB number): xTB control + read-out fix
Motivation: GEOM-QM9 references are GFN2-xTB geometries, so the head may be an amortised RDKit → xTB local relaxation.
The control that tests this is the same relaxation done explicitly.

- **TD + xTB (non-learned control):** each frozen-TD E1 conformer (seeds 0 and 1) optimised with GFN2-xTB 6.7.1
  (`--opt normal`, formal charge, one thread), every W dihedral harmonically constrained to TD's value (force constant
  1.0 Eh), then the W dihedrals restored exactly (as the MMFF control). Non-converged optimisations keep their last
  geometry and are counted; a conformer whose run raises stays TD's and is counted.
- **Read-out fix (deviation from the original run, reported as such):** the 12 molecules with no rotatable bond
  crashed TD's bond head. Fix: for such a molecule the TD forward is run with one bond marked rotatable, its torsion
  output discarded; the atom layers do not read `edge_mask`, so Δ̂ is the model's own (test `test_rigid_molecule_geometry_head`).
  Both T1 checkpoints are re-read out unchanged ("L(proj)fix", "L(noproj)fix"). The original T1 numbers stay the primary
  T1 result; the fixed ones are secondary.
- **Primary endpoint of this addendum:** MAT-R (TD + L(proj)fix) − (TD + xTB), seed 0, paired bootstrap 95 % CI.
  Secondary: the other E1 metrics, seed 1, xTB − TD, CPU time per conformer of each.
- **Rules:** "learned head beats xTB" only if the CI excludes 0 in its favour and seed 1 agrees. If xTB is better or
  indistinguishable, the head is reported as an approximation of constrained xTB relaxation, and its value as speed
  (CPU s / conformer ratio), not accuracy.
- Compute: CPU. The plafnet2 GPU-minute quota is still exhausted (2026-10-31 probe job 49881 held with
  `AssocGrpGRESMinutes`); xtb 6.7.1 has no GPU path anyway.
