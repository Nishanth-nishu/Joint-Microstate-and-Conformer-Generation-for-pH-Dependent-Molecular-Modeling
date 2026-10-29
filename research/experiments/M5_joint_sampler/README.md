# M5: joint-generation prototype sampler (demo on the 6 smoke parents, pH 7.4)

Code: `jmcg/joint/sample.py`. Job 45584 (2026-10-28 21:58 UTC / 2026-10-29 03:28 IST), seed 0, 20 (z, x) pairs per
parent. Arms:
- S0: released Uni-pKa;
- S1: calibrated a;
- J1: a + θ·φ, final coefficients from E2 `fit/e2_result.json`.

## Objective (what is sampled)
p_θ(z, x_k | G, h) = p(L | h) · p(z | L) · w_z(x_k), where:
- **Level free energy:** F_L = −LSE_{z∈L}(−G_U(z)). This is released Uni-pKa, so macro-pKa values are unchanged by
  construction.
- **Level probability:** p(L | h) = softmax_L(−F_L − L·ln10·(pH − PKA_MEAN)).
- **Within-level model (E2):** p(z | L) = softmax_{z∈L}(−(a·G_U(z) + ln10·θ·φ(X_z))).
  - Used only if every state of L has a validated ensemble; otherwise θ = 0 for that level, and the level is listed in
    `run_info.json` under `fallback_levels_theta0`.
- **Conformer:** x_k is one of the state's validated ETKDGv3+MMFF94s conformers, drawn with its Boltzmann weight w_z at
  298.15 K on MMFF94s energies within the state.

This is a conditional sampler with a coupled state marginal. It is not a normalised density over continuous x, and
MMFF energies are never compared across states.

## Results
- **S0 identity:** the S0 state probabilities equal softmax(fold-mean G_U) to 1.1e-16.
  - Against the sequential state stage they differ by up to 0.014, because the stage averages the per-fold
    populations while the sampler averages G first. The averaging differs; there is no error.
- **J1 vs S1:** within the levels where the coupling applies (acac, phenol, propylamine, histamine), J1 and S1 differ
  by < 1e-4. Example: histamine amine-cation tautomer 0.88769 vs 0.88764. In this demo, then, the coupling changes
  nothing visible.
- **Fallback (θ = 0) levels:** glycine 0 (zwitterion and neutral; mass 0.996), glycine −1, ibuprofen 0 and −1, and
  phenol −1.
  - Their states have no validated MMFF94s ensemble; all failures are `bonds_angles`.
  - So the dominant pH 7.4 states of glycine and ibuprofen have no conformer under this generator. 40 of 120 sampled
    pairs have no x.
- **Failure found and fixed in this run:** the first version masked out states without φ and renormalised within the
  level. That moved glycine's neutral level entirely onto the non-zwitterionic form, giving neutral glycine
  p = 0.99999 at pH 7.4, which is wrong. The θ = 0 level fallback replaced it.

## Status
- The prototype runs end to end and reduces exactly to the sequential baseline at θ = 0.
- Its conformer generator (ETKDGv3+MMFF94s + validator) fails on zwitterions and carboxylates/phenolates. Under the
  current thresholds that blocks the coupled model on the most common physiological charge states.
- Unblocking requires either validator bounds derived from the force field used, or an unminimised generator (TD or
  ETKDG) for φ. Either is a design change for the user to decide on; it is not made here.
