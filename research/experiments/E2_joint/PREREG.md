# E2/E3: conformer-coupled state free energies (joint prototype) vs the sequential baseline. Pre-registration

Written 2026-10-28 ~18:00 UTC (23:30 IST). This is before any conformer feature, fit or E2 metric exists, and while
E1 (job 45306) is still running. E1 fits nothing. The released Uni-pKa zero-shot tautomer MAE (3.21, T pilot) was
already known.

## What is being tested
The question is whether a state distribution that depends on the state-conditioned conformer ensemble predicts
measured tautomer ratios better than (a) the sequential baseline with the same calibration and (b) a capacity-matched
control that uses 2D information only.

## Model (joint-generation prototype; not a validated equilibrium model)
- **State free energy:**
  G_θ(z) = a·G_U(z) + θᵀ φ(X_z), with X_z ~ p(x | z).
  - G_U is the released Uni-pKa free energy in kT (fold mean, E1 protocol).
  - p(x | z) is the state-conditioned conformer generator: ETKDGv3, 30 conformers, per-conformer seeds from
    `conformer_seeds(0, 30)`, then MMFF94s minimisation (≤2,000 iterations). TD is not used here, because most
    Tautobase states fall outside TD-qm9's vocabulary and size.
- **Joint distribution:**
  p_θ(z, x | G, h) ∝ exp(−G_θ(z) − q_z ln10 (pH − PKA_MEAN)) · w_z(x), where w_z are the Boltzmann weights of the
  generated conformers within state z.
  - This is a conditional sampler with a coupled state marginal: θ ≠ 0 means p(z) depends on the conformer ensemble.
  - It is not a normalised density over continuous x; MMFF94s energies are not comparable across states and enter
    only through within-state quantities.
- **φ(X_z)**, all within-state, from validated, minimised conformers with Boltzmann weights at 298.15 K on
  MMFF94s relative energies:
  - φ1 = −ln N_eff, where N_eff = exp(entropy of the weights). A conformational-entropy proxy, dimensionless.
  - φ2 = Boltzmann-weighted count of intramolecular H-bonds: polar H (on N/O) to an N/O acceptor at 1.5–2.5 Å,
    acceptor at least 3 bonds from the donor heavy atom.
  - φ3 = Boltzmann-weighted polar SASA (N, O and polar H; RDKit FreeSASA, probe 1.4 Å), in units of 100 Å².
- **Tautomer task:** log10 K_pred = (G_θ(T1) − G_θ(T2)) / ln10, so swapping T1 and T2 negates the prediction.
  There is no intercept.

## Arms (identical eligible pairs and identical folds)
- **S0:** sequential baseline, released, no fit: a = 1, θ = 0.
- **S1:** sequential + calibration. Fit a only.
- **J1:** joint prototype. Fit a and θ (3 coupling parameters).
- **C1:** capacity control. Fit a and ψ (3 parameters), with ψ from the 2D graph only:
  - ψ1 = number of rotatable bonds;
  - ψ2 = number of topological donor–acceptor pairs 4–7 bonds apart (IMHB *potential*);
  - ψ3 = TPSA (2D) / 100 Å².

  C1 has the same parameter count as J1 and the same kind of information, minus conformers.

## Data, folds, fitting
- **Data:** measured aqueous Tautobase pairs (canonical_v1, label_status = measured; 221). Estimated and calculated
  pairs are not used.
- **Eligibility:** a pair is eligible if both tautomers have at least 1 valid minimised conformer. The same eligible
  set is used for every arm. Ineligible pairs and their reasons are reported.
- **Folds:** the T pilot's 5 scaffold-grouped folds:
  - Murcko scaffold, or the pair itself if acyclic;
  - `scaffold_folds(groups, 5, seed=0)` from `taut_headfit.py`;
  - recomputed on the eligible set.
- **Fitting:** ordinary least squares on log K with ridge λ = 1e-3 on the coefficients (features standardised on the
  training fold), fitted on 4 folds and predicted on the 5th. No hyperparameter is tuned, and no validation split is
  consumed.

## Endpoints and decision rule (fixed now)
- **Primary:** Δ_S = MAE(S1) − MAE(J1) on out-of-fold measured pairs, with a paired scaffold-cluster bootstrap 95% CI
  (10,000 resamples).
- **Co-primary:** Δ_C = MAE(C1) − MAE(J1), same bootstrap.
- **Verdict:**
  - **Coupling benefit SUPPORTED** only if both CIs have lower bound > 0.
  - **NOT SUPPORTED** if either CI includes 0 or is negative. The failing condition is reported.
  - If Δ_S > 0 but Δ_C is not, the gain is attributed to extra parameters or 2D information, not conformer coupling.
- **Secondary (descriptive):**
  - MAE, RMSE and the slope of predicted on measured per arm;
  - fitted coefficients per fold (sign stability);
  - S0 vs S1 (calibration alone).
- **E3 coupling ablation:** J1 with θ forced to 0 is S1 by construction. Additionally, J1 with each φ_k dropped, for
  attribution only.
- **pKa side effect (secondary):** apply G = G_U + θ·φ with the across-fold mean θ, without a, to the 356 external
  transitions. Report the MAE change against E1 oracle. This needs conformers for those states, and runs only if
  compute allows; otherwise it is reported as not run.

## Expected outcome (stated before running)
- S1 improves strongly over S0, because T1 showed compression (slope 0.10).
- The coupling term may help on IMHB-driven tautomers (e.g. enol/keto with a 6-ring H-bond).
- With about 200 pairs and 3 extra parameters, I expect Δ_S to be small, around 60% that its CI includes 0, and Δ_C
  to be smaller still.

## Limits declared now
- MMFF94s is a vacuum force field; there is no solvation model, and φ3 is only a crude solvent-exposure proxy.
- The conformer generator is ETKDGv3+MMFF, not TD, so this tests coupling to a classical ensemble, not to the
  learned generator.
- Tautobase has no conformer-resolved labels, so nothing here validates p(x | z).

## Amendment log
- **2026-10-28, before any feature or fit existed.** "Features standardised on the training fold" is implemented as
  division by the training-fold root-mean-square, with no centring. Centring would add an offset and break the
  antisymmetry required above, since swapping T1 and T2 must negate the prediction. The ridge penalty λ = 1e-3
  applies to the RMS-scaled coefficients.
- **2026-10-28 ~21:50 UTC (2026-10-29 03:20 IST): bug fix, before any fit or E2 metric existed.**
  - **Bug:** features job 45307 gave φ for only 246 of 402 states. `AllChem.MMFFGetMoleculeProperties` rewrites the
    bond types of the Mol it is given to MMFF aromaticity (job 45568). For example, the 2-pyridone ring becomes Kekulé
    single/double bonds. `features.py` called it on the same Mol it then validated, and built the bounds matrix from it.
    As a result, 102 states lost every conformer to `rt_bonds`; 2-pyridone, for example, had 0/30 valid.
  - **Fix:** MMFF now types and minimises a copy. Validation and φ use the RDKit-perceived state with the minimised
    coordinates, and the bounds come from the unmodified state. 2-pyridone now gives 30/30 valid conformers (job 45570).
  - **Features recomputed with the fix; no definition, threshold, arm or decision rule changed.**
  - **Declared limitation (found, not fixed):** the validator's 1–3 tolerance (±0.20 Å from UFF-based DG bounds)
    rejects MMFF94s minima whose angle closes onto an intramolecular contact. Example: 2-hydroxypyridine with the OH
    syn to the ring N has C–O–H ≈ 99°, so C···H is 1.79 Å against bounds of 2.05–2.13 Å; 26 of 30 conformers are
    dropped and only the anti rotamers survive. φ is therefore computed over validated conformers only, which can bias
    φ2 (intramolecular H-bonds) low. Thresholds are kept.
