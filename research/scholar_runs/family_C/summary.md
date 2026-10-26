# family_C summary — protonation, thermodynamics, and whether traditional methods already solve this
(2026-10-26 UTC; Google Scholar NOT searched — no tool available.)

## 1. The closest existing work, at mechanism level

**Uni-pKa (JACS Au 2024)** is the nearest neighbour of the proposed model. It enumerates the full
*protonation ensemble* (microstates = explicit protonation/tautomer graphs), predicts one standard
free energy per microstate with the 3D Uni-Mol encoder, and then obtains micro-pKa, macro-pKa and
pH-dependent distribution fractions analytically. That is exactly the discrete half, p(z | G, pH),
of p(z, x | G, h). What it does **not** do: it never models p(x | z). Conformers appear only as
training-time data augmentation ("swiftly generate multiple random conformations allows us to
incorporate data augmentation", Sec. 5.6), with no stated count and no Boltzmann aggregation. Its
derivatives — Starling (Rowan, ChemRxiv 2025) and AcepKa (PlayMolecule, arXiv 2604.00841) — reuse
the same architecture; AcepKa states "Accurate pKa prediction requires reasonable ensembles of 3D
conformers for the microstates" (Sec. 3.2) but gives neither a conformer count for the pKa pipeline
nor an aggregation rule, and reports no pKa RMSE.

**QupKake (JCTC 2024)** predicts micro-pKa per site from a GNN with GFN2-xTB-derived *scalar*
features (Fukui, Wiberg, charges). 3D is feature engineering; no ensemble, no populations.

**Epik 7 (Schrödinger)** is the production tool for task B. It is explicitly 2D: "The topological
nature of the ML approach means that Epik 7 ... is rapid but agnostic to 3D geometry and
stereochemistry." It emits state energies and "the most populated states for a user-supplied pH"
for Macro-pKa, so it answers "which microstates at pH h" without any conformer.

**Constant-pH MD** (Donnini λ-dynamics 2011, Amber/CHARMM/GROMACS CpHMD, phbuilder 2024, Chodera's
`protons` for OpenMM) is the mature *physics* answer to coupled protonation + conformation: it
samples z and x jointly by construction. But the mature benchmarks are protein residues
(RMSE ~0.8–0.9 pKa units); no drug-like small-molecule CpHMD benchmark was found in this session.
`protons` even obtains its ligand states from Epik and then *calibrates* them to reproduce
reference-state populations — i.e. the populations are an input, not an output.

## 2. Does any of this already solve task B (plausible microstates at a pH)?

Largely yes, for routine drug-like chemistry, and without 3D:
- macro-pKa RMSE 0.5–0.8 log units is now standard (QupKake 0.54–0.79; Uni-pKa MAE 0.49–0.81;
  Epik 7 median |err| 0.42, RMSE 0.72; SAMPL6 top entries RMSE 0.7–1.0).
- Epik 7 "recovers 95% of the most populated protonation states" relative to earlier versions.
So a generative model cannot claim novelty on "which states exist at pH h". The residual failure is
*which tautomer within a charge state*: SAMPL6 concluded "Prediction of dominant tautomers was a
major source of error for microscopic pKa predictions, especially errors in charged tautomers", and
Uni-pKa scores an MAE of 3.435 in log K on Tautobase. Tautomer/microstate *ratios* (task C) are
where every method is still weak — DFT RMSE ~3 kcal/mol (~2 pKa units) with 13% wrong dominant
tautomer (Wieder 2021).

## 3. Conformation-dependence of microstate populations — the magnitude

Evidence **for** a real effect:
- Over et al. 2014 (J. Med. Chem.): eight stereoisomers of one inhibitor; the diastereomers that can
  form an intramolecular NH→NR3 H-bond have tertiary-amine pKa 6.08–6.18 vs 7.07–7.16 for those that
  cannot — "a pKa difference of ∼1 unit", measured by Sirius GLpKa titration, with QM ΔpKa = 1.1.
  So a conformational feature (IMHB) is worth ≈1 pKa unit ≈ 1.4 kcal/mol ≈ 10× in state population.
  Caveat: the conformational preference is locked in by stereochemistry, so this is a between-molecule
  comparison, not within-molecule conformer coupling.
- Wieder 2021: moving from RRHO thermochemistry on minima to full configurational sampling changes
  tautomer ΔG by "a mean absolute error of ∼1 kcal mol−1".
- Uni-pKa Sec. 2.5: on one hard system, Jaguar's MAE improves "from MAE = 1.07 down to 0.20" with
  exhaustive conformational modelling; Table 2's "Jaguar++" is defined as "Jaguar with
  comprehensive conformational search, weighting 10 conformers".

Evidence **against** (and this is the honest counterweight):
- Wieder 2021's own ablation: "Much more relevant than including multiple conformations is locating
  the global minimum conformation" — only ~12 of 468 tautomer pairs changed by >1 kcal/mol when
  multiple minima were used instead of the global minimum.
- Uni-pKa reaches state-of-the-art macro-pKa accuracy **with no conformational sampling at all**,
  matching or beating conformer-searching DFT.
Net: for typical drug-like molecules the conformer-ensemble correction to microstate populations
is of order 0.5–1 kcal/mol (≈0.4–0.7 log units), i.e. real but usually smaller than current model
error, and it is dominated by *finding the right minimum* (e.g. the IMHB-forming one) rather than
by weighting many conformers correctly. Conformation-dependence is therefore best motivated as a
tail phenomenon (IMHB-capable, zwitterionic, macrocyclic, conformationally locked systems), not as
a generic claim.

## 4. Open unknowns
- QupKake's actual 3D/conformer protocol (SI and source not read).
- The number of conformers Uni-pKa generates for augmentation, and whether its free-energy head is
  conformer-invariant in practice (would require running the released code).
- Whether any published work measures *within-molecule* conformer-resolved microstate populations
  experimentally. None found. This is the gap that blocks a true population benchmark (see
  research/06_dataset_audit.md).
- Epik 7's peer-reviewed version and its exact population formula (abstract_only here).
