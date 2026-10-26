# Torsional Diffusion (Jing et al. 2022)

## Citation
Torsional Diffusion for Molecular Conformer Generation. Bowen Jing, Gabriele Corso, Jeffrey Chang,
Regina Barzilay, Tommi Jaakkola. NeurIPS 2022. arXiv:2206.01729. Code: https://github.com/gcorso/torsional-diffusion.
Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v2 HTML) + official code (utils/featurization.py) read.

## Problem
Learn p_G(C) for drug-like molecules faster and better than Cartesian diffusion.

## Method
Factorizes a conformer into local structure L (bond lengths, angles, ring conformations) and torsion angles τ:
"we use RDKit to provide approximate samples from p_G(L), and develop a diffusion-based generative model to learn
distributions p_G(τ | L)" (Sec. 3). Diffusion on the hypertorus T^m; extrinsic-to-intrinsic score model, an
SE(3)-invariant network over the 3D point cloud outputting one scalar per rotatable bond. Parity symmetry imposed.
Training uses "conformer matching": ground-truth conformers are replaced by RDKit-local-structure surrogates (Sec. 4).
Exact likelihoods via the probability-flow ODE (Sec. 3.5), used for energy-based training (Sec. 3.6).

## Data
GEOM-DRUGS and GEOM-QM9; GEOM-XL (102 molecules with >=100 atoms) as a transfer test set.
"The splits are random with train/validation/test of 243473/30433/1000 for GEOM-DRUGS" (Appendix, Splits).
NOT scaffold-separated.

## Inputs/outputs — protonation / charge / H handling
- Input is a single fixed molecular graph from the dataset SMILES; output is torsion angles only.
- Explicit hydrogens: `mol = Chem.AddHs(mol)` in `featurize_mol_from_smiles` (utils/featurization.py); no RemoveHs.
- Formal charge IS an atom feature: `atom_features.extend(one_k_encoding(atom.GetFormalCharge(), [-1, 0, 1]))`
  (utils/featurization.py) — only the bucket {-1, 0, +1} plus "other", so exotic charges collapse to one slot.
- Atom-type vocabulary: `qm9_types = {'H','C','N','O','F'}`, `drugs_types` = 35 elements incl. H.
- Other atom features: atomic number, aromaticity, degree, hybridization, implicit valence, ring membership (3-8),
  ring count. Bond features: single/double/triple/aromatic. Chirality is read into `chiral_tag` but NOT concatenated
  into the node features.
- No pH, no state variable, no alternative graphs. Hydrogen positions come entirely from the RDKit local structure.

## Experiments & baselines
RDKit ETKDG, OMEGA, GeoMol, GeoDiff on AMR/COV recall & precision. Boltzmann-weighted ensemble properties
(E, μ, HOMO-LUMO gap, E_min) computed POST HOC with GFN2-xTB after relaxation — the model does not predict them.

## Torsional Boltzmann generator (Sec. 3.6, 4.5)
Trained with MMFF energies on GEOM-DRUGS molecules with 3-7 rotatable bonds, using the local structures of the
first ground-truth conformer. "A separate model is trained for every temperature" (Appendix). Model density used
"only to reweight the samples to approximate the target density". Annealed training from T'=3000 K downwards.
This is the only part of TD that yields calibrated weights, and it is conditional on fixed L.

## Limitations
Local structure (hence all bond lengths/angles, ring puckers and every X-H geometry) is frozen at RDKit quality;
rings and non-rotatable degrees of freedom are not generated. Random (not scaffold) split.

## Relevance to p(z, x | G, h)
Best-understood candidate for the geometry half of a sequential z→x pipeline, and already reproduced here on
GEOM-QM9. Charged inputs are accepted in principle (formal charge feature, explicit H), but two concrete gaps:
(i) the charge one-hot only covers -1/0/+1, so a di-anion or a doubly protonated base falls in "other";
(ii) a protonation change that alters ring/local geometry or adds an H on a non-rotatable site is invisible to a
model that only moves torsions — the protomer must be fully specified BEFORE the RDKit local structure is built.
