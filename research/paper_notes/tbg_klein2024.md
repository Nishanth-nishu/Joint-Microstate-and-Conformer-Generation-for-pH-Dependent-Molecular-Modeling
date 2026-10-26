# Transferable Boltzmann Generators (Klein & Noé 2024)

## Citation
Transferable Boltzmann Generators. Leon Klein, Frank Noé. NeurIPS 2024. arXiv:2406.14426.
Code: https://osf.io/n8vz3/ (view-only link in the paper). Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v2 HTML; last ~3k chars not read).

## Method
Flow matching in Cartesian coordinates producing an invertible map with exact likelihoods, so samples can be
importance-reweighted to the target Boltzmann density. Transferability comes from conditioning on atom type plus
amino-acid identity and sequence position.

## Inputs/outputs — protonation / charge / H handling
- "The atom type embedding a_i is a one-hot vector of 54 classes."
- Amino acid identity (20 classes) and sequence position, "only used for the transferable experiments".
- Hydrogens are explicit particles; hydrogens bound to the same heavy atom are treated as indistinguishable.
- **Charges are not mentioned; protonation states are not discussed.** The target is a fixed Amber-14 topology,
  which fixes the protonation state implicitly in the force field setup.
- Generated samples are validated by bond-graph isomorphism (NetworkX) and chirality: "only generated samples with
  the correct configuration and chirality are considered valid samples" — i.e. a changed chemical graph is treated
  as a failure, not as a second state.

## Data / splits
Train: 200 dipeptides, 50 ns each. Test: 100 unseen dipeptides, 1 µs each (most evaluation on 8 covering all
20 amino acids). Alanine dipeptide separately (Amber ff99SBildn at 300 K, relaxed under GFN2-xTB).
Chemically disjoint train/test systems — the closest thing in this family to a held-out-chemistry split, though it
is a sequence split over dipeptides, not a drug-like scaffold split.

## Conformer populations
Confirmed: importance weights with Kish effective sample size; ESS ≈ 15.3% on 8 test peptides, 13.1% over all 100;
~98% correct-configuration rate.

## Relevance to p(z, x | G, h)
The reference point for "calibrated equilibrium populations, transferable across systems". Two lessons for this
project: (i) exact likelihood + reweighting is what makes populations trustworthy, and the TD torsional Boltzmann
generator is the small-molecule analogue; (ii) current transferable BGs are peptide-scale and assume one fixed
topology/protonation per system, so pH-dependent state populations are entirely outside their scope.
