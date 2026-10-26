# Uni-pKa
## Citation
Zheng H., Luo W., Zhou G., Zhu Z., Yuan Y., Ke G., Wei Z., Gao Z., "Uni-pKa: An Accurate and Physically Consistent
pKa Prediction through Protonation Ensemble Modeling", Research Square 2023, DOI 10.21203/rs.3.rs-3367941/v1
(also ChemRxiv 64e8da3879853bbd786ca4eb). Code: github.com/dptech-corp/Uni-pKa. Checked 2026-10-26.
## Fulltext status
abstract_only (plus secondary descriptions via web search and the AcepKa paper which re-implements it).
Do not treat mechanism details below beyond the stated ensemble formulation as verified.
## Problem
pKa prediction that respects the full protonation (micro)state ensemble rather than regressing a single number.
## Method (mechanism level, as stated in abstract/secondary sources)
A rule-based microstate enumerator builds the protonation ensemble from one input structure; an ML model predicts
the free energy of each microstate; macro-pKa follows from Boltzmann weights over microstates with net-charge and
pH corrections. Discriminative/energy-based, not generative.
## Inputs/outputs
Input: molecular graph (SMILES); pH enters only in the Boltzmann weighting step. Output: micro/macro pKa and
microstate populations. Conformer handling inside Uni-pKa is Unknown from the abstract; AcepKa states accurate
pKa requires reasonable 3D conformer ensembles for the microstates, implying an external conformer step.
## Relevance to p(z,x|G,h)
Defines the z-space (microstates) and the thermodynamic target our p(z|G,h) must reproduce. It is a scoring model
over enumerated z, not a generator, and does not model x jointly. Natural source of ground-truth z weights.
## Unknowns
Whether Uni-pKa itself generates/uses 3D conformers, and which conformer generator; resolvable by reading the
Research Square/ChemRxiv full text or the GitHub repo.
