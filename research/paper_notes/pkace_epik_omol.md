# Sequential-pipeline components: pKAce, Epik v7, OMol25
## Citation
Mayol-Llinas J. et al. (PlayMolecule), "pKAce: Small Molecule Protonation through Equivariant Neural Networks",
arXiv:2407.11103, 2024. Johnston R.C. et al., "Epik: pKa and Protonation State Prediction through Machine Learning",
ChemRxiv 2023 / JCTC. Levine D. et al., "The Open Molecules 2025 (OMol25) Dataset, Evaluations, and Models",
arXiv:2505.08762. Checked 2026-10-26.
## Fulltext status
abstract_only for all three.
## What they are
pKAce: equivariant NN predicting protonation sites; trained on several conformers per molecule, reported to reduce
prediction variability across conformations. Epik v7: ML pKa/protonation-state prediction, states enumerated at a
target pH (pH 7.4 +/- 2.0 default). OMol25: large QM dataset/benchmark whose evaluation suite explicitly includes
protonation and conformer optimisation tasks, with protonation reported as the hardest energy-based evaluation.
## Relevance to p(z,x|G,h)
These are the "predict/enumerate the state first" half of the sequential baseline the project must beat, and
together with RDKit/ETKDG or Torsional Diffusion they define that baseline concretely. pKAce's finding that
predictions vary across conformers is direct evidence that z and x are coupled - i.e. motivation for joint
modelling. OMol25 gives a protonation-aware energy benchmark that could score generated (z,x) pairs.
