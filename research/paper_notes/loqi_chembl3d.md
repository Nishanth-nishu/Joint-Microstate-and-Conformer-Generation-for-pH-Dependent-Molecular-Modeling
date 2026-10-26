# LoQI / ChEMBL3D
## Citation
Nikitin F., Anstine D.M., Zubatyuk R., Paliwal S., Isayev O., "Scalable Low-Energy Molecular Conformer Generation
with Quantum Mechanical Accuracy", ChemRxiv 2025, DOI 10.26434/chemrxiv-2025-k4h7v (v3 2026-07);
NeurIPS 2025 listing. Code/data: github.com/isayevlab. Checked 2026-10-26.
## Fulltext status
abstract_only (ChemRxiv abstract + ChEMBL announcement; PDF not text-extractable in this session).
## Problem
Scalable generation of low-energy 3D conformers with near-QM accuracy.
## Method (as stated in abstract)
Diffusion model trained on ChEMBL3D (>250M geometries, 1.8M drug-like compounds) optimized with AIMNet2 NNPs
in implicit solvent. Stereochemistry is handled by graph augmentation, enabling generation targeted at a specified
stereo configuration.
## Inputs/outputs
Input: molecular graph (stereo-augmented). Output: 3D conformers. pH is NOT an input (not mentioned).
## Relevance to p(z,x|G,h)
Closest "conditioning on a discrete chemical descriptor" analogue in conformer generation: stereochemistry is
injected as a discrete graph label the model is conditioned on, which is structurally the same trick one would use
to condition on a protonation label z. Its dataset reportedly "captures complex organic molecules in various
protonation states", i.e. protonated species are in the data distribution, but the model CONDITIONS on the given
(already protonated) graph rather than generating z. No pH conditioning.
## Unknowns
Exactly how protonation states appear in ChEMBL3D and whether any state label is exposed to the model; resolvable
from the ChemRxiv full text / dataset card.
