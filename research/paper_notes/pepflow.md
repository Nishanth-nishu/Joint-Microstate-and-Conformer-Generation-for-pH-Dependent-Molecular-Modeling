# PepFlow
## Citation
Li J. et al., "Full-Atom Peptide Design based on Multi-modal Flow Matching", ICML 2024 (PMLR 235), arXiv:2406.00735.
Checked 2026-10-26.
## Fulltext status
abstract_only.
## Method (as reported)
Multi-modal flow matching over three geometries at once: backbone frames on SE(3), side-chain torsions on
high-dimensional TORI, and residue identities as categorical distributions on the probability simplex. Partial
sampling of a subset of modalities gives fixed-backbone sequence design and side-chain packing from one model.
Limitations noted by the authors: ODE-based generation limits diversity; no property guidance.
## Relevance to p(z,x|G,h)
The only found model that combines a TORUS-valued continuous modality (torsions - the same manifold as Torsional
Diffusion) with a simplex/categorical discrete modality in one flow. This is the nearest existing mechanism to
p(z, tau | G, h), with residue identity in the role of z. Its "partial sampling" trick is directly the ablation
mechanism for sequential-vs-joint comparison.
