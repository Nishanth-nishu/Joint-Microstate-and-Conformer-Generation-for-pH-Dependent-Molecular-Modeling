# EC-Conf (Fan et al. 2024)

## Citation
EC-Conf: An ultra-fast diffusion model for molecular conformation generation with equivariant consistency.
Journal of Cheminformatics 16 (2024). DOI 10.1186/s13321-024-00893-2. arXiv:2308.00237.
Code: https://github.com/DeepLearningPS/EcConf. Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v6 HTML).

## Method
Consistency model over the probability-flow ODE with a modified Equiformer backbone: "a modified SE(3)-equivariant
transformer model was directly used to encode the Cartesian molecular conformations". Time and atom-type embeddings
are linearly projected and summed as input. EMA target network; 1-step sampling with optional refinement steps.

## Inputs/outputs — protonation / charge / H handling
The paper describes a 2D graph with node and edge features but does NOT specify the atom feature list. Whether
formal charge or explicit hydrogens are used is Unknown. What would resolve it: the dataset/featurization code in
DeepLearningPS/EcConf. The network input is explicitly only time + atom-type embeddings, which suggests element
identity alone reaches the score model — but this is not asserted here.

## Data / splits
ConfGF-style splits: GEOM-QM9 39860/4979/200 molecules; GEOM-Drugs 39852/4983/200 molecules. NOT scaffold-separated.

## Conformer populations
Framed as learning "the Boltzmann distribution of C under condition of G", but trained on low-energy samples only,
with no energy term and no weight output.

## Relevance to p(z, x | G, h)
Interesting only as a fast sampler. Its conditioning information may be too thin (atom type + time) to distinguish
protomers of the same heavy-atom skeleton; that is the first thing to check before considering it.
