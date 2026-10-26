# Multiflow / Discrete Flow Models (Campbell et al. 2024)
## Citation
Campbell A., Yim J., Barzilay R., Rainforth T., Jaakkola T., "Generative Flows on Discrete State-Spaces: Enabling
Multimodal Flows with Applications to Protein Co-Design", ICML 2024 (PMLR 235:5453-5512), arXiv:2402.04997.
Code: github.com/jasonkyuyim/multiflow. PDF: research/papers/arxiv_2402.04997.pdf. Checked 2026-10-26.
## Fulltext status
abstract_only.
## Method (as reported)
Discrete Flow Models: the discrete analogue of continuous flow matching realised with Continuous Time Markov
Chains; discrete diffusion appears as a special case. Multimodal generation combines a CTMC over the discrete
modality (amino-acid sequence) with a continuous SE(3) flow over structure in ONE model, trained jointly; the same
model can then generate sequence given structure, structure given sequence, or both.
## Relevance to p(z,x|G,h)
The canonical mathematical template for p(z,x): CTMC for z, continuous flow for x, one network, and - crucially -
the ability to *marginalise or condition* on either modality at inference, which is exactly how one compares joint
generation to the sequential enumerate-then-generate pipeline within a single trained model.
