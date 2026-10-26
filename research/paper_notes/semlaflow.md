# SemlaFlow
## Citation
Irwin R. et al., "SemlaFlow - Efficient 3D Molecular Generation with Latent Attention and Equivariant Flow Matching",
AISTATS 2025 (PMLR 258), arXiv:2406.07266. Code: github.com/rssrwn/semla-flow.
PDF: research/papers/arxiv_2406.07266.pdf. Checked 2026-10-26.
## Fulltext status
full_text_read (arXiv HTML v3 via WebFetch).
## Method (mechanism level)
Joint distribution over z=(x coordinates, a atom types, b bond orders, c formal charges) (Sec. 4). ONE shared Semla
denoiser emits multiple vector fields: "we parameterise a single Semla model to generate multiple vector fields".
Coordinates: equivariant OT flow matching, x_t ~ N(t x_1 + (1-t) x_0, sigma^2), sigma=0.2, t ~ Beta(2,1).
Atom and bond types: discrete flow models (Campbell et al. 2024) with categorical interpolants mixing data and a
uniform distribution (Eqs. 15-16). NOTABLE: formal charges are PREDICTED but "do not participate in the generative
flow matching process" - i.e. charge is an auxiliary head, not a sampled discrete state. Loss (Eq. 17) = MSE on
coordinates + cross-entropy on types/bonds/charges, with weights (1.0, 0.2, 0.5, 1.0) on QM9.
## Relevance to p(z,x|G,h)
Demonstrates shared-denoiser coupling with modality-specific interpolants (continuous ODE + discrete CTMC) and an
explicit design choice about which discrete variables are *sampled* vs *read off*. Protonation z would need to be
in the sampled set, not the auxiliary-head set. No protonation/pH/tautomer content.
