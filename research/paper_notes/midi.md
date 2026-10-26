# MiDi
## Citation
Vignac C., Osman N., Toni L., Frossard P., "MiDi: Mixed Graph and 3D Denoising Diffusion for Molecule Generation",
ECML-PKDD 2023, arXiv:2302.09048. Code: github.com/cvignac/MiDi. PDF: research/papers/arxiv_2302.09048.pdf. Checked 2026-10-26.
## Fulltext status
full_text_read (ar5iv HTML via WebFetch).
## Method (mechanism level)
Joint denoising diffusion over (coordinates, atom types, FORMAL CHARGES, bond types). Forward process is a product
of a zero-CoM Gaussian over positions and independent categorical (D3PM-style) noise over atom types, charges and
bond types (Sec. 4.1); "the 3D coordinates are corrupted with Gaussian noise, while the other components use
discrete diffusion" (Sec. 1). A single rEGNN/graph-transformer denoiser predicts the clean graph from the noisy one;
the loss is "a weighted sum of these components" (Sec. 4.4): MSE on coordinates, cross-entropy on each categorical.
Key design: PER-COMPONENT ADAPTIVE NOISE SCHEDULES - "we modify the noise schedule to vary according to the
component" with exponents nu_r=2.5, nu_y=1.5, nu_x=nu_c=1 on QM9 (Sec. 4.2), so geometry is resolved first and
atom types/charges are refined late.
## Inputs/outputs
Unconditional (or property-conditioned) generation of complete molecules; charges are an output categorical.
## Relevance to p(z,x|G,h)
The single most transferable mechanism for family A: it shows a discrete chemical label that is chemically close to
protonation (FORMAL CHARGE per atom) being denoised jointly with 3D coordinates in one network with coupled but
separately scheduled noise. The schedule ordering (coarse geometry -> then types/charges) is directly a design
question for joint z/x sampling. No protonation/pH/tautomer semantics anywhere in the paper.
## Exact supporting quotes
- Sec. 4: "the type and formal charge associated to all atoms".
- Sec. 2: "MiDi is the only model ... capable of handling formal charges".
- Sec. 4.2: "we modify the noise schedule to vary according to the component".
