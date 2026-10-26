# family_D summary - discrete+continuous joint generative models
Scouted 2026-10-26 (UTC). GOOGLE SCHOLAR WAS NOT SEARCHED (no tool available).

## Answer to the two scoped questions
(a) Coupling. Four distinct mechanisms exist, all using ONE shared denoiser over all modalities; they differ in how
the discrete part is noised:
  1. Gaussian noise on relaxed one-hots, same diffusion as coordinates - EDM (QM9 node features are atom type
     one-hot AND an integer atomic charge; one-hots are scaled before diffusion). Inherited by GeoLDM/DiffSBDD.
     FlowMol v1 ablated the sophisticated alternative (SimplexFlow) and reported the naive continuous treatment
     performed BETTER - a relevant negative result.
  2. True categorical (D3PM) transition matrices for discrete + zero-CoM Gaussian for coordinates, one network,
     loss = MSE + cross-entropies, with PER-COMPONENT ADAPTIVE NOISE SCHEDULES - MiDi (exponents nu_r=2.5, nu_y=1.5,
     nu_x=nu_c=1 on QM9, so coarse geometry is fixed first and atom types/charges are refined last). EQGAT-diff is
     the same family with x0-parameterisation, SNR-truncated loss weighting and PubChem3D pretraining. JODO uses one
     Diffusion Graph Transformer with relational attention and supports graph-only generation.
  3. CTMC discrete flow matching for the categoricals + ODE/OT flow for coordinates, shared denoiser, weighted
     modality-specific losses - SemlaFlow, FlowMol3, FlexiFlow. Theory: Campbell et al. 2024 (Multiflow, CTMC
     discrete flows, discrete diffusion as a special case) and Gat et al. 2024.
  4. Diffuse the DISCRETE variable only and read the continuous one off a frozen predictor at each step -
     ProteinGenerator (sequence-space diffusion with RoseTTAFold2 predicting structure each step). PepFlow is the
     outlier geometry-wise: SE(3) frames + TORUS-valued side-chain torsions + simplex-valued residue identity in one
     multi-modal flow, with partial sampling of modalities.
Independent noise schedules are the norm, not the exception: MiDi, MolDiff (bonds noised harder than atoms/
positions) and the loss weights in SemlaFlow all encode an ordering between modalities.

(b) pH / protonation. NONE of them conditions on pH or uses protonation or tautomer state as the discrete variable.
Verified by full-text reads of MiDi, SemlaFlow and FlexiFlow (not just abstract absence), plus abstracts elsewhere.

## Formal charge / H-count as a generated categorical - exactly which models
GENERATE formal charge as a sampled discrete variable: MiDi (atom types, formal charges, bond types under
categorical diffusion), EQGAT-diff (elements, formal charges, bonds, x0-parameterised), JODO (atom types, formal
charges, bonds, coordinates), FlowMol3 (atom types, charges, bond orders via CTMC), FlexiFlow (one-hot charge type
sampled by a discrete-flow update). EDM generates an integer atomic charge as a relaxed continuous node feature on
QM9. SemlaFlow is the informative exception: formal charges are part of the learned joint distribution and get a
cross-entropy loss, but "these do not participate in the generative flow matching process" - they are an auxiliary
head, not a sampled state. No model found generates an H-COUNT or a protonation label.

## Why this matters for p(z,x|G,h)
MolDiff gives the sharpest published argument for joint over sequential in this family: choosing the discrete part
without geometry can leave the geometry unrealisable (atom-bond inconsistency). The same argument transfers to
protonation - a state picked by rules or a pKa model without geometry may have no good conformers. Multiflow and
PepFlow supply the mathematical template (CTMC/simplex for z, continuous flow on the right manifold - including a
torus - for x, one network, partial sampling for conditioning), and MiDi supplies the schedule question (does z get
resolved before or after geometry?). None of them provides the chemistry: pH conditioning and a protonation-aware
z-space are genuinely unoccupied.

## Key unknowns
- Variable atom count under explicit hydrogens: none of the models read here was confirmed to sample the number of
  atoms; adding/removing a proton changes N, so this is an open design problem for z.
- Exactly whether FlowMol3's charges are sampled or auxiliary (as in SemlaFlow) - needs the full text.
- GeoLDM and DiffSBDD charge handling not independently verified; recorded Unknown, not absent.
- No source read here compares a joint discrete+continuous model against its own sequential ablation on a chemistry
  task with a thermodynamic ground truth - which is the evaluation gap this project would fill.
