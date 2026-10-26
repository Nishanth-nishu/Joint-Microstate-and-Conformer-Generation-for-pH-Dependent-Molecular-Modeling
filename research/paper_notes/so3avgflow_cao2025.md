# SO(3)-Averaged Flow Matching + Reflow (Cao, Geiger et al. 2025)

## Citation
Efficient Molecular Conformer Generation with SO(3)-Averaged Flow Matching and Reflow. Zhonglin Cao, Mario Geiger
et al. (NVIDIA). ICML 2025, PMLR 267. arXiv:2507.09785. Checked 2026-10-26. No code URL found in the text read.

## Fulltext status
full_text_read (arXiv v1 HTML).

## Method
Training objective averages the conditional flow over all rotations of the target rather than picking one
(cond-OT) or the Kabsch-optimal one (Fig. 1a). Reflow + distillation straighten trajectories to enable one- and
two-step generation. Two backbones compared: NequIP and a DiT.

## Data / splits
"We follow the exact splitting defined and used in previous works (Ganea et al., 2021; Jing et al., 2022;
Wang et al., 2024). The train/val/test set of GEOM-Drugs contains 243473/30433/1000 molecules" (App. B.1).
Random molecule split, NOT scaffold-separated.

## Inputs/outputs — protonation / charge / H handling
Appendix feature table lists per-atom "charge — Formal charge of atom" and "Total number of hydrogens".
Hydrogens are modelled explicitly (conformer figures colour H atoms). Post-hoc chirality correction by comparing
canonical SMILES and mirroring if mismatched (App.). No pH, no state variable.

## Conformer populations
The paper frames the orbit set with weights: "the integral ... can be written as Σ q̂(x̂) ... where q̂(x̂) is the
weight associated with each conformer" (Sec. 3), but the model is trained on samples and does not output weights.

## Relevance to p(z, x | G, h)
Attractive as the geometry generator if inference cost matters (1-2 steps). Accepts charged inputs via the formal
charge feature. Supplies no z, no pH conditioning, no calibrated populations.
