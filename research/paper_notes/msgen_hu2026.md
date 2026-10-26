# MSGEN / Hierarchical Multi-Scale Molecular Conformer Generation (Hu et al. 2026)

## Citation
Hierarchical Multi-Scale Molecular Conformer Generation. Jiapeng Hu, Weizhi Gao, Zhichao Hou, Xiaorui Liu
(North Carolina State University). ICLR 2026. https://iclr.cc/virtual/2026/poster/10006885. Checked 2026-10-26.

## Fulltext status
abstract_only. OpenReview PDF returned HTTP 403 (bot check) to both curl and WebFetch; the abstract was read from
the ICLR virtual page. NO methodological detail beyond the abstract is asserted below.

## Problem (abstract)
Existing deep generative conformer models ignore the hierarchical structure of molecules and "often fail to preserve
essential substructures, such as ring systems, leading to chemically implausible conformers".

## Method (abstract level only)
Multi-stage generation: coarse key substructures (e.g. scaffolds) are generated first and then used as conditional
guidance for finer scales; a "molecular upsampling" step aligns structural scales between stages. The authors state
the framework integrates with a range of existing molecular generative models.

## Inputs/outputs — protonation / charge / H handling
Unknown (abstract only). What would resolve it: the OpenReview PDF's featurization/appendix section, or released code.

## Relevance to p(z, x | G, h)
The only 2026 conformer generator found with an explicit coarse->fine conditional cascade. The coarse variable is a
substructure/scaffold, i.e. geometric-topological, not a chemical state; but the "coarse stage conditions the fine
stage" wrapper is the structural slot into which a state variable z could be placed in a sequential baseline.
