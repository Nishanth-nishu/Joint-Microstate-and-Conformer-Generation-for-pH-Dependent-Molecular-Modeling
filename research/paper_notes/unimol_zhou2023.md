# Uni-Mol (Zhou et al. 2023)

## Citation
Uni-Mol: A Universal 3D Molecular Representation Learning Framework. Gengmo Zhou, Zhifeng Gao, Qiankun Ding,
Hang Zheng, Hongteng Xu, Zhewei Wei, Linfeng Zhang, Guolin Ke. ICLR 2023.
Code: https://github.com/deepmodeling/Uni-Mol. Checked 2026-10-26.

## Fulltext status
abstract_only (ICLR listing + ChemRxiv record + repo description). Featurization was NOT resolved in this run.

## Method (abstract level)
Two pretrained SE(3)-Transformers sharing an architecture: a molecular model pretrained on 209M conformations and a
pocket model on 3M pocket candidates. Downstream tasks include molecular property prediction (14/15 SOTA),
protein-ligand binding pose prediction and molecular conformation generation.

## Inputs/outputs — protonation / charge / H handling
Unknown. Whether Uni-Mol's conformer pipeline keeps explicit hydrogens or uses a formal-charge channel was not
established; the search returned no authoritative statement. What would resolve it: `unimol/data/` preprocessing
and the conformation-generation task config in deepmodeling/Uni-Mol.

## Relevance to p(z, x | G, h)
Relevant as a large pretrained 3D encoder that a state-prediction head could sit on, and as a reminder that the
pretraining corpus (RDKit/ETKDG-derived conformations) carries the same single-protomer assumption as GEOM.
