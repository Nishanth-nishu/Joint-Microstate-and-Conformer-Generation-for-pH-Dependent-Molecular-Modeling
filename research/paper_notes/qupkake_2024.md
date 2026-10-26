# QupKake (JCTC 2024)

## Citation
Abarbanel O. D., Hutchison G. R. "QupKake: Integrating Machine Learning and Quantum Chemistry for
Micro-pKa Predictions." J. Chem. Theory Comput. 2024. DOI 10.1021/acs.jctc.4c00328 (PMC11325546).
Code: https://github.com/hutchisonlab/QupKake. Checked 2026-10-26.

## Fulltext status
full_text_read (PMC, first ~100k chars).

## Problem / Method
Three stages: (1) GFN2-xTB tautomer search with implicit water; (2) GNN+QM enumeration of
reaction sites; (3) GNN regression of micro-pKa per site. Node/edge/graph features "generated using
RDKit and GFN2-xTB" (Sec. 2.2.1) — coordination number, Fukui(0), Wiberg bond order, partial
charge, polarizability, xTB energy. Output: "several linear layers that outputs the pKa of that
site" (Sec. 2.3).

## 3D / conformers
3D geometry enters only as xTB-derived scalar atom features; the paper does not state how many
conformers are generated or whether the lowest-energy one is kept, and conformer dependence of
micro-pKa is **not discussed** (Unknown — resolving this needs the SI / source code).
Tautomers are enumerated with "the tautomer enumeration function of version 2022.03.3 of the RDKit
software package" (Sec. 2.1).

## Data
Training: ~2.5M ChemAxon-Marvin-**predicted** pKa values from ChEMBL ("consisting of ∼2.5 million
predicted acidic and basic pKa values", Sec. 2.3); finetuned on "the 5637 experimental micro-pKa
values" (Sec. 2.3).

## Outputs
A single micro-pKa per site. **No** macro-pKa, **no** pH-dependent populations, **no** microstate
free energies.

## Experiments
External sets: Novartis (280), Literature, SAMPL6/7/8. "RMSEs of 0.79 and 0.54 pKa units,
respectively" (Novartis, Literature; Sec. 3.2.1); SAMPL6 0.44, SAMPL7 1.04, SAMPL8 0.85 (Table 1).

## Relevance to p(z,x|G,h)
A strong site-level micro-pKa baseline for task B, but it never produces a joint state–conformer
object or a population; its "3D" is feature engineering, not a conformer ensemble.
