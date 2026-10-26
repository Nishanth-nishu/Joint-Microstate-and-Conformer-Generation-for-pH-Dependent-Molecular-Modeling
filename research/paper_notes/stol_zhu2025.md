# StoL (Zhu et al. 2025)

## Citation
Chemistry-Enhanced Diffusion-Based Framework for Small-to-Large Molecular Conformation Generation.
arXiv:2511.12182 (2025). Code: https://github.com/Yifei-Zhu/StoL. Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v1 HTML; last ~14k chars not read).

## Method
SMILES -> fragmentation at single bonds (adjacent fragments share >= 3 heavy atoms) -> a SchNet-backbone diffusion
model trained ONLY on small molecules generates ~50 candidates per fragment, clustered to ~10 -> sequential assembly
when shared-atom RMSD < 0.2 Å -> PoseBusters-style validity checks. A planarity-aware loss improves aromatic rings.

## Inputs/outputs — protonation / charge / H handling
"Hydrogen atoms are explicitly included, with element-wise loss weights adjusted by phase." Fragments are restricted
to C, N, O, F. Charge and protonation are not addressed anywhere in the text read => a charged protomer (e.g. a
protonated amine, a carboxylate) is outside the stated element/feature scope. No pH, no state variable.

## Data
Training: 438,142 QCDGE molecules with 6-10 heavy atoms (B3LYP/6-31G*/BJD3). Evaluation: StoL25-init, 200 drug-like
ChEMBL molecules with 16-25 heavy atoms, held out from training (a genuine size-extrapolation split, though not a
scaffold split). DFT comparison via Boltzmann-weighted RMSD; 101/200 molecules below 1.0 Å.

## Relevance to p(z, x | G, h)
Shows that a fragment-assembly decomposition transfers to larger molecules without large-molecule training data.
Not usable as-is for charged species (C/N/O/F neutral fragments only), but the fragment-wise factorization is
conceptually close to a per-site decomposition that a multi-site protonation model would want.
