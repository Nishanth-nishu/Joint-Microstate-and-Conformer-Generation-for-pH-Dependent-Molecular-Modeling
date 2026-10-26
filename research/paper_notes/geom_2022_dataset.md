# GEOM (Scientific Data 2022) — the baseline's data

## Citation
Axelrod S., Gómez-Bombarelli R. "GEOM, energy-annotated molecular conformations for property
prediction and molecular generation." Sci. Data 9:185 (2022). DOI 10.1038/s41597-022-01288-4
(PMC9023519). Data: doi:10.7910/DVN/JNGTDF (Harvard Dataverse, **CC0 1.0**).
Code/tutorials: github.com/learningmatter-mit/geom. Checked 2026-10-26.

## Fulltext status
full_text_read (PMC, first ~100k chars) + README (4.5 kB, downloaded to ext/geom_readme/) +
tutorial notebook 01 + Dataverse JSON API (file list, sizes, license).

## Protonation (the point that matters here)
- QM9 subset: "SMILES strings from the QM9 dataset were used as given" (Methods, SMILES
  pre-processing) → neutral closed-shell only → **no protonation variation**.
- Drugs subset: de-salting "identif[ies] the drug-like compound in each cluster and correct[s] its
  ionization state" (Methods) → charged species exist; per-molecule `charge` field. CREST was run
  with "default CREST arguments, except the molecular charge" (Methods, CREST simulation).
- MoleculeNet subset: the heaviest component was kept but they "did not perform protonation".
- **No pH anywhere; each SMILES appears in exactly one protonation state.**

## Computation
CREST 2.9 with xTB 6.2.3 (Code availability); default E_win 6.0 kcal/mol, E_thr 0.1, RMSD_thr
0.125 Å, B_thr 15.0 MHz. "vacuum conformer-rotamer ensembles (CREs) were generated for over 98% of
the molecules in each dataset" (Background & Summary); ALPB implicit water only for ~99.9% of BACE.

## Schema (from tutorial 01)
Per molecule: `totalconfs`, `uniqueconfs`, `lowestenergy`, `poplowestpct`, `ensembleenergy`,
`ensembleentropy`, `temperature` (298.15 K), `charge`, assay flags.
Per conformer: `xyz`, `geom_id`, `set`, `degeneracy`, `totalenergy` (Ha), `relativeenergy`
(kcal/mol), `boltzmannweight` = d·e^{−E/kT} / Σ d_i e^{−E_i/kT}, `conformerweights`.

## Sizes / files
QM9 133,258 species (13.5 conformers avg, max 1101); Drugs 304,466 species (102.6 avg, max 7451);
37M conformations overall. Dataverse files: rdkit_folder.tar.gz 50.1 GB, drugs_crude 42.7 GB,
crest_hess_and_orbs 30.1 GB, qm9_crude 1.0 GB, qm9_featurized 0.13 GB (not downloaded).

## Relevance to p(z,x|G,h)
GEOM-QM9 supervises x only; z and h are unidentifiable from it. GEOM's Boltzmann weights are
**computed, vacuum, within-one-state** conformer weights — usable as a computational reference for
p(x|z) but never as evidence about p(z|h) or about real solution populations.
