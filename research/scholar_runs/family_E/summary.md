# family_E summary — what the data actually supervises
(2026-10-26 UTC; Google Scholar NOT searched. Full table: research/06_dataset_audit.md)

## Headline
No public dataset pairs **protonation microstates** with **conformer ensembles** and
**experimentally referenced populations**. The three ingredients exist in three disjoint corpora:

1. **Conformer ensembles with computed Boltzmann weights, one fixed protonation state per
   molecule** — GEOM (QM9 + Drugs), QMugs, ChEMBL3D. GEOM stores per-conformer `totalenergy`,
   `relativeenergy`, `boltzmannweight` (= degeneracy-weighted Boltzmann factor), `degeneracy`, with
   a molecule-level `charge` and `temperature` = 298.15 K. These weights are **computed**
   (GFN2-xTB, CREST 2.9 / xTB 6.2.3, **in vacuum** for >98% of molecules), never measured.
2. **Protonation/tautomer microstate identities with experimental equilibrium constants, no 3D** —
   IUPAC Digitized pKa (24,211 rows, macroscopic pKa1/pKa2 only), DataWarrior pKaInWater (~7.9k),
   Tautobase (1,680 pairs, log K), SAMPL6/7/8.
3. **Geometries of many protonation states / tautomers, no populations** — OrbNet Denali
   (ChEMBL27, 2.3M geometries "including various protonation states, tautomers … and common
   salts"), SAMPL8 microstate structure files, ChEMBL3D.

## The critical GEOM finding for this project
GEOM-QM9 uses "SMILES strings from the QM9 dataset … as given" — QM9 is neutral closed-shell, so
**GEOM-QM9 contains no protonation variation at all**: h (pH) is unidentifiable from it and the
baseline Torsional-Diffusion-on-GEOM-QM9 setting supervises only x, never z. GEOM-Drugs is
slightly better: de-salting "identif[ies] the drug-like compound in each cluster and correct[s] its
ionization state", and CREST was run with "default CREST arguments, except the molecular charge",
so charged species exist and the per-molecule `charge` field records them. But each SMILES appears
in exactly **one** protonation state; there is no molecule for which two microstates of the same
parent both have ensembles, and there is no pH annotation anywhere.

## Hard constraint on any population benchmark
Boltzmann weights in GEOM are conformer weights *within* one fixed state, computed in vacuum at
GFN2-xTB. They cannot be combined across states (different numbers of atoms, no common reference),
and the number of enumerated microstates is never a probability. Any p(z) target must come from
free energies / measured equilibrium constants (Gunner et al. 2020), i.e. from corpus 2.

## Open unknowns
- SPICE 2 / OMol25: do they contain *paired* protonation variants of the same parent molecule?
  Unknown from the pages read; resolving step = read the SPICE generation scripts and the OMol25
  arXiv dataset section.
- Tautobase per-record solvent/temperature columns (the .dwar file was not opened; the DataWarrior
  format needs the DataWarrior reader).
- Whether OrbNet Denali's XYZ+ records can be mapped back to microstate SMILES (the comment line
  carries only charge and multiplicity).
