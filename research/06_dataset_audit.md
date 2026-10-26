# 06 — Dataset audit for pH-conditioned joint microstate + conformer generation p(z, x | G, h)
Compiled 2026-10-26 UTC. Sources: family_C and family_E run dirs. Google Scholar was NOT searched
(no tool in this session). Where a fact could not be verified from a page actually fetched, the
cell says `Unknown` plus what would resolve it. Label provenance is tagged
**E** = experimental, **P** = predicted by a model/tool, **C** = computed (QM/MM).

## A. Task definitions used throughout
- **A — valid state–conformer pairs**: generate (z, x) that are chemically valid, i.e. a legal
  protonation microstate graph with explicit H and a physically reasonable 3D geometry for it.
- **B — plausible microstates at pH h**: which microstates are populated at a given pH.
- **C — equilibrium populations of state–conformer pairs**: the actual probabilities p(z, x | G, h).

## B. Dataset table

| Dataset | Molecular identity | State identity | State-specific graph | H representation | Conformers per state | pH | pKa labels | Microstate probabilities | Conformer energies | Solvent | Temperature | Experimental reference populations | License | Size | URL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **GEOM-QM9** | SMILES (QM9, used "as given") | none — one neutral state per molecule | n/a | implicit in SMILES; explicit H in xyz | 13.5 avg, max 1101 | none | none | none | `totalenergy` (Ha), `relativeenergy` (kcal/mol), `boltzmannweight` **C** (GFN2-xTB) | **vacuum** (>98% of GEOM) | 298.15 K (`temperature`) | none | CC0 1.0 | 133,258 species | 10.7910/DVN/JNGTDF |
| **GEOM-Drugs** | SMILES (AICures/MoleculeNet), per-mol `charge` | ionization state fixed by de-salting; MoleculeNet subset "did not perform protonation" | one graph per molecule only | as above | 102.6 avg, max 7451 | none | none | none | same fields **C** | vacuum (ALPB water only for BACE subset) | 298.15 K | none | CC0 1.0 | 304,466 species / 37M confs total | same |
| **QM9** | 133,885 neutral closed-shell CHONF | none | n/a | explicit H | 1 (optimized) | none | none | none | B3LYP/6-31G(2df,p) **C** | gas phase | 0 K / 298 K thermo | none | CC0-type (figshare) | 133,885 | figshare 978904 |
| **QMugs 1.0** | ChEMBL, "Molecules were neutralized" and nonzero-charge species removed | none | n/a | explicit H | 3 (MTD + k-means cluster reps) | none | none | none | GFN2-xTB and ωB97X-D/def2-SVP energies **C** | not stated (no solvation model mentioned) | MTD 300 K; thermo 298.15 K | none | paper CC BY 4.0; data license not stated | 665,911 mols / 1,992,984 confs | 10.3929/ethz-b-000482129 |
| **ChEMBL3D (LoQI)** | ~1.8M ChEMBL compounds | "multiple protonation states and stereochemical configurations" | yes (per the abstract) | Unknown | ~250M geometries total | none | none | none | AIMNet2 energies **C** (implicit solvent) | implicit solvent | Unknown | none | Unknown | 250M geoms | ChemRxiv 689f756e… (separate audit) |
| **OrbNet Denali** | ChEMBL27 | "various protonation states, tautomers, … common salts" | geometry + charge only (XYZ+ comment line) | explicit H in xyz | NMS 300 K / AIMD 500 K samples | none | none | none | GFN1-xTB + ωB97X-D3/def2-TZVP **C** | gas phase | 300 K / 500 K sampling | none | CC-BY-4.0 (ColabFit mirror; not confirmed on figshare) | 2,337,230 geometries | figshare 14883867 |
| **SPICE 2** | PubChem (28,039 mols), dipeptides, monomers | charged and uncharged species present; assignment procedure Unknown | Unknown | explicit H | ~2M conformations total | none | none | none | ωB97M-D3(BJ)/def2-TZVPPD **C** | gas phase | MD-sampled | none | CC0 1.0 (MSKCC catalog) | ~2M confs | Zenodo 10835749 |
| **OMol25** | 83M systems, 83 elements | charge and spin varied; explicit protonation-state sampling **Unknown** | Unknown | explicit H | multiple conformations | none | none | none | ωB97M-V/def2-TZVPD **C** | Unknown (ORCA 6.0.1) | Unknown | none | CC-BY-4.0 | >100M calcs | huggingface.co/facebook/OMol25 |
| **MPCONF196** | 13 peptides/macrocycles | none | n/a | explicit H | 15–16 per compound (196 total) | none | none | none | CCSD(T)/CBS or DLPNO-CCSD(T)/CBS **C** | gas (solvMPCONF196 adds water) | n/a | none | Unknown | 196 conformers | JCTC 2018 |
| **SAMPL6 pKa** | 24 kinase-inhibitor-like fragments | microstate IDs + isomeric SMILES enumerated (e.g. 8 for SM01) | **yes** | explicit charges/H in SMILES | 0 (no 3D distributed) | titration curve (UV-metric, Sirius T3) | macroscopic **E**; 10 microstate-resolved pKa rows of which only **3 from direct NMR**, 7 "Estimated based on SM07/SM14 NMR measurement" | dominant microstate labels only, no populations | none | water (some psKa + Yasuda-Shedlovsky cosolvent extrapolation) | 25 °C | dominant-state identity for ≤8 molecules | MIT (repo) | 24 molecules | github.com/samplchallenges/SAMPL6 |
| **SAMPL7 pKa** | 22 molecules | microstate lists provided | yes | explicit | 0 | titration | macroscopic **E** | none | none | water | 25 °C | none | MIT | 22 | samplchallenges/SAMPL7 |
| **SAMPL8 pKa** | 23 molecules (SAMPL8-1…23) | **278 enumerated microstates** (micro/dortmund/dassault) | **yes** | explicit | **1 OEChem-generated 3D structure per microstate** (.sdf/.mol2/.pdb), no energies | titration (Sirius T3) | macroscopic **E** (23 rows, 1–2 pKa each, with SEM) | none | none | water | 25 °C | none | MIT | 23 mols / 278 states | samplchallenges/SAMPL8 |
| **IUPAC Digitized pKa v2.4a** | SMILES + InChI, 10,626 unique molecules | only `pka_type` (pKa1, pKaH1, pKa2 …) + `acidity_label` A/AH/B — **macroscopic, no microstate graph** | no | n/a | 0 | n/a (aqueous titration) | 24,211 pKa values **E**, with `method`, `assessment` (Reliable/Approximate/Uncertain), `remarks` | none | none | water; `cosolvent` non-blank for 2,786 rows | `T` mostly 25 °C (8,351) and 20 °C (6,524); 3,451 blank/not_stated | none | **CC BY-NC 4.0** | 24,211 rows (6.9 MB CSV) | github.com/IUPAC/Dissociation-Constants |
| **DataWarrior pKaInWater** | ~7,911 structures | strongest acidic / strongest basic site | partial (site, not full microstate) | n/a | 0 | n/a | **E** (but provenance "allegedly literature values", means over duplicates, references lost) | none | none | water | not recorded | none | ships with DataWarrior | ~7,911 | openmolecules.org |
| **Uni-pKa Dwar-iBonD** | DataWarrior + iBonD merge | enumerated into microstates by the Uni-pKa enumerator | **yes** (enumerated, not measured) | explicit | 0 | n/a | 8,232 experimental pKa **E** → 27,138 enumerated microstates (the expansion is **P**) | none | none | water | not recorded | none | Apache-2.0 (code); data via AISSquare | 8,232 → 27,138 | github.com/dptech-corp/Uni-pKa |
| **ChEMBL pKa** | ChEMBL compounds | site-level | partial | n/a | 0 | n/a | mixed; the ~1M-molecule / 3M-microstate ChEMBL pretraining set used by Uni-pKa and the ~2.5M set used by QupKake are **P** (ChemAxon Marvin / internal predictors), not **E** | none | none | water | n/a | none | CC BY-SA 3.0 (ChEMBL) | ~1–2.5M (predicted) | ebi.ac.uk/chembl |
| **Tautobase** | 1,680 unique tautomer pairs | tautomer pair (SMIRKS + structures) | **yes** | explicit | 0 | aqueous | log K / ΔG **E** (some entries estimated or calculated) | pairwise ratios = the closest thing to microstate populations | none | aqueous | not stated in README (Wieder used 298.15 K for QM) | **yes, as pair ratios** | not stated | 1,680 pairs (476 usable in Wieder 2021) | github.com/WahlOya/Tautobase |
| **Neutron-diffraction ligand structures (PDB)** | individual complexes (e.g. 5MOO trypsin + aniline) | protonation visible from D positions **E** | yes, per structure | explicit D | 1 crystal conformer | crystallization pH | none | none | none | crystal | ~100 K / RT | state identity only, no populations | PDB (CC0) | ~88 neutron entries quoted (undated); ligand-bound subset not enumerated here | rcsb.org |

## C. Which dataset supervises which task

- **Task A (valid state–conformer pairs)**: supervised, but only by *stitching*.
  - states: SAMPL8 (278 microstates with one 3D structure each), Uni-pKa's enumerated microstates,
    OrbNet Denali, ChEMBL3D.
  - conformers: GEOM, QMugs, ChEMBL3D.
  - GEOM-QM9 alone supervises **x only**, because QM9 SMILES are neutral and unvaried.
- **Task B (plausible microstates at pH)**: supervised by IUPAC pKa, DataWarrior, Dwar-iBonD,
  SAMPL6/7/8 — all 2D, all macroscopic, all without conformers.
- **Task C (equilibrium populations of state–conformer pairs)**: **not supervised anywhere.**
  The only experimental population-like data are (i) Tautobase pairwise log K (state ratios, no
  conformers, no 3D), and (ii) SAMPL6's dominant-microstate labels for ≤8 molecules, of which only
  3 come from direct NMR. GEOM's Boltzmann weights are computed vacuum conformer weights inside a
  single fixed state, so they are **not** a reference for p(z, x).

## D. Proposed benchmarks

### (a) Structural benchmark — "valid state–conformer generation" — BUILDABLE
Given G (parent molecule, 2D) and a target pH h, the model emits a set of (z, x).
Evaluation, none of which needs populations:
1. **State validity**: z parses, valences and formal charges are consistent, explicit H count
   matches the formal charge, and z is reachable from G by proton add/remove on a titratable atom.
2. **State coverage/precision vs an enumerator reference**: recall and precision of generated z
   against the SAMPL8 (278 states / 23 molecules) and SAMPL6 (68 state files / 24 molecules)
   curated microstate lists, and against Dimorphite-DL / Uni-pKa enumerations on a held-out
   drug-like set. Report both — a generator that emits every tautomer scores high recall and
   useless precision.
3. **Geometry quality per state**: for each generated (z, x), (i) GFN2-xTB single-point and
   relaxation RMSD (the standard "how far from a real minimum" proxy), (ii) coverage/matching
   RMSD against CREST ensembles **regenerated for the charged states** (GEOM gives this only for
   the one state it ships), (iii) sanity on the SAMPL8 microstate 3D structures, which are a free
   one-conformer-per-state reference set.
4. **Negative control**: Torsional Diffusion on GEOM-QM9 cannot do (1)–(2) at all; that is the
   baseline gap this benchmark is designed to expose.
Honest caveat: the "reference" in (2) is an *enumerator*, i.e. **P**, not experiment.

### (b) Population benchmark — "pH-dependent probabilities vs an independent reference"
**Only a weak version can be built honestly today.** There is no experimental dataset of
conformer-resolved microstate populations for drug-like molecules. What can be built, in
decreasing strength of evidence:

- **B1 — macro-pKa recovery (strong, experimental).** From the model's p(z | G, h) compute the
  titration curve, extract macro-pKa, compare to IUPAC v2.4a (24,211 **E** values with
  temperature, method and reliability flags; filter to `assessment == Reliable`, aqueous,
  20–25 °C) and to SAMPL6/7/8 blind sets. Target to beat: RMSE 0.5–0.8 log units (Epik 7,
  QupKake, Uni-pKa). This validates the *marginal over z*, not p(z, x).
- **B2 — tautomer ratio recovery (strong, experimental, and the real test of microstate
  populations).** Tautobase log K, 476 usable aqueous pairs. Current state of the art is poor:
  DFT RMSE ≈3.1 kcal/mol with 13% wrong dominant tautomer; Uni-pKa MAE 3.435 in log K. This is
  where a conformer-aware model could actually win, and it is the single most defensible
  population benchmark available.
- **B3 — dominant-microstate accuracy (weak, experimental).** SAMPL6's NMR-derived dominant-state
  sequence: 3 molecules directly, ≤8 with inference. Too small to train on, usable as a sanity set.
- **B4 — conformer populations within a state (computational reference only).** GEOM
  `boltzmannweight` at 298.15 K, vacuum GFN2-xTB. Report as agreement with a *computational*
  reference and label it **C**; it is not evidence about reality, and it exists only for the one
  protonation state GEOM shipped.
- **B5 — joint p(z, x) (does not exist).** No reference. If the project needs one, it must be
  *manufactured*: e.g. run CREST + CENSO (or an MLIP with implicit solvent) on every enumerated
  microstate of a few hundred molecules and build a free-energy-based p(z, x | h) via the
  Gunner-2020 microstate-free-energy parameterisation. That is a **computed** reference and must
  be reported as such, with its own validation against B1/B2.

**Verdict on (b):** a population benchmark *can* be built, but only as B1+B2 (+B3 as a sanity
check) — i.e. over z marginals and tautomer ratios, from IUPAC/SAMPL/Tautobase. A population
benchmark over the *joint* (z, x) cannot be built from experimental data at all; it can only be
synthesised computationally, and must never be derived from counts of generated or enumerated
states.

## E. Built stratified set (2026-10-26, job 40455; `data_manifests/sensitive_set_v2/manifest.json`)
IUPAC v2.4a filtered to Reliable+Approximate, aqueous (no cosolvent), 20–25 °C, pKa/pKaH types, single fragment,
≤50 heavy atoms, CHNOSPFClBrI: 8,956 rows → 6,395 (molecule, pKa-index) entries, 5,094 molecules (median of repeats;
24 entries with spread >1 unit flagged). Tautobase: aqueous + measured (blank Meas/Calc/Est) + numeric log K → 221 pairs.
Leakage flags vs Uni-pKa fine-tuning (Dwar-iBonD, protonation-invariant InChIKey block): clean pKa entries by
stratum_v2 = sensitive 2,229 / control 678 / other 1,542 (Reliable-only: 74 / 39 / 19); Tautobase clean 35 / 37 / 67.
Scaffold split (Murcko hash 70/15/15). Strata v1 rejected as non-discriminative (75 % sensitive) before any model run.
