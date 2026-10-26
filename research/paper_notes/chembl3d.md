# ChEMBL3D (dataset)

## Citation
- Title: ChEMBL3D: Quantum-Accurate 3D Conformers for ChEMBL at Scale (KiltHub title); readme.txt title "ChEMBL3D: A Large-Scale Collection of QM-Accurate Conformers for ChEMBL Molecules"
- Authors (dataset README): Filipp Nikitin, Roman Zubatyuk, Olexandr Isayev
- Year / venue: 2026 (KiltHub published 2026-03-10T14:28Z; README "posted on March 1, 2026"); described in the LoQI preprint (ChemRxiv 2025, doi 10.26434/chemrxiv-2025-k4h7v)
- DOI: 10.1184/R1/31428449 — https://kilthub.cmu.edu/articles/dataset/_b_ChEMBL3D_Quantum-Accurate_3D_Conformers_for_ChEMBL_at_Scale_b_/31428449
- License: MIT
- Code: loader `scripts/sgdataset.py` in the release; consumer code https://github.com/isayevlab/LoQI (`data_processing/process_chembl3d.py`)
- Date checked: 2026-10-26 08:26 UTC (13:56 IST)

## Fulltext status
Dataset documentation read in full: README.md (6172 B) and readme.txt (7010 B) from figshare files 62277919 / 62555893, plus loader scripts, saved to
`/home2/nishanth.r/jmcg_stage/ext/chembl3d_meta/chembl3d/`. The multi-GB archives (topologies 1.6 GB; zarr 11 parts, ~161 GB) were NOT downloaded. No topology SDF was inspected, so per-molecule protonation variants were not observed directly.
LoQI paper PDF not accessible (ChemRxiv Cloudflare/403; OpenReview verification page) — paper-level statements about ChEMBL3D come only from the abstract.

## Problem
Provide a very large, QM-quality (AIMNet2, implicit solvent) conformer dataset for ChEMBL v34 drug-like molecules, for ML conformer generation / ranking.

## Method (curation pipeline, README "Dataset Description")
1. Protomer generation with OpenEye **FixpKa**, "merge with original structures".
2. Stereocenter enumeration for unannotated centers (OpenEye Flipper).
3. Hydrogen-complete topology generation (Open Babel).
4. Conformer generation (OpenEye OMEGA Classic).
5. FIRE geometry optimisation with AIMNet2 potential ("toward DFT/CPCM quality").
6. Topology validation, filtering, deduplication; conformers with relative_energy > 6 kcal/mol removed.
Statistics: ~505M initial → ~280M processed conformers; ~2.4M unique structures for ~1.8M unique ChEMBL IDs; zarr holds 283,577,147 conformers.

pH: the README does not state a pH. OpenEye's FixpKa documentation says pKa normalisation sets "its most energetically favorable ionization state for pH=7.4" (docs.eyesopen.com/applications/quacpac/fixpka/fixpka_overview.html). So the protomer step is most likely a *single* rule-based pH-7.4 state per input, merged with the ChEMBL-deposited form — i.e. typically ≤2 protonation variants per compound. This is an inference from tool documentation, not a statement by the dataset authors (status: Partially supported).

## Data / schema (README "Zarr Format (Per Group)")
Per atom-count group NNN (007–142): `numbers (N,A) int8`, `coord (N,A,3) float32`, `energy (N,) float64` (absolute), `relative_energy (N,)` (relative to the lowest conformer in the `(mol_id, stereo_id)` block), `stereo_id (N,) int32`, `mol_id (N,) |S15` (e.g. `CHEMBL100101_0`), `charge (N,) int8` ("molecular charge value stored by the pipeline"). Topologies: `topologies/NNN.sdf` (+ .smi), one H-complete SDF record per `mol_id`; SDF `_Name` == zarr `mol_id`.
- Protonation state lives in the *topology* (explicit H count + formal charges in the SDF record) and the total `charge` array; there is no field for pH, pKa, microstate free energy, or state population/probability.
- Variant suffix `_k` on `mol_id`: LoQI's processor treats `CHEMBL123_0`, `CHEMBL123_1`… as variants of one parent (process_chembl3d.py:130-187). Stereoisomers are *not* the suffix (they are `stereo_id` within one `mol_id`), so suffix variants are distinct topologies — consistent with protomer (or other topology) variants. What exactly each suffix is (protomer vs original vs tautomer) is **Unknown**; resolve by inspecting `topologies/*.sdf` for one parent with ≥2 suffixes.
- Energies: absolute AIMNet2 energies are comparable across conformers of the same topology; comparing absolute energies across protonation states with different numbers of H/charge is not a free-energy of protonation (no proton chemical potential / pH term). The dataset does not offer such a comparison.

## Inputs/outputs
Not a model. Supplies: H-complete topology (with formal charges) + optimised conformers + energies per (mol_id, stereo_id).

## Evidence of charged states in LoQI's processed training split
From `train_charges_h.npy` (LoQI release 10.1184/R1/31441570; 944 B) parsed locally: per-element formal-charge distribution over charge bins {-2,-1,0,+1,+2,+3}. N: P(+1)=0.0816, P(-1)=0.0017; O: P(-1)=0.0410; S: P(-1)=0.0051, P(+1)=0.0063. Overall atom prior (`train_charges_prior_h.npy`): -1: 0.25%, 0: 99.20%, +1: 0.55%. Confirms charged (protonated amine / carboxylate-like) species are in the training data.

## Limitations
- Protonation enumeration appears limited to one rule-based state (+ deposited form); no multi-state ensemble at a pH, no pKa/populations.
- Relative energies reset per stereo block; no cross-protomer free-energy information.
- 14 elements (AIMNet2 set).

## Relevance to p(z,x|G,h)
Provides training pairs (z, x) for a *few* z per G (≈1.33 structures per ChEMBL ID on average, all variants included), but without labels p(z|G,pH). Could serve as a source of x|z data for multiple z; would need an external pKa/microstate-free-energy source (or a physics label) to supervise p(z|G,h). Its pH (implicitly 7.4) is a single condition, so it cannot by itself teach pH dependence.

## Exact supporting quotes
- README "Dataset Description": "Protomer generation (OpenEye FixpKa) and merge with original structures."
- README: "Unique structures: `~2.4M` (`~1.8M` unique ChEMBL compound IDs)"
- README "Zarr Format": "`charge`: `(N,)` `int8` - molecular charge value stored by the pipeline."
- README "Zarr Format": "`relative_energy` … energy shift relative to the lowest conformer in the corresponding `(mol_id, stereo_id)` block."
- KiltHub description: "a multi-stage workflow including protomer/stereoisomer handling, conformer generation, geometry optimization with an AIMNet2 potential"
- OpenEye FixpKa docs (overview): "the molecule is set to its most energetically favorable ionization state for pH=7.4."
- LoQI abstract (ChemRxiv): "This dataset captures complex organic molecules in various protonation states and stereochemical configurations."
