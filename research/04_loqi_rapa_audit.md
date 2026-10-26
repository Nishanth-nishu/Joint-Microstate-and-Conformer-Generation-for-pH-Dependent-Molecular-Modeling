# 04 — LoQI / ChEMBL3D and RAPA audit vs. the proposed joint model p(z, x | G, h)

Date: 2026-10-26 08:26 UTC (13:56 IST). Run dir: research/scholar_runs/loqi_rapa/. Notes: paper_notes/loqi.md, paper_notes/chembl3d.md, paper_notes/rapa.md.
Google Scholar: not searched (no tool in this session).

Notation: G = constitution graph (heavy-atom skeleton / parent compound), z = protonation microstate (explicit-H graph + formal charges), x = 3D coordinates, h = conditions (pH; possibly solvent).

## Evidence base and its limits
- LoQI: **code read** (isayevlab/LoQI HEAD 04f3caef, 2026-10-01); **paper abstract only** (ChemRxiv PDF blocked by Cloudflare/403, OpenReview verification page; not bypassed). Mechanism claims come from code, configs and dataset README.
- ChEMBL3D: KiltHub README + readme.txt + loader scripts + small processed-split statistics read; topologies/zarr (≥1.6 GB / ~161 GB) not downloaded, so no per-molecule inspection of protonation variants.
- RAPA: full text read (PMC / Europe PMC XML); SI and code not read.
- TD, sequential pipelines, constant-pH MD columns: not re-audited in this family; labels rest on well-known primary descriptions and are marked as such. Re-verify against their own paper notes before citing.

## Direct answers — LoQI / ChEMBL3D
1. **Representation of protonation states:** implicitly, as the explicit-H molecular graph plus per-atom formal charges. Charges are a 6-class one-hot (−2…+3) concatenated to the 17-element atom one-hot → 23 atom classes (loqi_flow.yaml:71-81; module.py:270-287; featurize.py:74-75).
2. **Multiple states per molecule in training data:** Partially supported. ChEMBL3D README step 1: "Protomer generation (OpenEye FixpKa) and merge with original structures"; ~2.4M structures for ~1.8M ChEMBL IDs. LoQI keeps every variant `mol_id` (`CHEMBL123_k`) as a separate example, grouped by parent only for splitting (process_chembl3d.py:130-187). Processed train split contains charged atoms (N +1 8.2%, O −1 4.1%). Not verified for individual molecules.
3. **Tool / pH:** OpenEye FixpKa. pH not stated by the authors; OpenEye documents FixpKa as setting "its most energetically favorable ionization state for pH=7.4" → at most ~one rule-based state per input plus the deposited form. pH value used by ChEMBL3D = Unknown (resolve via paper SI or authors).
4. **Inference requires a supplied state:** Confirmed. Input SMILES → RDKit AddHs → fixed graph/charges (validation.py:61; featurize.py:53-97). Output conformers have exactly the input topology (featurize.py:239-253).
5. **Model predicts/selects a state:** Not supported by available evidence. Atom types, bonds and charges are `discrete_null` ("pass in and not noise/denoise", builder.py:188-190); loss on coordinates only (loqi_flow.yaml:31-38). No protonation/pKa/tautomer code in the repository (grep).
6. **pH conditioning:** Not supported by available evidence (no pH argument in api.py:160-172, CLI, or configs).
7. **State probabilities / equilibrium populations:** Not supported by available evidence; ChEMBL3D has no population/free-energy field, and its `relative_energy` resets per (mol_id, stereo_id).
8. **H / charges / bond orders:** all H explicit and generated as nodes; formal charges fixed per atom; Kekulé bond orders 1/2/3 fixed, plus E/Z and tetrahedral stereo edge types. AIMNet2 scoring uses total formal charge (molecule_metrics_aimnet2.py:156).
9. **Exact model input:** {element one-hot ⊕ formal-charge one-hot per atom (incl. H)}, fully connected edges with {none/1/2/3/aromatic-unused + stereo types 5–8}, time t, noised coordinates x_t; output velocity/denoised x.
10. **Dataset contains several states vs. model selects a state:** the *dataset* (and LoQI's training set) contains more than one protonation variant for part of the compounds; the *model* never selects among them. It learns p(x | z) amortised over a mixture of z.

## Direct answers — RAPA
- **Problem:** H-invisible, N/O-ambiguous X-ray protein structures → ambiguous His/Asn/Gln/Ser/Thr/Tyr (and acid-dyad Asp/Glu) rotamer+protonation states; returns all energetically degenerate full-protein configurations instead of one (Abstract; §2.2.1).
- **Coupling:** each ambiguous residue has a discrete RP state = (rotamer flip × protonation); scored jointly with neighbours' RP states via an H-bond lookup table; branch on states within 1 kcal/mol (§2.2.1.3–2.2.1.5, §2.2.3). Not PDB altlocs (no altloc/occupancy mention); heavy-atom scaffold fixed.
- **Protein-specific:** residue templates, Lys/Arg always charged, acid-dyad rule, inter-residue H-bond networks on a crystal scaffold.
- **Transferable:** "set of accessible states, not one"; joint scoring of discrete protonation + local geometric alternatives; branching/self-consistent pruning; evidence that accessible sets are small (62/77 ≤8).
- **Generative 3D sampling:** No (enumerates discrete H positions/flips on fixed coordinates).
- **pH:** No (pKa corrections listed as future work, Discussion).
- **Establishes:** single-configuration preparation discards plausible states; joint protonation/rotamer assignment is tractable for proteins. **Does not establish:** pH dependence, population weights (authors say entropy needed), anything about small-molecule conformational ensembles or learned generative models.

## Comparison table
Labels: C = Confirmed, P = Partially supported, N = Not supported by available evidence, NA = Not applicable, U = Unknown. "bg" = background knowledge, not audited in this family.

| Aspect | Torsional Diffusion (TD) | LoQI (+ChEMBL3D) | RAPA | Sequential protonation→conformer (Epik / Dimorphite-DL / MolGpKa / Uni-pKa → RDKit/OMEGA/TD/LoQI) | Constant-pH MD/MC | Proposed joint p(z,x\|G,h) |
|---|---|---|---|---|---|---|
| Input representation | C: RDKit mol incl. H, torsion space on RDKit local structure (bg; TD paper arXiv 2206.01729) | C: explicit-H graph, element⊕formal-charge one-hot, Kekulé bonds + stereo edges (featurize.py:53-97) | C: one prepared protein PDB (§2.1) | C: SMILES/2D → enumerated microstate(s) → conformer tool (bg) | C: full-atom system + titratable-site models (bg) | Planned: G (parent graph) + h |
| Protonation-state handling | N: state fixed by input SMILES (bg) | C: fixed, user-supplied state; charges as features (loqi_flow.yaml:71-76) | C: discrete His/acid-dyad states per residue (§2.2.2) | C: explicit enumeration step (bg; tool-dependent) | C: sampled dynamically (bg) | Planned: z sampled |
| Several states in training data | U (GEOM inputs as deposited; not audited here) | P: FixpKa protomer + original, ~2.4M structures/1.8M IDs (ChEMBL3D README) | NA (no training) | NA / tool-specific (pKa models trained on pKa data) | NA | Required |
| State selection by model | N | N (discrete_null; builder.py:188-190) | C (deterministic within 1 kcal/mol window; §2.2.1.5) | C/P: rule- or pKa-based selection outside the conformer model (bg) | C (bg) | Planned |
| pH conditioning | N | N (no pH input; api.py:160-172); dataset implicitly ~pH 7.4 (FixpKa docs; U for ChEMBL3D itself) | N (Discussion: future work) | C for Epik/Dimorphite-DL/Uni-pKa pH inputs (bg) | C (bg) | Planned |
| State probabilities | N | N | N (degeneracy cutoff, no weights; Discussion) | P: Epik state penalties, Uni-pKa microstate free energies→populations; Dimorphite-DL none (bg) | C via sampling (bg) | Planned |
| Equilibrium populations of (z,x) | N (TD targets Boltzmann-like conformer ensembles for fixed z; bg) | N (trained on single lowest-energy conformer per mol_id; process_chembl3d.py:78-121) | N (authors: entropy term needed) | P: z-populations and x-ensembles from separate models, not coupled (bg) | C in principle, costly (bg) | Planned; must be validated |
| Generative 3D sampling | C (torsions) (bg) | C (all-atom coordinates, flow/diffusion) | N (discrete H/flip choices on fixed scaffold) | C in the conformer stage (bg) | C (MD/MC trajectories) (bg) | Planned |
| H / charge handling | P: explicit H from RDKit, charges as given (bg) | C: explicit H generated; formal charges fixed; AIMNet2 total charge = sum formal charges | C: H positions enumerated for Ser/Thr/Tyr; His/Asp/Glu protonation (§2.2.2) | C: H/charges set by protonation tool (bg) | C (bg) | Planned: H count/charges are part of z |
| Coupling of z and x | NA | N: x\|z only | P: discrete protonation coupled to local rotamer flips and neighbour H-bonding | N: z chosen without x (except 3D-aware pKa modes; U) | C (bg) | Planned core claim |
| Data | GEOM-DRUGS/QM9 (bg) | ChEMBL3D, AIMNet2 (ωB97M/CPCM-trained) optimised, ~283.6M conformers (readme.txt) | 77 DUD-E proteins; MD validation | pKa datasets + conformer libraries (bg) | force fields (bg) | TBD |
| Domain | small molecules | small/drug-like molecules (14 AIMNet2 elements) | proteins | small molecules | proteins and small molecules (bg) | small molecules |
| Code availability | C (bg) | C, MIT | C, MIT | C/varies; Epik commercial (bg) | C (Amber/GROMACS variants) (bg) | — |

## Overlap verdict
- **LoQI overlaps the x | z factor substantially:** charge-aware, explicit-H, QM-informed all-atom conformer generation conditioned on a given microstate, trained on a dataset that includes more than one protonation variant for many compounds. Any claim that "conformer generators ignore charged states" is false for LoQI.
- **LoQI does not cover z | G, h:** no state variable, no pH, no probabilities, no populations. Its training set treats each variant as an independent molecule and keeps only its lowest-energy conformer. The dataset offers no pH-dependent supervision for z (single FixpKa state, implicitly pH 7.4).
- **RAPA** is protein-side motivation (joint discrete protonation + structural alternative, multiple accepted states). It has no pH, no populations, and no generative sampling, and does not overlap a generative small-molecule model.
- **Novelty claim that survives:** a single learned model that samples z and x *jointly*, conditioned on pH, with calibrated p(z|G,h) and state-coupled conformer ensembles. **Mandatory baselines:** (i) FixpKa/Epik/Dimorphite-DL/Uni-pKa → LoQI, and (ii) the same → TD. The proposal must beat or match these on population calibration and conformer quality. Using ChEMBL3D for the x|z part is possible (MIT), but pH labels for z must come from elsewhere.
- No citing work extending LoQI or RAPA toward pH or state selection was found (OpenAlex 0 citing works for both; 4 web searches). AcepKa (arXiv 2604.00841, abstract only) is the nearest sequential pKa-ensemble + separate conformer-generator system and should be audited in the pipeline family.

## Open unknowns (and how to resolve them)
- ChEMBL3D pH value and what each `_k` suffix is: inspect one `topologies/NNN.sdf` (topologies.tar.zst, 1.6 GB) for a parent with ≥2 suffixes, or read the paper SI.
- Whether the LoQI paper reports per-protonation-state performance: needs paper PDF access (ChemRxiv via browser).
- Fraction of ChEMBL IDs with >1 protonation variant: count `parent_mol_id` multiplicity in `train_smiles.pickle` (238 MB) or the split manifest.
