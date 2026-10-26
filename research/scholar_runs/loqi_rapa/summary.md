# loqi_rapa — family summary (2026-10-26 08:26 UTC / 13:56 IST)

Backends: OpenAlex (lookup, search_works, citation_neighbours), WebSearch, WebFetch, git clone. Google Scholar: not searched (no tool).
Full audit: research/04_loqi_rapa_audit.md. Notes: paper_notes/loqi.md, chembl3d.md, rapa.md.

## Closest papers and what they actually do
- **LoQI** (Nikitin et al., ChemRxiv 2025; code read, paper abstract-only). Conditional all-atom conformer generator p(x | z): the input explicit-H graph has per-atom formal charges one-hot concatenated to element type, and bond orders and stereo edges are fixed. Only coordinates are noised and learned (charges/bonds/atoms `discrete_null`). There is no pH, no state selection, no probabilities. Training uses one lowest-energy conformer per ChEMBL3D `mol_id` and keeps all protonation variants as separate examples.
- **ChEMBL3D** (KiltHub 10.1184/R1/31428449). OpenEye FixpKa protomer + original structure, Flipper stereo, OMEGA conformers, AIMNet2 FIRE optimisation. About 283.6M conformers and ~2.4M structures for ~1.8M ChEMBL IDs. The schema has mol_id, stereo_id, energy, relative_energy and total charge, with no pH or population field. FixpKa is documented by OpenEye as a single most-favourable state at pH 7.4; the pH for ChEMBL3D itself is not stated.
- **RAPA** (Ghattas et al., JCIM 2025; full text). For proteins, it does deterministic, self-consistent enumeration of His/Asn/Gln/Ser/Thr/Tyr/acid-dyad rotamer+protonation states on a fixed crystal scaffold. Scoring uses an H-bond lookup table, and every configuration within 1 kcal/mol is kept. It has no pH, no populations, and no 3D generation.
- AcepKa (arXiv 2604.00841; abstract only): pKa protonation ensemble plus a separate conformer generator. It belongs to the pipeline family.

## Verdict
LoQI covers x|z, including charged states. Nothing found covers pH-conditioned z|G,h jointly with x. A protonation tool → LoQI pipeline is the mandatory strong baseline.

## Open unknowns
ChEMBL3D pH and suffix semantics; LoQI paper body (per-state results); multiplicity of variants per parent. Resolution paths are in 04_loqi_rapa_audit.md.
