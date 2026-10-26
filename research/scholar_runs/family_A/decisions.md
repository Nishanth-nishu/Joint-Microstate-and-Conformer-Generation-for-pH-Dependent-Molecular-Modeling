# family_A decisions (2026-10-26, UTC)

- Google Scholar not searched: no tool in this session (recorded as `not_searched_no_tool` in queries.csv). Never claim otherwise.
- OpenAlex `search_works` matches on TITLE WORDS only. Queries containing "protonation"/"proton" returned proton
  therapy, proton-transport MD and particle physics: effectively 0 relevant hits for 3 separate query phrasings.
  Decision: stop using OpenAlex free-text for family A concepts; rely on WebSearch + targeted WebFetch, and use
  OpenAlex only to resolve known titles (Uni-pKa, Torsional Diffusion) and for citation_neighbours.
- `citation_neighbours` on Torsional Diffusion (90 citations, top 25 by citation count) returned no protonation/pH
  descendant. Decision: no second citation-expansion round on that node; the top-cited tail is dominated by generic
  diffusion surveys and CV papers.
- Searches for "tautomer-aware 3D generation", "protomer conformer generation", "protonation state docking
  diffusion", "co-folding ligand protonation" all returned pipelines or unrelated 3D generators. Decision: record as
  a negative result rather than keep re-phrasing; 4 distinct phrasings is enough evidence of absence to report
  honestly (absence in these indexes, not proof of non-existence).
- AcepKa promoted to relevance 3 and read in full: it is the only found work that takes pH as an explicit input and
  emits both microstates and conformers. Decision: treat it as the concrete sequential baseline, not as prior art
  for joint modelling.
- LoQI/ChEMBL3D kept at relevance 2 despite being conformer-only, because its stereochemistry-conditioning trick and
  its protonation-state-containing training data are the nearest conditioning precedent on the x-side.
- Dropped as off-topic after one look: protein CpHMD/pKa predictors (DeepKa, PypKa, Titr-DMD), protein conformational
  ensemble generators (ConfDiff, OpenComplex2), MolGpKa-style pKa regressors with no 3D generation.
- CORDIS project 101210115 (flow matching for protonation + multiconformer PROTEIN structures) recorded at
  relevance 1: a funded project description, not a published method, and protein- not ligand-scoped. Worth
  re-checking for outputs later; a competitor on the adjacent protein problem.
- No PDF text extractor exists on this host (no pdftotext, no pypdf, no pip). Decision: full text read via arXiv
  HTML through WebFetch; PDFs still archived to research/papers/ for later. fulltext_status is recorded honestly:
  only AcepKa and FlexiFlow were read in full for family A.
