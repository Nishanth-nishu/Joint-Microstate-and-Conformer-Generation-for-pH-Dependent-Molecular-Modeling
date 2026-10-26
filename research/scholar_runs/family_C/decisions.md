# family_C decisions (2026-10-26 UTC)

- Google Scholar not searched: no tool in this session (see CONVENTIONS.md). Recorded as
  `not_searched_no_tool` in queries.csv. All coverage is OpenAlex + general web + direct fetches.
- OpenAlex title search is poor for tool papers: Epik 7, MolGpKa, MolTaut, pKAI, Dimorphite-DL
  (partially), Uni-pKa's published JACS Au version and Swails/Roitberg CpHMD were NOT retrieved by
  title words. Backtracked to WebSearch + publisher/PMC/GitHub fetches for each of these.
- ACS (pubs.acs.org) and Nature return 403/login to WebFetch. Backtracked to PMC mirrors
  (PMC11423309 Uni-pKa, PMC11325546 QupKake, PMC7904668 SAMPL6, PMC8409483 Wieder, PMC3968888 Over,
  PMC9023519 GEOM). ChemRxiv article pages also 403 → used ChemRxiv listing snippets only, and
  marked those papers abstract_only.
- Epik 7 kept at relevance 3 but fulltext_status=abstract_only: method details come from the
  abstract and Schrödinger's own technical page. Flagged in the note so no peer-reviewed
  methodological claim is attributed to it.
- Dropped as off-topic after inspection: protein pKa / CpHMD benchmark literature (nitrophorin,
  staphylococcal nuclease, Hv1, cysteine pKa). Protein CpHMD is the mature form of
  "coupled protonation + conformational sampling", but the systems are not drug-like small
  molecules; kept only three representative entries (Donnini, phbuilder, Henderson) at relevance 2.
- Searched specifically for a drug-like *small-molecule* constant-pH MD benchmark (multiple
  queries). None found beyond dicarboxylic acids in a thesis. Recorded as an open gap, not as
  evidence of absence.
- Could not isolate a quantitative "conformer search on/off" ablation for any ML pKa model. The
  only quantitative ablations found are (a) Jaguar++ 10-conformer weighting in the Uni-pKa paper
  and (b) Wieder 2021's multiple-minima vs global-minimum comparison. Both recorded with numbers.
