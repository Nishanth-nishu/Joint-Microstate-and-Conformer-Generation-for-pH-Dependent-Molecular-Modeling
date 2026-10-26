# family_B decisions (backtracking log) — all times UTC, 2026-10-26

- Google Scholar was NOT searched: scout_cli.py / Scholar access is unavailable on this Linux session. One
  `google_scholar,not_searched_no_tool` row is recorded in queries.csv. No claim of Scholar coverage anywhere.
- "COSMIC" is unresolvable by general web search (name collides with cosmology). Backtracked to OpenAlex
  search_works, which returned W4395663930 (JCIM 2024, Insilico Medicine). Kept at relevance 2.
- COSMIC full text not obtained: pubs.acs.org returns HTTP 403 and no preprint exists. Marked not_accessible;
  every mechanism row for COSMIC is Unknown rather than a negative label.
- "PDCG" likewise absent from web search; resolved via OpenAlex to Chemistry (MDPI) 2026, doi 10.3390/chemistry8020029.
  MDPI blocks automated fetching (403) despite being open access, so PDCG stays metadata-only, relevance 1.
- "MSGEN" has no arXiv entry. Backtracked to the ICLR 2026 poster page, which gave the abstract and authors.
  OpenReview PDF blocked (403 bot check). Abstract-only; no method details asserted.
- No PDF text extractor exists in this environment (no pdftotext, no pypdf, no PyMuPDF). Switched to arXiv HTML
  (arxiv.org/html/<id>) + a local HTML-to-text script rather than downloading unreadable PDFs into research/papers/.
  research/papers/ is therefore empty for this family; all full texts were read as HTML.
- GEOM full text: nature.com redirects to an IdP authorize URL and PMC9050668 resolves to an unrelated article.
  Backtracked to arXiv:2006.05531v4 HTML, which carries the full Methods including the de-salting/ionization text.
- GeoDiff: raw.githubusercontent timed out on utils/datasets.py twice; only utils/chem.py (bond-type constants) was
  retrieved. Refused to infer GeoDiff's node features from ConfGF's; left Unknown with the resolving file named.
- ET-Flow: arxiv.org/html/2410.22388v2 returns 404 and the abstract page carries no featurization. Backtracked to
  the official repo, where etflow/commons/featurization.py confirms AddHs + GetFormalCharge. Paper itself stays
  abstract_only, so dataset/split claims for ET-Flow are taken from DiTMC and SO(3)-AvgFlow, which reproduce it.
- Uni-Mol featurization not resolved (search returned venue/code but no feature statement). Left Unknown rather
  than guessed; the resolving path (unimol/data/ preprocessing) is recorded in the note.
- LoQI listed in candidates.csv and as a column of feature_matrix.csv (all Unknown) but deliberately NOT audited —
  another agent owns it. Existing notes loqi.md / loqi_chembl3d.md were not opened or relied on.
- A note tordiff.md already existed in research/paper_notes/ from another family. Wrote tordiff_jing2022.md rather
  than overwriting it; the two may overlap and should be reconciled before the write-up.
- Searched explicitly for 2025-2026 pH- or protonation-conditioned conformer generation ("conformer generation flow
  matching 2026 protonation state charge aware") and found none. SemlaFlow (AISTATS 2025) surfaced as the nearest
  miss — it generates formal charge jointly with coordinates, but it is de novo 3D molecule generation, not
  conformer generation for a given graph, so it is recorded here and left to family A/C rather than audited.
