# Research workspace conventions (read before writing anything here)

Staging copy: /home2/nishanth.r/jmcg_stage/research (shared NFS). It is rsynced into
/scratch/nishanth.r/jmcg/research on gnode118 and pushed from there after each phase.

## Tooling limitation (recorded 2026-10-26)
The scholar-scout skill (scout_cli.py at C:/Users/HP/...) is NOT installed on this Linux session. Substitute:
- OpenAlex via reference-lookup MCP: search_works (title words), citation_neighbours (refs + citing), lookup_reference.
- Web search (general web; NOT Google Scholar) and WebFetch for arXiv/ChemRxiv/PMC/publisher/GitHub pages.
- Google Scholar was NOT searched directly by any query. Never claim it was.
Every query row records backend statuses honestly.

## One run dir per search family
research/scholar_runs/<family>/ containing:
- queries.csv  : header `timestamp_utc,family,query,backend,status,n_results,notes`
                 backend ∈ {openalex,websearch,webfetch,citation_neighbours,google_scholar}
                 add one google_scholar row per family with status `not_searched_no_tool`.
- candidates.csv : header `id,title,year,venue,doi_or_url,code_url,relevance(0-3),fulltext_status,checked_utc`
                 fulltext_status ∈ {full_text_read, abstract_only, not_accessible}
- decisions.md : backtracking decisions (few candidates, weak relevance, off-topic) one line each.
- summary.md   : family findings — closest papers, what each actually does (mechanism level), open unknowns.

## Paper notes
research/paper_notes/<short_id>.md with sections: Citation (title, authors, year, venue, DOI/URL, code, date checked),
Fulltext status, Problem, Method (mechanism level; equations/sections cited with section/page),
Data, Inputs/outputs (what the model takes and returns — esp. protonation/charge/H handling, pH),
Experiments & baselines, Limitations, Relevance to p(z,x|G,h), Exact supporting quotes (≤40 words each, with section/page).
Abstract-only papers: say so at the top; do NOT state methodological details from them.

## PDFs
Download open-access PDFs to research/papers/<short_id>.pdf (git-ignored: copyright). Never bypass paywalls/CAPTCHAs.

## Honesty
Do not label a method as lacking a feature because the abstract doesn't mention it; write `Unknown` + what would resolve it.
Do not invent details. Dates in UTC.
