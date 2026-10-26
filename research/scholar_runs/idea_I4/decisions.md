# idea_I4 — backtracking decisions (2026-10-26, UTC; IST = UTC+5:30)

- 10:17 UTC: Google Scholar not searched (no tool); one `not_searched_no_tool` row recorded.
- 10:18 UTC: Two OpenAlex title searches hit HTTP 429; by ~10:22 UTC the IP's free daily OpenAlex budget was exhausted (resets 00:00 UTC / 05:30 IST). Citation expansion of SAMPL6 switched to the Semantic Scholar API (logged as backend webfetch). GEOM citers only available as top-25-by-citations from the earlier MCP call — low-cited 2025–26 GEOM citers NOT enumerated (open gap).
- 10:20 UTC: AcepKa (arXiv + JCIM) has 0 references/citers recorded in OpenAlex → citation_neighbours uninformative; relied on its full text instead.
- 10:20 UTC: "ligand preparation protonation tautomer benchmark" in OpenAlex returned mostly organometallic tautomerism (off-topic) → moved to web search.
- 10:22 UTC: PEACE (Zheng & Heid 2026) ChemRxiv full text 403; kept as abstract_only. Its abstract describes a METHOD (tautomer/protomer populations with conformational effects) with "practical benchmarks", not a (state × conformer) benchmark → classed same-objective-different-mechanism for the I2/I5 methods, "related" for I4.
- 10:25 UTC: Keserű pKaHub benchmark PDF unreadable by WebFetch; full text obtained from Europe PMC fullTextXML (PMC13126629). No 3D/conformer content → overlaps only the pKa-recovery component of I4.
- 10:30 UTC: Starling (Rowan) found via web search; full text is HTML on rowansci.com → read. Changes earlier belief that "no ML pKa model aggregates conformers thermodynamically" (see summary).
- 10:35 UTC: Uni-pKa conformer handling resolved from code instead of SI (SI not read).
- 10:45 UTC: Charged/zwitterion conformer-generator evaluation: no hit in OpenAlex (429), web, or arXiv API. Recorded as absence-of-evidence, not evidence of absence.
- 10:50 UTC: LoQI, Epik 7 full text still inaccessible (Cloudflare/challenge/paywall) — not bypassed.
- 11:05 UTC: Stopped expanding: last 6 query families returned no new candidate with relevance ≥2 for the benchmark (as opposed to method) claim.
