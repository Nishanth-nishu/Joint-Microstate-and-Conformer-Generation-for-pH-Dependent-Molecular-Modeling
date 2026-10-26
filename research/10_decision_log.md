# Decision log (UTC)

- 2026-10-26 08:18 — Phase 0: TorDiff-NAX copied with fresh single-author history (user asked to remove all contributors; no co-author trailers).
- 2026-10-26 08:20 — scholar-scout skill unavailable on the Linux login node (path given is a Windows path). Substituted OpenAlex (reference-lookup MCP), web search, WebFetch. Google Scholar recorded as not searched for every family.
- 2026-10-26 08:25 — Baseline frozen at /scratch/nishanth.r/tordiff; new work goes in /scratch/nishanth.r/jmcg only.
- 2026-10-26 — Phase I finding: TD's torsion set is state-dependent (utils/torsion.py:9-38 rotatable iff the cut fragment has >1 atom) → p(z,τ|G,h) is trans-dimensional; joint model must use per-state torsion sets or a heavy-atom torsion parameterisation.
- 2026-10-26 ~08:50 — Literature round 1 complete (157 queries, 78 registry papers: 27 full text, 48 abstract-only, 3 inaccessible). No PDF text extractor on login node → full text read via arXiv/PMC HTML. Citation expansion thin (4 citation_neighbours calls) — must be extended in Phase III/VIII.
- 2026-10-26 — AcepKa (arXiv 2604.00841) identified as the closest system: same interface (SMILES+pH → weighted microstates + conformers), sequential. Becomes mandatory Baseline A.
- 2026-10-26 — Task B (p(z|G,h)) judged largely solved by 2D tools (Epik 7, Uni-pKa); task C has no experimental supervision. Initial direction shifted from "joint generator" (I1) toward geometry-aware microstate populations on a conformation-sensitive stratum (I2/I5) + benchmark (I4). Provisional, pending Phase VIII.
- 2026-10-26 — LoQI audit based on code + dataset docs; paper PDF inaccessible (ChemRxiv Cloudflare 403). Must re-attempt full text before citing methodological details from the paper.
