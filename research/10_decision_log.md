# Decision log (UTC)

- 2026-10-26 08:18 — Phase 0: TorDiff-NAX copied with fresh single-author history (user asked to remove all contributors; no co-author trailers).
- 2026-10-26 08:20 — scholar-scout skill unavailable on the Linux login node (path given is a Windows path). Substituted OpenAlex (reference-lookup MCP), web search, WebFetch. Google Scholar recorded as not searched for every family.
- 2026-10-26 08:25 — Baseline frozen at /scratch/nishanth.r/tordiff; new work goes in /scratch/nishanth.r/jmcg only.
- 2026-10-26 — Phase I finding: TD's torsion set is state-dependent (utils/torsion.py:9-38 rotatable iff the cut fragment has >1 atom) → p(z,τ|G,h) is trans-dimensional; joint model must use per-state torsion sets or a heavy-atom torsion parameterisation.
