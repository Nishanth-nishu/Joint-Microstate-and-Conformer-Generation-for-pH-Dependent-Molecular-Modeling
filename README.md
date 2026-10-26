# Joint Microstate and Conformer Generation for pH-Dependent Molecular Modeling

Research programme that asks whether a generative model can jointly learn chemical microstate selection
(protonation/tautomer state `z`) and 3D conformer sampling (`x`, or torsions `τ`), conditioned on pH `h`:
`p(z, x | G, h)`, and whether it beats sequential "pick a state, then generate conformers" pipelines.

The starting point is a completed reproduction of Torsional Diffusion on GEOM-QM9 (prior work, see `prior_work/`
and the round reports carried over from TorDiff-NAX). Nothing here assumes the joint idea is novel or correct:
every phase is literature- and data-audited first.

- `research/` — evidence system: search logs, paper registry, audits, candidate ideas, novelty claims, decisions.
- `prior_work/TORDIFF_NAX_HANDOFF.md` — handoff of the Torsional Diffusion baseline work.
- Baseline code, data and checkpoints live on gnode118 at `/scratch/nishanth.r/tordiff` (frozen, not modified).
