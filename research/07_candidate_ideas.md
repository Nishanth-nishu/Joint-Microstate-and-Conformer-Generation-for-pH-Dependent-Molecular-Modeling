# Phase VII (initial) — Candidate research ideas

Status 2026-10-26: **initial set after the first literature/data pass (Phases II, IV, VI round 1). None is advanced.**
Phase VIII novelty loops (intent file + targeted search + adversarial judge) have NOT yet been run on any idea.
Evidence pointers: `scholar_runs/<family>/summary.md`, `paper_notes/*.md`, `04_loqi_rapa_audit.md`, `06_dataset_audit.md`.

## What the first pass established (evidence → constraint on ideas)

| Finding | Source | Consequence |
|---|---|---|
| No paper found that jointly *samples* protonation/tautomer state and 3D geometry conditioned on pH | family_A, family_B summaries (absence in OpenAlex + web; Google Scholar not searched) | Joint p(z,x\|G,h) not yet ruled out — **provisional** |
| AcepKa (JCIM 2026, arXiv 2604.00841): SMILES + pH → rule-enumerated microstates → Uni-Mol ΔG → Boltzmann weights at pH → RDKit/MMFF conformers per state | paper_notes/acepka.md (full text) | The **sequential pipeline already exists and is packaged**; it is the baseline, and the joint idea must beat it on something measurable |
| Uni-pKa / Epik 7 / QupKake solve p(z\|G,h) (task B) well without real 3D coupling (Epik 7 explicitly 2D; MAE ≈ 0.4–0.8 pKa) | family_C summary, paper_notes/epik7_2023.md, unipka_2024.md | No novelty or likely gain on task B alone |
| LoQI generates conformers for a *given* (possibly charged) state; no pH, no selection, no probabilities; ChEMBL3D mixes FixpKa protomers (≈pH 7.4) + deposited forms | 04_loqi_rapa_audit.md, paper_notes/loqi.md (code-based; PDF inaccessible) | p(x\|z,G) with charged states is done; ChEMBL3D gives no pH labels |
| RAPA: protein side-chain flip × protonation enumeration within 1 kcal/mol; no pH, no populations, no generation | paper_notes/rapa.md | Establishes value of returning *sets* of coupled state+structure alternatives; protein-specific |
| Conformation dependence of microstate populations is real but modest for typical drug-like molecules (~0.5–1 kcal/mol; ~1 pKa unit in IMHB cases) and dominated by finding the right minimum | paper_notes/over_2014_imhb_pka.md, wieder_2021_tautomer_ratios.md | Geometry→state coupling matters only in a **tail** (IMHB-capable, zwitterions, macrocycles, multi-site) — must be shown on a stratified subset |
| No experimental supervision for task C (equilibrium populations of (z,x) pairs); macro-pKa (IUPAC, SAMPL) and Tautobase pairwise ratios are the only experimental population-like labels | 06_dataset_audit.md | Population claims must be restricted to pKa/tautomer-ratio recovery; any p(z,x) reference is computed |
| GEOM-QM9 neutral as given; GEOM-DRUGS one protomer/SMILES, no pH; QMugs neutralised | 06_dataset_audit.md, paper_notes/geom_2022_dataset.md | QM9 TD reproduction is geometry reference only |
| TD's torsion set changes with protonation state | 01_baseline_audit.md §4.4 | τ-space is trans-dimensional across z |
| FlexiFlow / MiDi / Multiflow / PepFlow already provide discrete+continuous (incl. formal-charge categorical, torus+simplex) joint generation | family_D summary | A "discrete–continuous diffusion for protonation" mechanism is a **known mechanism in a new variable**; novelty can't rest on the mechanism |

## Candidate ideas (distinct hypotheses)

### I1 — Joint pH-conditioned state–conformer generator (the original idea)
- Hypothesis: a single model sampling (z, τ) jointly, conditioned on pH, gives better state–conformer consistency than AcepKa-style sequential pipelines.
- Closest: AcepKa (sequential, same interface), FlexiFlow/MiDi (mechanism), Uni-pKa (p(z|G,h)), LoQI (p(x|z)).
- Risk: with available data, z-labels come from the *same* predictors the baseline uses (Uni-pKa/FixpKa/Epik), so the joint model can at best imitate them; no geometry-dependent state labels exist. **Reject condition:** joint ≤ sequential on every metric that does not reduce to the external state predictor's own accuracy.
- Initial assessment: **weakly supervised; likely "combination of known components" unless I2/I3 supply a coupled signal.**

### I2 — Geometry-aware microstate free energies (conformer-ensemble-averaged state energies)
- Hypothesis: computing per-microstate free energy as an ensemble average over *generated* state-specific conformers (TD/LoQI samples + MLIP/implicit-solvent energies) improves pKa/tautomer-ratio prediction over single-conformer or 2D models **on the conformation-sensitive tail** (IMHB, zwitterions, flexible polyamines).
- Data: IUPAC pKa (CC BY-NC), SAMPL6/7/8, Tautobase (experimental); generator = existing TD/LoQI; energies = AIMNet2 / GFN2-xTB + ALPB.
- Closest: Uni-pKa (conformers only as augmentation), QupKake (xTB features), Wieder 2021 (tautomers, configurational sampling with MD), Jaguar++.
- Difference: generative-model conformer ensembles per microstate as the integration measure for a learned/physics free-energy at pH; measurable on experimental labels.
- Reject: no improvement over Uni-pKa on the stratified tail with matched data; or gain also obtained by RDKit+MMFF ensembles (then the generator is irrelevant).

### I3 — State-consistent conformer generation for charged/ionisable drug-like molecules (structural benchmark)
- Hypothesis: conformer generators trained/evaluated on one protomer misrepresent geometry of the pH-relevant protomer (e.g. salt bridges/IMHBs in zwitterions), and a state-conditioned generator trained on state-correct ensembles (CREST/xTB-ALPB regenerated per microstate) closes that gap.
- Closest: LoQI (already handles charged inputs — strong overlap), DiTMC/ET-Flow (charge features).
- Reject: TD/LoQI given the correct state already match state-correct references (then it's purely a data issue, not a method).

### I4 — Evaluation/methodology contribution: a pH-aware state–conformer benchmark
- Hypothesis: no benchmark evaluates (state, conformer) pairs at pH; existing conformer benchmarks silently fix one protomer, and pKa benchmarks ignore geometry. A benchmark with (a) structural metrics per microstate, (b) experimental pKa/tautomer-ratio recovery via ensemble free energies, (c) conformation-sensitive stratum, is itself a contribution.
- Closest: GEOM benchmarks, SAMPL challenges, AcepKa evaluation.
- Reject: if AcepKa or a 2025–26 benchmark already defines this (to search in Phase VIII).

### I5 — Hybrid: external enumerator + learned geometry-conditioned state re-weighting
- Hypothesis: keep a reliable enumerator (Dimorphite-DL / AcepKa rules / Uni-pKa) for z-candidates, generate state-specific conformers, and learn a small geometry-conditioned correction ΔG(z | ensemble) that is trained on experimental pKa — "learn when coupling matters".
- This is the learned version of I2; Candidate C ("coupled discrete–continuous generation") only if I5 shows geometry signal.
- Reject: learned correction ≈ 0 (no geometry signal) or overfits the small experimental set (scaffold split).

## Working ranking (provisional, before Phase VIII)
1. I2/I5 (geometry-aware microstate populations on the conformation-sensitive tail) — has experimental supervision, testable, cheap on one 3090 (generation + MLIP energies, small learned head).
2. I4 (benchmark) — likely publishable regardless of method outcome; packages negative results honestly.
3. I1 (full joint generator) — retained only as Candidate C if I5 finds geometry→state signal.
4. I3 — high overlap with LoQI; keep as a sub-experiment.
