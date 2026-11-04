# FlexiTors V2: research decision (2026-11-04)

Inputs:
- `FINAL_RESULTS.md`: evaluation 57155, run once;
- `FINAL_SCIENTIFIC_ASSESSMENT.md`;
- `LITERATURE_COMPARISON.md`;
- the decision rule in DESIGN.md addendum B, fixed 2026-10-31 19:50 UTC before any T2 test number.

## 1. The pre-declared rule
| addendum-B condition | value | met? |
|---|---|---|
| endpoint 1 (improves on TD) | "yes": −63 mÅ, 3 / 3 seeds, both TD seeds | yes |
| endpoint 4 (vs constrained xTB) | "worse": +9.3 mÅ [+3.6, +14.1] (TD seed 1: +14.2) | — |
| head no more expensive than xTB per conformer | 1.23 vs 0.30 CPU s (like for like): 4.1× more | no |
| a demonstrated use xTB cannot serve | none shown (torsion preservation is matched by restrained xTB + restoration) | no |

How each rule option fares:
- **A (continue developing the head):** needs endpoint 4 "beats" or "indistinguishable". Not met.
- **B (optional component):** needs the head cheaper than xTB, or a use xTB cannot serve. Not met.
- **C (stop the head, redirect to microstate–conformer):** endpoint 4 is "worse" and the head is not cheaper.
  **Met.**
- **D (one more narrow experiment):** must be able to change the choice among A–C. It cannot:
  - **Longer training:** endpoint 5 is a clear but small gain (−3.9 mÅ per 10 epochs), and validation still falls.
    Even if the 9–14 mÅ gap to xTB closed, the head would still cost 4× more than xTB per conformer, so neither A nor
    B would be met.
  - **Generalisation to GEOM-DRUGS:** not feasible here (no data or model), and the rule excludes it.

**Rule outcome: C.**

## 2. The brief's three options, on its criteria
- **A:** methods paper.
- **B:** benchmark / evaluation paper ("does reference-consistent physical optimisation change conclusions about
  learned refiners?").
- **C:** freeze FlexiTors, prioritise joint microstate–conformer generation.

| criterion | A | B | C |
|---|---|---|---|
| Scientific novelty | Low. Refining TD samples exists (DiSCO). The new constraint (P_W) is not rewarded: projector effect ≤ 1 mÅ, H2 failed. | Moderate but narrow. A torsion-restrained physics control on the same samples was not found in the checked works, and it changes the conclusion on QM9. | The novelty lies in the microstate–conformer project. FlexiTors contributes reusable tools. |
| Strength of evidence | Against A: dominated by xTB on accuracy, validity and cost (H4, H5 failed). | Strong for one refiner, one dataset, one TD checkpoint. A benchmark claim needs ≥ 2–3 refiners and a second dataset. | The evidence is complete and documented; nothing is pending. |
| Publication potential | Not supportable: no advantage to report. | A workshop paper or short note at best now. A full paper needs DiSCO and ≥ 1 other refiner plus GEOM-DRUGS. | FlexiTors is a negative result and protocol, citable as a technical report or appendix. |
| Remaining experiments | None would change the verdict (see D above). | Run DiSCO (public code) and ideally ET-Flow-class refiners through the same restrained-xTB control; add GEOM-DRUGS (data and models not on the cluster); vary k. | None for FlexiTors. |
| Compute | — | DRUGS-scale xTB and model inference; GPU quota currently exhausted. Weeks of CPU at QM9 scale, far more for DRUGS. | Frees gnode118 CPU for the main project. |
| Time to completion | — | 1–3 months, dominated by obtaining DRUGS data and models and GPU time. | Immediate. |
| Risk of inconclusive results | — | High: the QM9 finding may be an artefact of xTB references and may not transfer to DRUGS (whose references are CREST/GFN2 too, but molecules are larger, where cost trade-offs differ). | Low for this decision. |
| Relation to long-term direction | Off-direction. | Tangential: about evaluating refiners, not about protonation-state–conformer coupling. | On-direction. The projector, pure-L data pipeline, constrained-xTB protocol and audit tooling carry over. |

## 3. Recommendation
**Option C: freeze FlexiTors V2 as a completed investigation and prioritise the joint microstate–conformer work.**

Why:
1. It is what the pre-declared rule gives. The rule was written before the results, and there is no basis for
   overriding it.
2. A methods paper (A) is not supported. The head is dominated by a cheaper physical step on every pre-declared
   measure. The one architectural novelty, the projector, has no measurable benefit.
3. A benchmark paper (B) is the only publishable framing, but its claim currently rests on one refiner and one
   dataset whose references are xTB geometries. Making it a general claim needs other refiners and DRUGS-scale data,
   which are not available here. The brief also forbids generalising beyond QM9 without that evidence.

**What to keep:**
- **Reports:** the FlexiTors records (design, results, this decision), as a self-contained negative result.
- **Reusable pieces:**
  - the torsion-restrained xTB control (the baseline any geometry method must beat);
  - the projector and its tests;
  - the pure-L target construction;
  - the run-once evaluation and audit tooling.
- **A possible later note:** B remains open as a short note if DiSCO is ever run through the same control. That
  would be a separate, pre-declared experiment, not part of this decision.

**What not to claim:**
- superiority over xTB;
- any benefit of the projector;
- new torsional modes;
- anything beyond GEOM-QM9.
