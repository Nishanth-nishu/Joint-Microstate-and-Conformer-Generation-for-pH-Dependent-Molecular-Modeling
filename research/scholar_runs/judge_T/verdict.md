# Judge verdict T1/T2 (independent judge, 2026-10-27 10:36 UTC / 16:06 IST; transcribed by orchestrator)

T1 — B Revise (demote to diagnostic). The softmax-gradient argument is elementary; Uni-pKa's authors already report the
symptom: "the model tends to attach similar free energies to different tautomers and fails to distinguish the tautomers with
drastic energy differences" and "Tautomerism has not been explicitly addressed in our framework" (PMC11423309). Stronger
form: an LSE macro-pKa loss is exactly invariant to any within-macrostate redistribution that keeps the LSE fixed;
identification only via parameter sharing. Narrowed: "Released Uni-pKa's aqueous Tautobase errors are dominated by
compression of within-macrostate tautomer free-energy gaps (regression slope of predicted on measured log K_T well below
1), the signature expected because a log-sum-exp macro-pKa loss leaves within-macrostate free-energy differences
unconstrained except through parameter sharing." Falsifier: zero-shot MAE < 1.0; or slope CI includes 1; or MAE not
≥ predict-zero + 0.3 while slope ≥ 0.7.

T2 — B Revise, then pilot. Zankov 2019 JCIM (abstract + refs only; KPFU copy) did conjugated pKa/log K_T learning with
the Kabachnik relation (twin NN); joint-vs-pKa-only alone = obvious multitask. Narrowed: "On held-out (cluster-split)
measured aqueous Tautobase pairs, a Uni-pKa free-energy head fitted jointly to Dwar-iBonD macro-pKa and Tautobase log K_T
achieves lower tautomer MAE than the same head fitted to Tautobase alone, while macro-pKa MAE on SAMPL6/7/8 + Novartis
rises by ≤ 0.05 and SAMPL6 micro-pKa MAE does not rise." Falsifier: CI of Δ(joint − tautomer-only) includes 0; pKa rise
> 0.05; joint MAE > 1.0 log K. For a paper (not just a repair): within 0.2 log K of sPhysNet-Taut on the same pairs.
Risks: frozen features may not separate tautomers (head-only dead → backbone fine-tune); leakage of pKa-derived 'E' pairs.
Uni-pKa pretraining labels: ChEMBL predicted pKa through the same LSE loss (no within-macrostate signal).
Unknowns flipping the verdict: Zankov full text (macro labels + ensemble → T2 toward C); Starling/Epik 7 tautomer data;
Tautobase licence.
