# idea_I4 summary — written by orchestrator from the search agent's hand-back, 2026-10-26
No paper found that scores (microstate, conformer) pairs at pH, tests conformer generators on protomers/charged forms, or stratifies ML pKa error by conformational sensitivity. Overall class: RELATED.
1. AcepKa (full text) — builds the pairs; reports only macro-pKa RMSE and molecule-level conformer RMSD (Platinum 2017). Conformers/microstate and aggregation rule not stated. No licence; web app only.
2. Keserű pKaHub benchmark (JCIM 2026, Europe PMC full text) — 90k pKa, 7 predictors, no 3D. Same objective, different mechanism (pKa only).
3. Agbaglo/Janesko JPCA 2025 (abstract only) — conformational-sampling protocols for pKa, 20 SAMPL7 molecules, QM.
4. SAMPL6/euroSAMPL1 — states, no geometry.  5. GEOM — one fixed protomer.
Also: PEACE (Zheng & Heid 2026, abstract only), FlexiSol (Grimme 2025).
Attack: publishable alone only at JCIM / J Cheminform / LiveCoMS (NeurIPS D&B if ML generators central), and only with a measured failure (generators degrade on charged states, or pKa error rises with conformational sensitivity). Redundant if PEACE / LoQI / GEOM-Drugs-Revisited full texts already stratify by charge. Low-cited GEOM citers not checked (OpenAlex budget).
Full-text retries: Uni-pKa code read (Apache-2.0; arithmetic mean over 10 conformers); Starling uses Boltzmann-style LSE over 10 conformers/microstate; LoQI and Epik 7 still inaccessible.
