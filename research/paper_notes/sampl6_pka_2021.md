# SAMPL6 pKa challenge overview

## Citation
Işık M., Rustenburg A. S., Rizzi A., Gunner M. R., Mobley D. L., Chodera J. D.
"Overview of the SAMPL6 pKa challenge: evaluating small molecule microscopic and macroscopic pKa
predictions." J. Comput. Aided Mol. Des. 35:131 (2021). DOI 10.1007/s10822-020-00362-6 (PMC7904668).
Measurements: Işık et al., JCAMD 2018, DOI 10.1007/s10822-018-0168-0. Checked 2026-10-26.

## Fulltext status
full_text_read (PMC, first ~100k chars).

## Data
24 kinase-inhibitor-like fragments. "Experimental macroscopic pKa measurements were collected using
a spectrophotometric method with the Sirius T3" (Sec. 2.1) — i.e. **macroscopic only**.
"For a subset of eight molecules, the dominant microstate sequence could be inferred from NMR
experiments" (Sec. 2.3) — microstate information exists only as a *dominant-state label* for 8
molecules, not as populations.

## Evaluation
Macro-pKa by RMSE after Hungarian matching; microstates by "dominant microstate prediction
accuracy" where "The dominant tautomer of each charge state was determined as the microstate with
the lowest free energy" (Sec. 2.3.2).

## Findings
"Top-performing submissions for macroscopic pKa predictions achieved RMSE of 0.7–1.0 pKa units"
(Abstract). "Prediction of dominant tautomers was a major source of error for microscopic pKa
predictions, especially errors in charged tautomers" (Abstract). Successful QM entries credited
"capturing the chemical ensemble for each molecule including conformers and tautomers" (Sec. 3.1).

## Relevance to p(z,x|G,h)
This is the only blind, experimentally-anchored benchmark near the target task, and it supervises
B only weakly (8 molecules with dominant-microstate labels) and C not at all (no measured
state populations, no conformer-resolved data).
