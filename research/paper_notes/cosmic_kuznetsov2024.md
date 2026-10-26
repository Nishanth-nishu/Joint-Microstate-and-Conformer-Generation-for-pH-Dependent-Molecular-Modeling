# COSMIC (Kuznetsov et al. 2024)

## Citation
COSMIC: Molecular Conformation Space Modeling in Internal Coordinates with an Adversarial Framework.
Maksim Kuznetsov, Fedor Ryabov, Roman Schutski, Rim Shayakhmetov, Yen-Chu Lin, Alex Aliper, Daniil Polykovskiy
(Insilico Medicine). J. Chem. Inf. Model. 2024. DOI 10.1021/acs.jcim.3c00989. PMID 38668753. Checked 2026-10-26.

## Fulltext status
not_accessible. pubs.acs.org returned HTTP 403; no arXiv/ChemRxiv preprint was found. Only title, venue, authors
and third-party summaries were obtained. Essentially all mechanism-level questions are Unknown.

## What is known (metadata + third-party description, flagged as such)
Generative model of conformation space in internal coordinates using an adversarial framework (GAN + adversarial
autoencoder, consistent with the related patent application US 2022/0406404). Described second-hand as a two-step
procedure: build a conformation from internal coordinates, then predict distances between neighbouring atoms to
refine it. An energy-based evaluation metric named RED is introduced. No code repository was located.

## Inputs/outputs — protonation / charge / H handling
Unknown. What would resolve it: the JCIM methods section / Supporting Information (institutional access), or a code
release.

## Relevance to p(z, x | G, h)
Relevant only as the internal-coordinate (non-diffusion) point of comparison, and as the one method in this family
that proposes an energy-based ensemble metric. Not auditable from open sources at present.
