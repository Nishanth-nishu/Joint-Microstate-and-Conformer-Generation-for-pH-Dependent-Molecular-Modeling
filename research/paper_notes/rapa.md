# RAPA

## Citation
- Title: A Self-Consistent Approach to Rotamer and Protonation State Assignments (RAPA): Moving Beyond Single Protein Configurations
- Authors: Mossa Ghattas, Prerna Gera, Steven Ramsey, Anthony Cruz-Balberdy, Nathan Abraham, Vjay Molino, Daniel McKay, Tom Kurtzman
- Year / venue: 2025, J. Chem. Inf. Model. 65(14):7639–7650 (received 2025-04-15, accepted 2025-06-03); earlier ChemRxiv preprint (article 675c85b2085116a13353ce8b)
- DOI: 10.1021/acs.jcim.5c00859; PMC12308802 (CC-BY 4.0)
- Code: https://github.com/Ventus-Therapeutics/rapa/tree/v1.0 (MIT); data https://github.com/KurtzmanLab/RAPAData
- Date checked: 2026-10-26 08:26 UTC (13:56 IST)

## Fulltext status
full_text_read (PMC HTML via WebFetch + Europe PMC fullTextXML, saved at research/papers/rapa_europepmc.xml). Supporting Information not read. Code not read.

## Problem
X-ray structures (≤1.5 Å) do not resolve H and cannot distinguish O/N/C, so side-chain rotamer (Asn/Gln amide flips, His ring flips) and protonation states (His HID/HIE/HIP, Asp/Glu in acid dyads, hydroxyl H direction) are ambiguous. Existing prep tools return one configuration; RAPA returns *all* configurations energetically consistent with the crystal structure.

## Method (Sections 2.2.1–2.2.4)
- Step 1: mark unambiguous polar atoms known; residues with unknown polar heavy atoms or alternate rotamers are "ambiguous".
- Step 3: for each ambiguous residue, enumerate its RP (rotamer+protonation) states jointly with all RP states of proximal neighbours (polar heavy atoms ≤3.8 Å); fix the residue if one state is best.
- Iterate to self-consistency (Steps 2–4).
- Step 5 branching: states within a user cutoff (1 kcal/mol here) of the lowest are "degenerate"; the system is copied once per degenerate state; each branch recursed (Step 6) and written as PDB (Step 7).
- Residue rules (2.2.2): His 6 RP states (2 rotamers × HIE/HID/HIP); Asn/Gln 2 rotamers; Asp/Glu deprotonated unless in an acid dyad (carboxylate O ≤3.8 Å), then all singly protonated configurations; Ser/Thr 9 H positions; Tyr 2 planar H positions; Lys/Arg always charged; Cys by disulfide distance.
- Energy (2.2.3): H-bond energy lookup table from Jaguar water-dimer PES indexed by heavy-atom distance and angle; donor–donor/acceptor–acceptor clash = negative of H-bond energy; RP state energy = sum of local polar terms.
- Coupling of protonation with structural alternatives: protonation and side-chain flip are a single discrete RP variable per residue, scored jointly with neighbours' RP variables. These are *not* PDB altloc records (no mention of altloc/occupancy in the text); the heavy-atom skeleton is the input crystal structure (flip = swap of N/O labels or ring orientation).

## Data
77 DUD-E proteins without missing loops, covalent ligands or heme; prepared with OpenEye Spruce; ligands and waters removed (2.1). MD validation: Amber ff99SB + OPC, 10 ns, restrained (3DS) complex simulations and unrestrained rotamer-state simulations.

## Inputs/outputs
Input: one prepared protein PDB. Output: a set of full-protein PDB configurations, each a fixed assignment of RP states. No pH input; no new heavy-atom coordinates; no probabilities/weights attached to configurations.

## Experiments & results
- 84% Asn/Gln and 89% His fixed; 4% Asn/Gln and 10% His degenerate (Table 1).
- 63/77 systems with >1 configuration; 469 configurations; 62/77 have ≤8; max 24 per monomer.
- 467/469 configurations stable (RMSD-90% ≤2 Å) in 10 ns MD; 56% of degenerate Asn/Gln and 92% of degenerate His had ≥2 stable states in unrestrained MD.
- No baseline comparison against other prep tools' accuracy in the main text (other tools only cited in Introduction).

## Limitations (authors')
Simple residue–residue energetics; no pKa corrections; no entropy → no free energy of relative populations; limited hydroxyl/thiol rotational treatment; 44% of degenerate Asn/Gln had only one stable state in MD; prioritisation of configurations unresolved.

## What transfers to small molecules
- The idea of returning a *set* of state assignments within an energy window instead of one, and of coupling a discrete protonation/tautomer-like choice with a discrete local geometric alternative (e.g. amide/ring flip ↔ H position) under one score.
- Self-consistent / branching enumeration as a cheap alternative to full combinatorics, and the empirical finding that the accessible set is small.
- What does not transfer: residue-template rules (His/Asn/Gln/Asp/Glu/Ser/Thr/Tyr), dependence on a fixed crystal heavy-atom scaffold and on inter-residue H-bond networks (intramolecular-only for a free ligand), and the protein-specific hard rules (Lys/Arg always charged).

## Relevance to p(z,x|G,h)
RAPA is a precedent that protonation state and local structure should be assigned jointly and that a single assignment is insufficient. It does **not** establish: pH dependence (no pH), population weights (explicitly missing entropy/free energy), generative 3D sampling (only enumerates discrete H/flip choices on a fixed scaffold), or anything about small-molecule conformational ensembles. It is evidence for the motivation, not overlap with a generative model.

## Exact supporting quotes
- Abstract: "we posit that there are multiple rotamer and protonation states that are consistent with the resolved structure of the proteins"
- 2.2.1.5: "We term these as “degenerate” states which consist of the RP state with the lowest assessed energetic environment … and all other RP states that have an environmental energy within a user defined cutoff"
- 2.2.2.2: "Local energetics are evaluated for the six RP states which are all combinations of the 2 possible rotamer states and 3 protonation states (HIE, HID, and HIP)."
- 2.2.2.4: "They are considered negatively charged (unprotonated) unless they are part of an acid dyad."
- 2.2.2.1: "Lysine (LYS) and Arginine (ARG) side chain nitrogen atoms are always assigned the positively charged protonation state."
- Discussion: "calculating residue specific corrections based on environmental estimates of pKa, inclusion of an entropy term to obtain a free energy estimate of the relative populations of alternate RAPA states"
- Footnote 1: "We use the term “configuration” to denote a complete description of the positions of every atom in the protein."
