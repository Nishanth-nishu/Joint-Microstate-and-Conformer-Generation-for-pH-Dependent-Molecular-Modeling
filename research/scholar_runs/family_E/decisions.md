# family_E decisions (2026-10-26 UTC)

- Google Scholar: not searched, no tool (queries.csv row `not_searched_no_tool`).
- No multi-GB file was downloaded. Dataset facts come from (a) the Harvard Dataverse REST API
  (file list, sizes, license), (b) the GEOM README (4.5 kB) and tutorial notebook 01, (c) PMC
  full texts, (d) GitHub REST listings and small raw files. Small samples kept under
  /home2/nishanth.r/jmcg_stage/ext/: geom_readme/, iupac_pka/ (6.9 MB CSV), sampl6_pka/, sampl8/.
- Nature and ACS pages 403/redirect to login for WebFetch → used PMC mirrors throughout.
- Harvard Dataverse HTML page returned empty to WebFetch → backtracked to its JSON API, which gave
  the authoritative license (CC0 1.0) and file sizes.
- OMol25 and SPICE could not be verified for *protonation-state* content from the pages fetched.
  Recorded as `Unknown` with the resolving action (read arXiv 2505.08762 dataset section /
  SPICE generation scripts), per CONVENTIONS honesty rule — NOT recorded as "no protonation states".
- QMugs was initially a promising conformer+drug-like candidate; dropped to relevance 2 after the
  paper stated "Molecules were neutralized, and salts and solvents were removed" and nonzero-charge
  species were removed. It therefore cannot supervise protonation at all.
- MPCONF196 kept at relevance 1: high-accuracy conformer energies but 13 neutral compounds, no
  protonation variation.
- Neutron-diffraction ligand protonation: searched; found only individual case studies (e.g. 5MOO
  trypsin/aniline) and an undated count of ~88 neutron PDB entries. No curated, protonation-resolved
  ligand dataset found. Recorded as a gap, not as non-existence.
- ChEMBL3D/LoQI examined only to the depth needed for the table (another agent audits it).
