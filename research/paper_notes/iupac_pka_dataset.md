# IUPAC Digitized pKa Dataset v2.4a (direct inspection)

## Citation
Zheng J. W., Lafontant-Joseph O., IUPAC Digitized pKa Dataset v2.4a, Zenodo
DOI 10.5281/zenodo.7236452; paper: Zheng, Lafontant-Joseph, Green, "Digitized Dataset of Aqueous
Acid Dissociation Constants", RSC Advances 16(23), 2026. Repo: github.com/IUPAC/Dissociation-Constants.
**License CC BY-NC 4.0** ("Reproduced by permission of International Union of Pure and Applied
Chemistry"). Checked 2026-10-26; CSV copied to /home2/nishanth.r/jmcg_stage/ext/iupac_pka/.

## Fulltext status
full_text_read (README) + the data file itself profiled (6.9 MB, 24,211 rows).

## Provenance
Hand digitization of three IUPAC compilations: Serjeant & Dempsey (organic acids, 1979), Perrin
(organic bases, 1965) and its 1972 supplement. The README calls it an early release with possible
remaining errors.

## Profile (measured, not quoted)
- 24,211 rows / 10,626 unique molecules (by InChI).
- `pka_type`: pKaH1 11,155; pKa1 8,716; pKa2 1,674; pKaH2 1,134; pKa3 390; pKa4 174; pKb1 158;
  pKaH3 93; pKb 89; pKa5 46 … → **macroscopic, index-ordered constants. No microstate identity,
  no state-specific graph.** (A few `pka_type` strings name a site, e.g. "pK(indole-ring)".)
- `acidity_label`: AH 12,602; A 11,141; B 302; other 166.
- `T` (°C): 25 → 8,351; 20 → 6,524; blank → 2,085; not_stated → 1,366; 30 → 905; 40 → 509.
- `cosolvent` non-blank for 2,786 rows (so ~11% are not pure water).
- 2,957 rows have a charged SMILES.
- Quality flag `assessment` ∈ {Reliable, Approximate, Uncertain}; `method` codes are translated by
  `method_translation.csv`.

## Relevance to p(z,x|G,h)
The best-curated **experimental** supervision available for the z-marginal (task B), with
temperature, cosolvent and reliability metadata — but zero 3D, zero microstate resolution, and a
**non-commercial** license that must be respected in any released benchmark.
