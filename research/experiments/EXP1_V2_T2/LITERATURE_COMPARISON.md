# Literature comparison: refining generated conformers (EXP1-V2-T2)

Written 2026-11-02 09:17 UTC / 2026-11-02 14:47 IST. Sources were checked the same day.

**Our method (FlexiTors V2), for the comparison row only.** It takes conformers sampled by a frozen Torsional Diffusion
model (TD, `qm9_default`) and refines only the local structure (bond lengths, angles, ring puckers). A learned per-atom
equivariant displacement head sits on TD's encoder. Each step is projected by P_W = I − Cᵀ(CCᵀ)⁺C. The rows of C are
the dihedral gradients of one dihedral per rotatable bond plus 6 rigid-motion rows, so torsions are kept to first order.
They are then restored exactly. Training target: the GEOM-QM9 reference conformer with its torsions set to the RDKit
seed's ("pure local-structure target"). Loss: projected endpoint MSE. Physics baselines on the same TD samples:
GFN2-xTB with torsions harmonically restrained to TD's values (then restored), and MMFF94s. Data: GEOM-QM9 only.
Metrics: COV/MAT recall and precision on the E1 subset of 200 test molecules.

**Access labels.** "full text" = paper body read (arXiv HTML). "code only" = official repository read, paper body
not read. "abstract only" = abstract or listing page only. "not stated" = the source we read does not say.
"inaccessible" = we could not open the source.

---

## 1. Torsional Diffusion (TD)
Jing, B., Corso, G., Chang, J., Barzilay, R., Jaakkola, T. "Torsional Diffusion for Molecular Conformer Generation".
NeurIPS 2022. arXiv:2206.01729 (v2; metadata checked with lookup_reference). Code: github.com/gcorso/torsional-diffusion.
**Access: full text** (arXiv v2 HTML, all of it, including App. H). Local notes: `paper_notes/tordiff_jing2022.md`.

1. **Input.** RDKit conformers provide the local structure. Sec. 3: "we use RDKit to provide approximate samples from
   p_G(L)". Algorithm 3 starts by generating RDKit conformers. TD is a generator, not a refiner of an outside sample.
2. **Prediction.** Scores on torsion angles only: one scalar per rotatable bond, with diffusion on the hypertorus (Sec. 3).
3. **Torsions preserved?** No. Torsions are the only thing TD changes. Local structure L stays at RDKit values.
4. **All coordinates?** No. Bond lengths, angles and ring conformations are frozen at RDKit values. Sec. 5 lists
   relaxing "the rigid local structure assumption" as future work.
5. **Objective and targets.** Denoising score matching on the torus. "Conformer matching" (Sec. 4.1) builds the
   targets: RDKit generates L̂, and the torsions are set "to minimize RMSD(C, Ĉ)" to each ground-truth conformer C.
6. **Physics comparison on the same samples.** For ensemble properties (Sec. 4.4, Table 3), generated conformers are
   relaxed with unconstrained GFN2-xTB: a random 100-molecule GEOM-DRUGS subset, min(2K, 32) conformers per molecule.
   App. H, Table 10 gives TD's errors **without** relaxation: median E 36.91 kcal/mol, μ 0.92 D, Δε 4.93 kcal/mol,
   E_min 36.94 kcal/mol. The authors say these errors are "far too large ... to be chemically useful". The text we read
   gives no xTB settings (crude/tight). TD does not compare a learned refiner with xTB. No constrained xTB was found.
7. **Data, metrics, QM9 headline.** GEOM-DRUGS, GEOM-QM9, GEOM-XL. Metrics: COV and AMR, recall and precision.
   2K samples per molecule. **GEOM-QM9 (App. H, Table 7, δ = 0.5 Å, 1,000 test molecules):** COV-R 92.8/100.0,
   AMR-R 0.178/0.147, COV-P 92.7/100.0, AMR-P 0.221/0.195 (mean/median). The RMSD metrics use no relaxation.

## 2. DiSCO
Lee, Danyeong; Lee, Dohoon; Bang, Dongmin; Kim, Sun. "DiSCO: Diffusion Schrödinger Bridge for Molecular Conformer
Optimization". AAAI 2024, vol. 38(12), pp. 13365–13373. doi:10.1609/aaai.v38i12.29238. Verified via OpenAlex
(lookup_reference) and mlanthology. No arXiv version was found. Code: github.com/DanyeongLee/DiSCO (branches `main` and
`reproduce`; the URL given as Danyeong-Lee/DiSCO resolves to the same project).
**Access: abstract only for the paper** (the AAAI page and PDF return a BotStopper "Access Denied" page; OpenReview
does not host it). **Code only** for everything below. We did not see the paper's tables, so we report **no headline
numbers**.

1. **Input.** Conformers from other generators. Abstract (via mlanthology/Semantic Scholar): DiSCO refines conformers
   from any generation method. The `reproduce` README gives four starting sources: RDKit ETKDG, Torsional Diffusion,
   DMCG, and RDKit + clustering. Configs `230709-fromtorsion-qm9-7steps-final.yaml` and `dataset/qm9-torsion.yaml` read
   TD samples from `../torsional-diffusion/geodiff_split_qm9_run/{train,val}_qm9_20steps_50confs.pkl` and
   `test_qm9_20steps.pkl`. So a TD-specific DiSCO is **trained on TD outputs** for the training molecules. The deployed
   `main` model (`fromrdkit-qm9.pt`) starts from raw ETKDG (`dataset.py`: `EmbedMolecule(mol, ETKDGv3())`, MMFF off).
2. **Prediction.** Cartesian coordinates for all atoms. It is an I2SB-style diffusion bridge (`diffusions/i2sb.py`) in
   SE(3)-equivariant Euclidean space (abstract). The network predicts a scaled displacement, `(x_t − aligned_x0)/std_fwd_t`.
   Sampling starts at `pos1`, the input conformer, and runs 7 steps (from-TD) or 10 steps (from-RDKit).
3. **Torsions preserved?** Not by construction (code). We found no torsion constraint or projection.
   `dataset_torsion.py` pairs a random TD conformer of a molecule (`random.choice(torsion_rdmol_list)`) with the stored
   reference record, after rigid alignment only. It does no torsion matching, so the target can differ from the input in
   torsions. The paper's own statement on torsions: inaccessible.
4. **All coordinates?** Yes. The bridge acts on every atom's position, with no masking. `local_structure_error.py` compares
   TD vs TD+DiSCO against references on all-pairs distances and index-ordered triples. That points to a local-structure
   analysis in the paper, but we cannot cite it.
5. **Objective and targets.** MSE between the prediction and the bridge target (`F.mse_loss(pred, target)`). The target
   endpoint x0 is a dataset reference conformer, superposed on the input by `align`. Which reference conformer is
   paired with which input: per stored record (`train_data_40k.pkl`), with no per-pair matching found. The source of the
   references is not named in the code (GEOM, per the README's data links).
6. **Physics comparison on the same samples (code).** (a) MMFF: `scripts/mmff_iters_performance.sh` applies MMFF94s
   (`MMFFOptimizeMolecule(..., mmffVariant='MMFF94s', maxIters=k)`, k = 5 … 1000) to TD's DRUGS samples, then runs
   DiSCO. With `test_starting_point=true` it scores the MMFF-relaxed TD samples on their own. So "TD + MMFF(k)" and
   "TD + MMFF(k) + DiSCO" are evaluated on the same samples. `mmff_four_models_performance` does the same with 200
   iterations for all four sources. The MMFF here is **unconstrained** (no torsion restraint in the code). (b) xTB:
   `xtb_prop.py` relaxes generated and starting-point conformers with xTB (`xtb_optimize(conf_mol, level='crude')`),
   then computes E, gap and dipole. References get single points only. `traj_energy_xtb.py` gives xTB energies along
   the DiSCO trajectory. In this code xTB is used for property evaluation; it is not a competing refiner. **No
   constrained xTB or constrained force field was found.** How the paper presents the MMFF results: inaccessible.
7. **Data, metrics, QM9.** GEOM-QM9 and GEOM-DRUGS. The QM9 test file is `test_data_200.pkl` (200 molecules) on the
   GeoDiff split (path `geodiff_split_qm9_run`). The training file is `train_data_40k.pkl`. Metrics (`test.py`,
   `evaluation.py`): COV and MAT, recall and precision, mean and median, δ = 0.5 Å (QM9) and 1.25 Å (DRUGS), via
   `get_best_rmsd`. Ensemble properties use xTB. **Headline QM9 numbers: inaccessible.** We did not establish whether
   their 200 molecules overlap our E1 subset.

## 3. Flow-Matching Based Refiner (Xu & Gao 2025)
Xu, Xiangyang; Gao, Hongyang. "Flow-Matching Based Refiner for Molecular Conformer Generation". arXiv:2510.04878 v1,
6 Oct 2025 (15:00 UTC / 20:30 IST). No venue on arXiv. OpenAlex W4414972626.
**Access: full text** (arXiv v1 HTML).

1. **Input.** At inference, samples from upstream generators: MCF, ET-Flow and DMT (Sec. 5.1). "Given an upstream
   sample, we set x₀ = x̃" (Algorithm 2 / Sec. 4). TD and RDKit appear only as baselines, never as refiner inputs.
   **Training never uses upstream samples.** Algorithm 1 samples x₁ ∼ p_data and builds x₀ = x₁ + σε with σ = 1
   (Eq. 4), a base that "retains the ground-truth signal".
2. **Prediction.** A time-conditioned Cartesian vector field u_θ(x_t, t, G) on all atoms. It is SE(3)-equivariant and
   fine-tuned from the open-source ET-Flow architecture and weights (Sec. 5.1; lr 0.007 → 0.001).
3. **Torsions preserved?** Not stated. Nothing in the method fixes or projects torsions.
4. **All coordinates?** Yes. The ODE acts on x ∈ ℝ^{3N} (Eq. 10, Alg. 1). No masking is described.
5. **Objective and targets.** Flow-matching regression (Eq. 9) on the interpolant x_t = (1−t)x₀ + t x₁ + s(t)z,
   s(t) = √(t(1−t)) (Eq. 6). x₀ is Kabsch-aligned to x₁, and the target is u_t = x₁ − x₀ (Alg. 1). The target is the
   GEOM ground-truth conformer. Integration starts at a calibrated time t* = 1 − σ*/σ (Eq. 5) to skip the low-SNR
   phase. 20 generation + 20 refinement steps, against 50 steps for the baselines.
6. **Physics comparison on the same samples.** No. GFN2-xTB appears only in Sec. 5.3 (Table 3): min(2K, 32)
   conformers per molecule on "the same 100-molecule subset" as Jing et al. (GEOM-DRUGS) are relaxed before computing
   Boltzmann properties. Both the upstream and the "+ Refiner" ensembles are relaxed, so xTB is part of the evaluation,
   not a baseline. Table 3 (median): MCF-L EE 0.68 → 0.20 kcal/mol with the refiner; ET-Flow 0.23 → 0.21; ET-Flow Δε
   0.38 → 0.39. No MMFF or DFT relaxation of the same samples was reported.
7. **Data, metrics, QM9.** GEOM-QM9 and GEOM-DRUGS, Ganea et al. splits (106,586 / 13,323 / 1,000 for QM9). 2K samples
   per molecule. Recall and precision COV/AMR. **QM9 Table 2 uses δ = 0.05 Å, not the usual 0.5 Å.** The caption calls
   it "the more challenging COV threshold" because 0.5 Å is saturated. So these COV values are **not comparable** to TD's
   Table 7 or to ours. Table 2 (mean/median): ET-Flow COV-R 75.72/87.23, AMR-R 0.083/0.031, COV-P 70.32/75.00,
   AMR-P 0.114/0.053. ET-Flow + Refiner: 78.40/88.89, 0.076/0.028, 77.36/89.94, 0.103/0.031. DMT-B + Refiner:
   79.50/89.44, 0.070/0.026, 80.37/97.92, 0.076/0.021. DRUGS: Table 1, δ = 0.75 Å.

## 4. CoarsenConf
Reidenbach, Danny; Krishnapriyan, Aditi S. "CoarsenConf: Equivariant Coarsening with Aggregated Attention for Molecular
Conformer Generation". arXiv:2306.14852 (v2, 19 Oct 2023). Shown at NeurIPS 2023 as a **workshop** poster
(GenBio@NeurIPS2023, per neurips.cc/virtual/2023/74750). Code: github.com/ASK-Berkeley/CoarsenConf.
**Access: full text** (arXiv v2 HTML, all of it).

1. **Input.** An RDKit approximate conformer R: ETKDG + MMFF when MMFF converges, otherwise plain ETKDG (App. G). The
   learned conditional prior "only uses the RDKit" structure (Sec. 3). Training also encodes the ground truth X.
2. **Prediction.** Cartesian coordinates for all atoms, from a hierarchical SE(3)-equivariant VAE. Coarse beads are
   the fragments left after cutting rotatable bonds. The decoder predicts the offset, X = R + ΔX ("predict the difference
   between the RDKit reference and ground truth conformations", Sec. 3).
3. **Torsions preserved?** No. The model "can adjust atomic coordinates, distances, and torsion angles freely"
   (Conclusion). It "leave[s] explicit angle optimization to future work" (App. B).
4. **All coordinates?** Yes.
5. **Objective and targets.** ELBO (Eq. 1). Implemented as Kabsch-aligned MSE + β₁·KL + β₂·distance loss over 1- and
   2-hop edges (App. B, Eq. 5). Targets are GEOM ground-truth conformers, 5 per molecule. "OT" adds an optimal-transport
   assignment between generated and ground-truth conformers. We found no torsion-matched pairing of RDKit inputs to
   targets.
6. **Physics comparison on the same samples.** Not as a refiner baseline. For properties (Sec. 4.2, Table 2), generated
   conformers are relaxed with GFN2-xTB (BFGS) before properties are computed. Rigid docking minimizes all poses with
   Vina BFGS (Sec. 4.4). Flexible docking uses RDKit+MMFF as one baseline (Sec. 4.3). We did not see a comparison of
   CoarsenConf against xTB or MMFF relaxation of its own inputs.
7. **Data, metrics, QM9.** GEOM-QM9 / GEOM-DRUGS, splits 106,586 / 13,323 / 1,000 (Sec. 4.1). **QM9 Table 1 reports
   precision only, at δ = 0.5 Å:** CoarsenConf-OT COV-P 80.2/100.0, AMR-P 0.149/0.107. TD (ℓ=2) is listed in the same
   table at 78.4/100.0, 0.222/0.197, from TD's public samples. TD's own Table 7 gives COV-P 92.7, but AMR-P 0.221
   matches. We did not resolve the difference in COV-P.

## 5. Uni-Mol (conformation generation)
Zhou, G., Gao, Z., Ding, Q., Zheng, H., Xu, H., Wei, Z., Zhang, L., Ke, G. "Uni-Mol: A Universal 3D Molecular
Representation Learning Framework". ICLR 2023 (mlanthology bibtex; OpenAlex also lists ChemRxiv
doi:10.26434/chemrxiv-2022-jjm0j-v4). Code: github.com/deepmodeling/Uni-Mol (formerly dptech-corp).
**Access: abstract only for the paper** (OpenReview forum/PDF behind a browser check; ChemRxiv 403). **Code only**
for 1–7 (`unimol/utils/conf_gen_cal_metrics.py`, `unimol/losses/conf_gen.py`, README).

1. **Input.** RDKit conformers. Per molecule with K references, `clustering()` pools M = min(20K, 2000) ETKDG+MMFF
   conformers plus M/4 plain ETKDG conformers. It keeps N = 2K K-medoids centres as inputs. Hydrogens are removed.
2. **Prediction.** The network outputs both a pairwise distance matrix and Cartesian coordinates
   (`distance_predict, coord_predict = net_output[0], net_output[-1]`). Metrics read the coordinates directly. We found
   no distance-to-coordinate fit in the metric script.
3. **Torsions preserved?** Not stated, and no constraint in the code.
4. **All coordinates?** Yes (heavy atoms; RemoveHs).
5. **Objective and targets.** L1 coordinate loss after rotational alignment, plus L1 distance loss (`conf_gen.py`;
   README weights coord_loss = 1, distance_loss = 1). In `gen_data` mode every input record of a SMILES gets the same
   target, the first reference conformer. We did not resolve how training pairs are built.
6. **Physics comparison on the same samples.** Not found. The metric script has an optional `use_ff` (MMFF on outputs
   before RMSD), default False.
7. **Data, metrics, QM9.** GEOM-QM9 / GEOM-DRUGS. The QM9 reference file is `test_data_200.pkl` (200 molecules). COV
   and MAT, recall only in this code, heavy-atom `GetBestRMS`, δ = 0.5 Å (QM9) and 1.25 Å (DRUGS). **Headline numbers:
   inaccessible** (paper body not read).

## 6. Optional works found
- **RINGER** (Grambow et al., arXiv:2305.19800, macrocyclic peptides; metadata checked). **Abstract and search snippet
  only.** According to the search result, all samples from every method (RDKit, OMEGA, GeoDiff, TD, RINGER) were
  re-optimized with GFN2-xTB "to enable a level comparison". This is unconstrained relaxation of the same samples, used
  as an equaliser rather than as a refiner baseline. We did not verify the paper's sections.
- **GeoOpt-Net** (Liu, Cai, Zhao, Mo, arXiv:2601.22723, 30 Jan 2026). **Abstract only.** A single-pass SE(3)-equivariant
  refiner from force-field conformers to B3LYP/TZVP geometries. It is benchmarked against RDKit, xTB, Auto3D and UMA.
  Torsion handling: not stated. It does not use GEOM-QM9 generator outputs.
- **LoQI / ChEMBL-3D** (Nikitin et al., NeurIPS 2025): a generator trained on AIMNet2-optimized geometries, not a
  refiner. Listed only because the search returned it.

---

## Comparison table

| Work (access) | (1) Input | (2) Predicts | (3) Torsions kept | (4) All coords | (5) Objective / target | (6) Physics on same samples | (7) Data / metrics / QM9 headline |
|---|---|---|---|---|---|---|---|
| TD (full text) | RDKit local structure | torsion scores | no (only torsions move) | no (L frozen) | torus DSM; RDKit L + torsions RMSD-matched to GT | unconstrained GFN2-xTB relaxation for properties only (DRUGS, 100 mols); raw errors Table 10 | QM9/DRUGS/XL; COV/AMR R+P; QM9 δ=0.5: COV-R 92.8, AMR-R 0.178, COV-P 92.7, AMR-P 0.221 (Table 7) |
| DiSCO (abstract + code) | RDKit, TD, DMCG, RDKit-clustering samples; TD-specific model trained on TD samples | Cartesian (bridge displacement) | no constraint in code | yes | MSE to bridge target; endpoint = rigidly aligned dataset reference, no torsion matching | MMFF94s (unconstrained, 5–1000 iters) on TD samples, then DiSCO; xTB crude opt for properties; constrained: not found | QM9 (GeoDiff split, 200 test) + DRUGS; COV/MAT R+P δ=0.5/1.25; numbers inaccessible |
| Xu & Gao 2025 (full text) | MCF / ET-Flow / DMT samples (inference only); training on GT + Gaussian noise | Cartesian vector field | not stated / no constraint | yes | flow matching, Kabsch-aligned, target = GEOM GT | xTB relaxation of both upstream and refined sets for properties only; no FF/xTB baseline | QM9 δ=0.05 Å (non-standard); ET-Flow+Refiner COV-R 78.40, COV-P 77.36 (Table 2) |
| CoarsenConf (full text) | RDKit ETKDG(+MMFF) conformer | Cartesian offset ΔX (VAE) | no ("freely") | yes | ELBO: aligned MSE + KL + distance; GT conformers, OT variant | GFN2-xTB relaxation for properties; MMFF only as a docking baseline | QM9 δ=0.5, precision only: COV-P 80.2, AMR-P 0.149 (Table 1) |
| Uni-Mol (code only) | K-medoids of ETKDG(+MMFF) conformers, heavy atoms | distances + Cartesian | no constraint in code | yes (heavy atoms) | L1 coord (aligned) + L1 distance | not found (optional MMFF on outputs, off) | QM9 200-mol test file; COV/MAT recall δ=0.5; numbers inaccessible |
| FlexiTors V2 (ours) | frozen TD (qm9_default) samples | per-atom Cartesian displacement, projected | yes: first order via P_W, then restored exactly | no: torsions + rigid motion removed | projected endpoint MSE; GT with torsions set to RDKit seed's | torsion-restrained GFN2-xTB and MMFF94s on the same TD samples | GEOM-QM9 only; COV/MAT R+P; E1 = 200 test mols |

## What this implies for FlexiTors V2 novelty (evidence-limited)
- **Refining TD's samples is not new.** DiSCO's code trains and tests a refiner on TD outputs for GEOM-QM9 and DRUGS.
  It also runs MMFF94s on the same TD samples (alone and before DiSCO). The FlexiTors V2 write-up should cite DiSCO as
  the closest prior work and compare with it where possible. DiSCO's own numbers could not be read; its code and TD
  sample archives are public.
- **Torsion-preserving refinement.** None of the five works or the searches below showed a refiner that keeps torsions
  fixed by construction (projection or restraint) while changing only local structure. DiSCO, Xu & Gao, CoarsenConf and
  Uni-Mol move all coordinates. TD keeps local structure fixed and moves only torsions. We say "not found in the works
  and searches listed", not "absent from the literature".
- **Torsion-matched local-structure target.** TD's conformer matching sets RDKit torsions to match the ground truth.
  Our target does the reverse (ground-truth local structure with the RDKit seed's torsions). We found no work that uses
  this target for a refiner.
- **Constrained physics baseline.** Xu & Gao, CoarsenConf and TD use xTB only to relax conformers before computing
  properties. DiSCO's code uses unconstrained MMFF on the same samples. No torsion-restrained xTB or force-field baseline
  for a learned refiner was found in the works and searches listed.
- **Comparability caveats.** Xu & Gao's QM9 COV uses δ = 0.05 Å. CoarsenConf's QM9 table is precision only and gives TD
  a different COV-P than TD's paper. DiSCO and Uni-Mol use 200-molecule QM9 test files. Our E1 averages over paired
  molecules (184 of 200). Do not put headline numbers from these works next to E1 numbers without re-evaluating under
  one protocol.

## Searches run
- WebSearch: "DiSCO Diffusion Schrödinger Bridge Molecular Conformer Optimization Lee pdf"
- WebSearch: "\"DiSCO\" \"Schrödinger bridge\" conformer \"torsional diffusion\" MMFF \"local structure\""
- WebSearch: "learned refinement of generated molecular conformers compared to GFN2-xTB optimization same samples torsional diffusion local structure"
- WebSearch: "conformer refinement preserve torsion angles refine bond lengths angles neural network \"torsional diffusion\" local structure correction"
- WebSearch: "machine learning conformer refinement baseline \"constrained\" dihedral xTB optimization generated conformers comparison"
- WebSearch: "Uni-Mol universal 3D molecular representation learning chemrxiv conformation generation RDKit distance matrix"
- WebSearch: Nikitin/Isayev "Scalable Low-Energy Molecular Conformer Generation with Quantum Mechanical Accuracy"
- OpenAlex search_works: "molecular conformer refinement diffusion RDKit conformers optimization" (15 results)
- OpenAlex citation_neighbours of DiSCO (6 citing works; none is a conformer refiner)
- lookup/verify_reference: 10.1609/aaai.v38i12.29238, 2206.01729, 2306.14852, 2305.19800, Uni-Mol title

Findings from these searches: RINGER (xTB re-optimization of all methods' samples), GeoOpt-Net (learned refiner
benchmarked against xTB on force-field inputs) and ParametrizANI (constrained dihedral scans, not a refiner). None
compared a learned refiner against **torsion-constrained** xTB or force-field relaxation of the same generated samples.

## Inaccessible sources
- DiSCO paper body: ojs.aaai.org article page and PDF (`/article/download/29238/30337`) returned a BotStopper
  "Access Denied" page. scholar.archive.org needed a browser check. Underline has only the title and event details. No
  arXiv version was found. Everything about DiSCO above comes from the abstract and the GitHub code.
- Uni-Mol paper body: OpenReview forum and PDF (browser check); ChemRxiv (HTTP 403).
- RINGER, GeoOpt-Net, LoQI: abstract or search snippet only (not pursued further; optional).
