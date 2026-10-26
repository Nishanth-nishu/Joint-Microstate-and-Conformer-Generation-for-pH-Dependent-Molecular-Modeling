# family_B — existing conformer-generation architectures: what they actually assume

Run date 2026-10-26 (UTC). 20 methods noted (cap respected). Companion files: queries.csv, candidates.csv,
feature_matrix.csv, decisions.md; per-paper notes in research/paper_notes/.
Google Scholar was not available and was never queried.

## 1. The headline finding

**No conformer generator in this family conditions on pH, carries a categorical chemical-state variable, or
generates more than one chemical graph.** Every one of them learns p(x | G) for a single, fixed, fully-specified
input graph G. The labels "pH conditioning", "state prediction", "multiple protonation states", "categorical
diffusion", "explicit state transition", "tautomer support", "multi-site protonation" and "unseen-pH evaluation"
are negative or Unknown across the entire matrix. That is the gap the proposed p(z, x | G, h) would fill, and the
audit does not find it pre-empted — but see §4 for why the data, not the architecture, is the binding constraint.

## 2. Charge and hydrogen handling (what an input graph is allowed to be)

Confirmed from methods sections or official code:

| Method | explicit H | formal charge feature | evidence |
|---|---|---|---|
| Torsional Diffusion | yes (`Chem.AddHs`) | yes, bucket **{-1, 0, +1}** + "other" | utils/featurization.py |
| GeoMol | yes (`Chem.AddHs`) | yes, bucket {-1, 0, +1} | model/featurization.py |
| ET-Flow | yes (`Chem.AddHs`) | yes (`GetFormalCharge`) | etflow/commons/featurization.py |
| DiTMC | yes (H in atom-type vocab; num-H 0-4 feature) | yes, **-5 … +5** + other | Table A8 |
| SO(3)-AvgFlow | yes (total-num-H feature) | yes ("charge — formal charge of atom") | appendix feature table |
| MCF | Unknown | yes, {-1 … +1} | Table 6 |
| CoarsenConf | Unknown | "atomic charge" among node features | Sec. on featurization |
| ConfGF | yes (`AddHs`) | **no** — stored node feature is atomic number only | confgf/dataset/dataset.py |
| StoL | yes, explicitly | not addressed; fragments limited to C, N, O, F | Sec. on method |
| Transferable BG | yes (H indistinguishable within a heavy atom) | not mentioned; topology fixed by Amber-14 | Sec. 3/4 |
| GeoDiff, EC-Conf, Uni-Mol, COSMIC, PDCG, MSGEN, GeoMFormer | Unknown | Unknown | not resolved — see decisions.md |

Two practical consequences. (a) The field's default charge encoding is a three-way bucket inherited from GeoMol by
Torsional Diffusion and MCF; a dianion or a doubly protonated base collapses into "other". DiTMC and SO(3)-AvgFlow
are the exceptions with a usable charge range. (b) Formal charge is present as a *feature* everywhere it appears —
never as a variable the model reasons about, predicts or changes.

## 3. Mechanism-level notes worth keeping

- **Torsional Diffusion** factorizes p_G(L, τ) = p̂_G(L)·p_G(τ|L) and takes L (bond lengths, angles, **ring
  conformations, and every X-H geometry**) from RDKit. A protonation change that adds an H to a non-rotatable site
  or re-puckers a ring is therefore not something TD can express: the protomer must be fully fixed before RDKit
  builds L. This is the single most important architectural caveat for using TD downstream of a state model.
- **Torsional Boltzmann generator** (TD Sec. 3.6/4.5) is the only small-molecule result in this family that
  produces calibrated conformer populations: trained with MMFF energies, exact likelihood from the probability-flow
  ODE, samples reweighted to the target density. Cost: "a separate model is trained for every temperature", fixed
  local structures, and only 3-7 rotatable bonds tested.
- **Transferable Boltzmann Generators** (Klein & Noé) give transferable, reweightable equilibrium sampling, but on
  dipeptides with a fixed Amber topology, and they treat a changed bond graph as an invalid sample rather than as a
  second state — the exact opposite of what a protonation-microstate model needs.
- **Particle Guidance** improves TD's precision/coverage (~19% / ~8% median error reduction) by making the sample
  set deliberately non-i.i.d. — so its sample frequencies must not be read as populations.
- **CoarsenConf** and **MSGEN** are the two coarse→fine cascades. Both coarse variables are geometric/topological
  (rigid fragment beads; scaffolds), never chemical, so neither is a latent-state precedent — but MSGEN's "coarse
  stage conditions fine stage, and wraps around existing generators" is the structural slot a sequential z→x
  baseline would occupy.
- **GeoMFormer** is not a generator at all; it is the closest published template for one backbone serving an
  invariant head (classification, e.g. a state) and an equivariant head (coordinates).
- Population/Boltzmann weights: every GEOM-benchmarked generator (TD, ET-Flow, DiTMC, MCF, SO(3)-AvgFlow,
  CoarsenConf, EC-Conf) evaluates Boltzmann-weighted ensemble properties only *post hoc*, by relaxing samples with
  GFN2-xTB and comparing to CREST. None predicts a weight.

## 4. Datasets and splits — the real obstacle

- **Splits are not scaffold-separated anywhere in this family.** TD states it plainly: "The splits are random with
  train/validation/test of 243473/30433/1000 for GEOM-DRUGS". MCF: "We split them randomly based on molecules".
  DiTMC, SO(3)-AvgFlow, CoarsenConf and ET-Flow all reuse the same Ganea et al. split. StoL (small→large size
  extrapolation) and Transferable BG (unseen dipeptides) are the only held-out-chemistry evaluations found, and
  neither is a drug-like scaffold split. Any scaffold-separated claim this project makes will have no comparable
  published baseline number and must be re-run.
- **GEOM contains charged species but exactly one protonation state per compound, with no pH label.** From the GEOM
  methods: 3.9% of drug molecules arrived as salt clusters, and "for the acid/base SMILES strings, used reaction
  SMARTS in RDKit to generate the protonated molecule and counterion"; de-salting "identif[ies] the drug-like
  compound in each cluster and correct[s] its ionization state". Charge then enters physics only as a scalar:
  "Default values were used for all CREST arguments, except for the charge of each geometry."
  Tautomers were deliberately not enumerated: "CREST has specific commands to generate tautomers, but they were not
  used here." The MoleculeNet-derived part of GEOM was de-salted but "did not perform protonation", so the
  convention is not even uniform within GEOM.
  **There is therefore no supervision for z in GEOM.** It would have to be manufactured — enumerate protomers with
  a pKa/microstate predictor and re-run CREST per protomer — and that is a dataset-construction project, not a
  modelling one. GEOM does supply per-conformer energies and statistical weights (with the authors' warning that
  CREST weights are inaccurate and DFT weights should be used for probability benchmarks), which is the one piece a
  population head could train on.
- Most GEOM ensembles are **vacuum**, which sits awkwardly with a pH-conditioned (aqueous) model; only the 1,513
  BACE species have implicit-water CREs.

## 5. Which method should be the conditional geometry generator in a sequential z → x baseline

Ranked, with what each would need to accept charged / variant-H inputs:

1. **DiTMC** — best fit. Formal charge already encoded over -5…+5, explicit H in the atom vocabulary plus a
   num-H feature, full Cartesian output (so ring and X-H geometry respond to the protomer), current SOTA precision,
   code released. Needed changes: none to accept a charged graph; add a conditioning token for z (and h) alongside
   the existing atom-wise/pair-wise tokens — the AdaLN conditioning path is already the natural insertion point.
   Caveat: O(N²) attention; GEOM-only training distribution.
2. **ET-Flow** — same input capability (AddHs + GetFormalCharge confirmed in code), cheaper inference, released
   code, chirality machinery already present. Needed: verify the charge encoding width; add z conditioning.
3. **SO(3)-AvgFlow** — if 1-2 step sampling matters. Formal charge and total-H features confirmed; no code URL was
   found, which is the practical blocker.
4. **Torsional Diffusion** — already reproduced here, cheapest to extend, but structurally the worst fit for a
   protonation pipeline: torsions only, local structure frozen from RDKit, and a {-1,0,+1} charge bucket that would
   need widening. Usable as the *baseline* precisely because it is already running; should not be the end state.
5. **MCF** — graph-agnostic (Laplacian eigenvectors), so a different protomer is simply a different graph with no
   architectural friction; held back by the narrow charge range and unconfirmed H handling.

Not recommended: ConfGF (no charge channel), StoL (neutral C/N/O/F fragments only), GeoMFormer (not a generator),
EC-Conf / COSMIC / PDCG / Uni-Mol / MSGEN (featurization unresolved).

## 6. Open unknowns (and what would resolve each)

- GeoDiff node features — read `utils/datasets.py` / `models/epsnet/*` in MinkaiXu/GeoDiff.
- EC-Conf atom features; the paper only says time + atom-type embeddings reach the network — read the EcConf repo.
- MCF and CoarsenConf hydrogen treatment — their preprocessing code.
- Uni-Mol conformer featurization (explicit H? charge?) — `unimol/data/` and the conformation-generation config.
- COSMIC (ACS 403) and PDCG (MDPI 403) — need institutional access or a different network path.
- MSGEN beyond the abstract — OpenReview PDF (bot-checked from this session).
- Whether ET-Flow or any of these papers reports a scaffold-split ablation anywhere in its appendix.
- Whether GEOM's `uncleaned_smiles` field could be mined to recover the *original* ionization states of the 11,886
  salt-derived drug molecules, which would give a small, real, multi-protomer subset to test against.
