# Phase I — Baseline audit: Torsional Diffusion (TD) as used in this project

Audited 2026-10-26 (UTC). Sources: upstream `gcorso/torsional-diffusion` @ `5f713b4` (read-only clone), the
project checkout on gnode118 `/scratch/nishanth.r/tordiff/torsional-diffusion` (upstream + `patches/0001-ablation-hooks.patch`,
HEAD `dd354f5`, all hooks default-off), and prior-work docs carried over in this repo (`prior_work/`, `notes/`,
`RESULTS_QM9.md`, `EXPERIMENTS.md`). **Nothing in the baseline was modified.**

## 1. Status of the baseline (prior work, not redone)

| Item | Value | Evidence |
|---|---|---|
| Reproduction | AMR-R 0.176 Å (paper 0.178), COV-R@0.5 88.6 % (paper 92.8 %; gap = 65 RDKit embedding failures counted as 0) | `RESULTS_QM9.md` §1–2 |
| Own retrain | CTRL s0 0.179 Å, 3 training seeds (`workdir/qm9_CTRL_base_100ep_e100_s{0,1,2}`) | `RESULTS_QM9.md`, job-0 inventory |
| Released checkpoint | `workdir/pretrained/` (qm9_default: 44 node feats, 2nd-order irreps) | `slurm/common.sh` `QM9_MODEL_ARGS` |
| Env | Python 3.9.23 venv, torch 1.13.1+cu117, PyG 2.0.4, e3nn 0.5.1, RDKit 2022.9.5, spyrmsd 0.5.2 | Phase-0 job 40060 log |
| Data | `data/QM9` 14 GB (GEOM-QM9 pickles, standardized pickles, split.npy) | Phase-0 inventory |
| Cost | ≈ 10.7–11 h per 100-epoch QM9 training run on one RTX 3090 | `RESULTS_QM9.md` §4 |
| GEOM-DRUGS | `drugs.tar.gz` truncated (18.35/25.31 GB) — not usable yet | `prior_work/TORDIFF_NAX_HANDOFF.md` §6 |

## 2. Repository map (TD, files that matter)

| File | Role |
|---|---|
| `utils/featurization.py` | RDKit Mol → PyG `Data`: node features x, `edge_index`, `edge_attr` (bond type one-hot), `z` |
| `utils/dataset.py` | `ConformerDataset`: GEOM pickles → one `Data` per molecule with list of conformer `pos` + Boltzmann weights; `TorsionNoiseTransform` adds wrapped-normal torsion noise |
| `utils/torsion.py` | `get_transformation_mask` (rotatable bonds + which side rotates), `modify_conformer`/`perturb_batch` (apply Δτ) |
| `diffusion/torus.py` | Wrapped-normal score tables on the torus (`.p.npy`, `.score.npy`) |
| `diffusion/score_model.py` | `TensorProductScoreModel` (e3nn), outputs one scalar score per rotatable bond |
| `diffusion/sampling.py` | RDKit ETKDG seeds → reverse SDE/ODE on torsions (20 steps default) |
| `utils/training.py`, `train.py` | denoising score-matching loop, Adam 1e-3, plateau scheduler |
| `generate_confs.py`, `evaluate_confs.py` | inference + AMR/COV evaluator (heavy-atom, symmetry-aware RMSD) |
| `utils/standardization.py` | GEOM conformer ↔ RDKit graph matching (“standardized pickles”) |
| `utils/boltzmann.py`, `utils/xtb.py` | torsional Boltzmann generator (energy-based training; xTB energies) |

## 3. Baseline architecture

```
SMILES (one fixed graph G) ─RDKit AddHs─► Mol(explicit H)
   │                                         │
   ├─featurize_mol: x ∈ R^{N×44} (QM9)       └─ETKDG seed conformer: bond lengths/angles/rings L (frozen)
   │   atom type one-hot, Z, aromatic, degree, hybridization,
   │   implicit valence, formal charge one-hot{-1,0,+1,other}, ring sizes, #rings
   ├─edge_attr: bond type one-hot {single,double,triple,aromatic}
   └─get_transformation_mask: rotatable bonds m (side to rotate)

reverse diffusion on τ ∈ T^m (σ: π → 0.01π, log schedule, 20 steps)
  score_model(x ⊕ emb(σ), radius graph(pos) ⊕ bond graph) ─e3nn TP convs (4 layers, ns=32, nv=8)─►
  bond-centred TP readout → one pseudoscalar score per rotatable bond (scaled by torus score norm)
  → Δτ applied rigidly by modify_conformer → conformer x = (L_RDKit, τ)
```
Training objective: denoising score matching on the torus, weighted by the wrapped-normal score norm (`utils/training.py`,
`diffusion/torus.py`). No discrete variables, no conditioning other than σ.

## 4. State assumptions in the input — what protonation would break

1. **One fixed graph per model input.** `ConformerDataset.featurize_mol` (`utils/dataset.py:190-229`) keeps a conformer only
   if its canonical (non-isomeric) heavy-atom SMILES equals the molecule's SMILES. A different protomer is a different graph
   and would be filtered out or treated as a different molecule. No notion of “same molecule, different state”.
2. **Hydrogens are explicit nodes.** `Chem.AddHs` (`featurization.py:117`) ⇒ protonation changes N (node count), the
   degree and implicit valence of the (de)protonated heavy atom, and the edge set. Formal charge is already a feature
   (`featurization.py:80`, one-hot over {-1,0,+1}+other) but **changing it alone is chemically wrong**: the H node and its
   bond must appear/disappear, and resonance/bond orders can change (carboxylate, amidinium/guanidinium, imidazolium,
   enolate/phenolate, tautomers move a double bond).
3. **Aromaticity/hybridization can change** (e.g. pyridine vs pyridinium keeps aromaticity; 2-hydroxypyridine ↔ 2-pyridone
   tautomer changes it; amine N sp3 ↔ amide-like sp2 is not a protonation effect but N-inversion sensitivity is).
4. **The torsion space is state-dependent.** `get_transformation_mask` (`utils/torsion.py:9-38`) marks a bond rotatable iff
   cutting it leaves a fragment with > 1 atom. Thus C–O(H) of –COOH, Ar–OH, –SH are torsions, but C–O⁻ of –COO⁻ and Ar–O⁻
   are not; –NH2 → –NH3⁺ keeps the torsion but changes its symmetry. **The dimension m and ordering of τ differ between
   states**, so `p(z, τ | G, h)` is a trans-dimensional distribution — τ cannot be a single shared vector across z.
   Torsion order follows `mol.GetBonds()` order, which changes when H atoms are added/removed.
5. **Local structure L comes from ETKDG of the given state.** Bond lengths/angles near the ionisable site differ by state
   (e.g. C–O in COO⁻ ≈ symmetric; N–C in amidinium); the prior work showed L error dominates TD's error on QM9
   (A1 oracle −54 %), so state-specific L matters at least as much as τ.
6. **Evaluation**: `evaluate_confs.py` matches generated vs GT conformers by molecule SMILES; RMSD is heavy-atom; a
   state-changing generator needs evaluation keyed on (molecule, state) and must also count wrong-state outputs.
7. **QM9/GEOM-QM9 is neutral small molecules** stored in one (neutral) state each; it supplies no protonation-state or pH
   supervision (see `06_dataset_audit.md`). TD-QM9 is therefore only the geometry reference, not a state benchmark.

## 5. Files

Could require modification (in a *new* package, not in place): featurization (state-aware graphs; pH/state features),
dataset (group conformers by parent molecule, attach state labels), torsion mask/indexing (per-state, with a map between
shared heavy-atom torsions across states), score model (conditioning inputs; optional state head), sampling (state
selection step / alternation), evaluator (state-aware matching + validity checks).

Must remain unchanged: `/scratch/nishanth.r/tordiff/**` — TD checkout (SHA-pinned by prior tests), `workdir/pretrained`,
`workdir/qm9_CTRL_*` checkpoints, `model_parameters.yml`, QM9 data and splits, `results/**`. The new work imports TD
read-only (as `flexitors/td.py` already does) or vendors a copy.

## 6. Minimal implementation plan (deferred until the literature/data audits justify a direction)

1. `jmcg/states.py`: state enumeration wrapper (Dimorphite-DL/other tool chosen in Phase VI) returning RDKit Mols with explicit
   H, sanitised, with a heavy-atom index map to the parent.
2. Per-state featurization = TD's `featurize_mol` unchanged (it already encodes formal charge, valence, degree), plus an
   optional global pH / state-context vector concatenated at `node_embedding` (`score_model.py`, the `torch.cat([data.x,
   node_sigma_emb])` site) — this is the cheapest conditioning point.
3. Torsion bookkeeping: per-state masks + heavy-atom-keyed torsion IDs so “shared torsions” are identifiable across states.
4. Sequential baseline first (state from external tool/oracle → TD-style generator per state), joint variants only if the
   audit supports them.

## 7. Open items for Phase I
- GPU peak memory / parameter count of the qm9_default model: **not yet measured here**; to be measured in the first smoke
  test (do not invent).
- Whether GEOM-DRUGS SMILES contain charged species (relevant: a TD-DRUGS model has seen formal charges) — family E agent.
