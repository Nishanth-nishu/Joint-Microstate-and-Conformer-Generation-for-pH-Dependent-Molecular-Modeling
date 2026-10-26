"""Generation CLI for every Stage A inference arm (Exp 3 and Exp 5); writes the ``{corrected_smiles: [Mol]}`` pickle of
``generate_confs.py:210-213`` so the unchanged ``evaluate_confs.py`` scores every arm (EVAL §0 rule 1).

Per molecule it follows ``generate_confs.py:125-180`` (``get_seed`` from the csv's corrected SMILES, ETKDG seeds,
``perturb_seeds`` unless m = 0, the reverse SDE unless m = 0, ``pyg_to_mol``), with per-conformer random streams
(``seeding.py``) and these arms:
* ``O``, ``K1``, ``K1_wrong``, ``K1_canonical``, ``K2_inf``: EXP-3 §5 / §9.2 operators (``stereo_guard.py``);
* ``F_filter``: O, then conformers whose geometric E/Z label differs from a SMILES target are replaced by fresh
  conformers (further conformer indices) for at most 3 extra rounds, extra NFE counted; if fewer than 2K valid remain,
  the output is filled to 2K with the least-violating invalid ones, flagged and counted (EXP-3 §3.4, §9.2; E-CHEM-064);
* ``PG``: Particle Guidance through TD's own ``sample()`` with one of the four README settings
  (``torsional-diffusion/README.md:80-91``, E-SAMPLE-058) and ``--batch_size 128`` (EXP-5 §12.5, E-SAMPLE-033); its
  global RNG is seeded per molecule.
Sample counts (``--n``): ``2K`` (= 2 × csv n_conformers, TD's rule), ``4K``, or ``pool`` = max(128, 4K) (EXP-5 §12.1).
A conformer's trajectory depends only on (seed, molecule, conformer index), so a 2K run is the prefix of a pool and a
subset run (``--keys``) reproduces the full run on those molecules (EXP-3 §9.1).
Outputs (``--out_dir``): ``confs.pkl``, ``seeds.pkl`` (ETKDG seed positions per conformer: the same-label validity
reference, EVAL §3.1), ``meta.json`` (per-molecule counters, replacements, fills), ``budget.json`` (NFE, wall time;
EVAL §8, EXP-5 §4.1), ``population.txt`` and ``provenance.json``. Generation is resumable (``confs.partial.pkl``).
"""
from __future__ import annotations

import copy
import json
import pickle
import time
from argparse import ArgumentParser, Namespace
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from flexitors.data.ez_targets import load_targets
from flexitors.geometry.internal import dihedral
from flexitors.provenance import provenance, write_provenance
from flexitors.sampling.reverse import reverse_diffusion
from flexitors.sampling.seeding import conformer_seed, etkdg_seed, numpy_rng
from flexitors.sampling.stereo_guard import ARMS, StereoHooks, constrainable_bonds
from flexitors.td import import_baseline

RDLogger.DisableLog('rdApp.*')
F_FILTER_ROUNDS = 3   # EXP-3 §3.4 "resample until K = 2L valid or 3 extra rounds"
PG_BATCH = 128        # EXP-5 §12.5 / E-SAMPLE-033
PG_SETTINGS = {       # torsional-diffusion/README.md:80-91 (tuned for drugs_default; chosen on QM9 validation, §12.5)
    'inv_recall': dict(pg_invariant=True, pg_kernel_size_log_0=1.7565691770646286,
                       pg_kernel_size_log_1=1.1960868735428605, pg_langevin_weight_log_0=-2.2245183818892103,
                       pg_langevin_weight_log_1=-2.403905082248579, pg_repulsive_weight_log_0=-2.158537381110402,
                       pg_repulsive_weight_log_1=-2.717482077162461, pg_weight_log_0=0.8004013644746992,
                       pg_weight_log_1=-0.9255658381081596),
    'inv_precision': dict(pg_invariant=True, pg_kernel_size_log_0=-0.9686202580381296,
                          pg_kernel_size_log_1=-0.7808409291022302, pg_langevin_weight_log_0=-2.434216242826782,
                          pg_langevin_weight_log_1=-0.2602238633333869, pg_repulsive_weight_log_0=-2.0439285313973237,
                          pg_repulsive_weight_log_1=-1.468234554877924, pg_weight_log_0=0.3495680598729498,
                          pg_weight_log_1=-0.22001939454654185),
    'noninv_recall': dict(pg_invariant=False, pg_kernel_size_log_0=2.35958, pg_kernel_size_log_1=-0.78826,
                          pg_langevin_weight_log_0=-1.55054, pg_langevin_weight_log_1=-2.70316,
                          pg_repulsive_weight_log_0=1.01317, pg_repulsive_weight_log_1=-2.68407,
                          pg_weight_log_0=0.60504, pg_weight_log_1=-1.15020),
    'noninv_precision': dict(pg_invariant=False, pg_kernel_size_log_0=1.29503, pg_kernel_size_log_1=1.45944,
                             pg_langevin_weight_log_0=-2.88867, pg_langevin_weight_log_1=-2.47591,
                             pg_repulsive_weight_log_0=-1.01222, pg_repulsive_weight_log_1=-1.91253,
                             pg_weight_log_0=-0.16253, pg_weight_log_1=0.79355),
}


def n_samples(mode: str, n_csv: int) -> int:
    """Conformers per molecule: TD's 2K, 4K, or the K-free-block pool max(128, 4K) (EXP-5 §12.1)."""
    return {'2K': 2 * n_csv, '4K': 4 * n_csv, 'pool': max(128, 4 * n_csv)}[mode]


def load_model(model_dir: Path, ckpt: str, device):
    """TD's released model exactly as ``generate_confs.py:96-112`` builds it (yaml arguments, strict load)."""
    import torch
    import yaml
    with open(Path(model_dir) / 'model_parameters.yml') as f:
        args = Namespace(**yaml.full_load(f))
    model = import_baseline('utils.utils').get_model(args)
    model.load_state_dict(torch.load(Path(model_dir) / ckpt, map_location='cpu'), strict=True)
    return model.to(device).eval(), args


def embed(mol: Chem.Mol, data, key: str, seed: int, indices: range) -> list | None:
    """ETKDG seeds, one embedding per conformer index with its own seed (``sampling.py:63-102`` without seed files);
    None if any embedding fails, as TD treats a short embedding (``sampling.py:71-73``)."""
    from rdkit.Chem import AllChem
    import torch
    out = []
    for j in indices:
        m = Chem.Mol(mol)
        m.RemoveAllConformers()
        if AllChem.EmbedMolecule(m, randomSeed=etkdg_seed(seed, key, j)) != 0:
            return None
        d = copy.deepcopy(data)
        d.pos = torch.from_numpy(m.GetConformer().GetPositions()).float()
        d.seed_mol = copy.deepcopy(m)
        out.append(d)
    return out


def seed_stage(confs: list, key: str, seed: int, first: int, hooks: StereoHooks, tors) -> None:
    """``perturb_seeds`` (sampling.py:105-113) with per-conformer uniform draws, the arm's freeze mask and seed check,
    then the K2 fold of the randomised seeds."""
    for j, d in enumerate(confs):
        m = int(d.edge_mask.sum())
        uni = numpy_rng(seed, key, first + j, 'perturb').uniform(low=-np.pi, high=np.pi, size=m)
        pos = d.pos.numpy().astype(np.float64)[None]
        upd = hooks.seed_updates(pos, uni[None])[0]
        d.pos = tors.modify_conformer(d.pos, d.edge_index.T[d.edge_mask], d.mask_rotate, upd)
        fold = hooks.fold_updates([j], d.pos.numpy().astype(np.float64)[None], None, None)[0]
        if np.any(fold != 0):
            d.pos = tors.modify_conformer(d.pos, d.edge_index.T[d.edge_mask], d.mask_rotate, fold)
            upd = upd + fold
        d.total_perturb = upd


class ConformerNoise:
    """Per-conformer step noise z ~ N(0, 1) from each conformer's own stream (seeding.py)."""

    def __init__(self, seed: int, key: str, first: int, m: int):
        self.seed, self.key, self.first, self.m = seed, key, first, m
        self.rngs: dict[int, np.random.Generator] = {}

    def __call__(self, idx: list[int], step: int) -> np.ndarray:
        rows = []
        for i in idx:
            if i not in self.rngs:
                self.rngs[i] = numpy_rng(self.seed, self.key, self.first + i, 'noise')
            rows.append(self.rngs[i].standard_normal(self.m))
        return np.stack(rows).astype(np.float32)


def run_conformers(mol, data, key: str, a: Namespace, model, device, targets: dict[int, float], arm: str,
                   first: int, n: int) -> tuple[list | None, dict]:
    """Generate conformer indices [first, first + n) of one molecule under ``arm``; returns (mols, info)."""
    sampling = import_baseline('diffusion.sampling')
    tors = import_baseline('utils.torsion')
    confs = embed(mol, data, key, a.seed, range(first, first + n))
    if confs is None:
        return None, {'failure': 'embed'}
    m = int(data.edge_mask.sum())
    bonds, n_slots = constrainable_bonds(mol, targets)
    assert n_slots == m, f'{key}: TD rotatable count {m} != flexitors {n_slots}'
    hooks = StereoHooks(arm if arm in ARMS else 'O', bonds, m, a.steps)
    info = {'m': m, 'nfe': 0, 'pg_wall_s': 0.0}
    if m > 0:  # generate_confs.py:150-151, 153: molecules without rotatable bonds keep their seeds
        seed_stage(confs, key, a.seed, first, hooks, tors)
        if arm == 'PG':
            import random
            import torch
            h = conformer_seed(a.seed, key, first, 'pg')
            random.seed(h)
            np.random.seed(h % 2 ** 32)
            torch.manual_seed(h)
            t0 = time.time()
            # the invariant kernel reads a conformer of `mol` (spyrmsd from_rdkit, sampling.py:148); TD's `mol` keeps
            # its first embedding (sampling.py:101-102), here the first seed of this run
            sampling.sample(confs, model, a.sigma_max, a.sigma_min, a.steps, PG_BATCH, mol=Chem.Mol(confs[0].seed_mol),
                            **PG_SETTINGS[a.pg_setting])
            info['pg_wall_s'] = time.time() - t0
            info['nfe'] = n * a.steps
        else:
            info['nfe'] = reverse_diffusion(confs, model, a.sigma_max, a.sigma_min, a.steps, a.batch_size, device,
                                            noise=ConformerNoise(a.seed, key, first, m), hooks=hooks)
    info['counts'] = {k: v.tolist() for k, v in hooks.counts.items()}
    mols = [sampling.pyg_to_mol(mol, d, mmff=False, rmsd=False) for d in confs]
    info['seeds'] = np.stack([d.seed_mol.GetConformer().GetPositions() for d in confs]).astype(np.float32)
    return mols, info


def violation(mol: Chem.Mol, quads: np.ndarray, phi: np.ndarray) -> float:
    """Largest |δ_b| − π/2 over the specified bonds (> 0: E/Z label differs from the target)."""
    tau = dihedral(mol.GetConformer().GetPositions(), quads)
    return float(np.max(np.abs(np.arctan2(np.sin(tau - phi), np.cos(tau - phi))) - np.pi / 2))


def f_filter(mol, data, key, a, model, device, targets, n) -> tuple[list | None, dict]:
    """EXP-3 §3.4 / §9.2 F-filter on top of O for one molecule."""
    mols, info = run_conformers(mol, data, key, a, model, device, targets, 'O', 0, n)
    if mols is None or not targets:
        return mols, {**info, 'replaced': 0, 'filled': 0, 'rounds': 0}
    bonds, _ = constrainable_bonds(mol, targets)
    spec = [b for b in bonds if b.kind == 'specified']
    quads, phi = np.array([b.quad for b in spec]), np.array([b.phi for b in spec])
    pool = [(violation(x, quads, phi), x, s) for x, s in zip(mols, info['seeds'])]
    valid = [p for p in pool if p[0] <= 0]
    invalid = [p for p in pool if p[0] > 0]
    nxt, rounds, extra_nfe = n, 0, 0
    while len(valid) < n and rounds < F_FILTER_ROUNDS:
        need = n - len(valid)
        new, inf = run_conformers(mol, data, key, a, model, device, targets, 'O', nxt, need)
        nxt, rounds, extra_nfe = nxt + need, rounds + 1, extra_nfe + inf.get('nfe', 0)
        if new is None:
            break
        for x, s in zip(new, inf['seeds']):
            v = violation(x, quads, phi)
            (valid if v <= 0 else invalid).append((v, x, s))
    filled = max(0, n - len(valid))
    out = valid[:n] + sorted(invalid, key=lambda p: p[0])[:filled]
    info.update(nfe=info['nfe'] + extra_nfe, extra_nfe=extra_nfe, replaced=nxt - n, filled=filled, rounds=rounds,
                fill_columns=list(range(n - filled, n)), seeds=np.stack([p[2] for p in out]))
    return [p[1] for p in out], info


def main() -> None:
    ap = ArgumentParser(description='Stage A generation (Exp 3 arms, Exp 5 pools)')
    ap.add_argument('--model_dir', required=True)
    ap.add_argument('--ckpt', default='best_model.pt')
    ap.add_argument('--csv', required=True, help='test-format csv: smiles, n_conformers, corrected_smiles')
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--arm', required=True, choices=list(ARMS) + ['F_filter', 'PG'])
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--n', default='2K', choices=['2K', '4K', 'pool'])
    ap.add_argument('--keys', default=None, help='file with corrected smiles to generate (affected subset)')
    ap.add_argument('--ez_table', default=None, help='round5/data/ez_targets.csv (Exp 3 arms)')
    ap.add_argument('--pg_setting', default=None, choices=sorted(PG_SETTINGS))
    ap.add_argument('--steps', type=int, default=20)
    ap.add_argument('--batch_size', type=int, default=32)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--checkpoint_every', type=int, default=25)
    a = ap.parse_args()
    if (a.arm == 'PG') != (a.pg_setting is not None):
        raise SystemExit('--pg_setting is required for, and only for, --arm PG')
    if a.arm in ('K1', 'K1_wrong', 'K2_inf', 'F_filter', 'K1_canonical') and a.ez_table is None:
        raise SystemExit(f'--arm {a.arm} needs --ez_table')
    import torch
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    model, margs = load_model(Path(a.model_dir), a.ckpt, device)
    a.sigma_max, a.sigma_min = margs.sigma_max, margs.sigma_min
    sampling = import_baseline('diffusion.sampling')
    targets = load_targets(Path(a.ez_table) if a.ez_table else None)
    df = pd.read_csv(a.csv)
    if a.keys:
        keys = {ln.strip() for ln in open(a.keys) if ln.strip()}
        df = df[df.corrected_smiles.isin(keys)]
    if a.limit:
        df = df.iloc[:a.limit]
    run_id = json.dumps({k: v for k, v in vars(a).items() if k not in ('limit', 'checkpoint_every')}, sort_keys=True)
    part = out / 'confs.partial.pkl'
    state = {'run_id': run_id, 'confs': {}, 'seeds': {}, 'meta': {}}
    if part.exists():
        with open(part, 'rb') as f:
            old = pickle.load(f)
        if old['run_id'] == run_id:
            state = old
    t0 = time.time()
    for row_i, (raw, n_csv, smi) in enumerate(df[['smiles', 'n_conformers', 'corrected_smiles']].values):
        if smi in state['meta']:
            continue
        n = n_samples(a.n, int(n_csv))
        mol, data = sampling.get_seed(smi, dataset=margs.dataset)
        if not mol:
            state['meta'][smi] = {'failure': 'seed'}
            continue
        tg = targets.get(smi, {})
        if a.arm == 'F_filter':
            mols, info = f_filter(mol, data, smi, a, model, device, tg, n)
        else:
            mols, info = run_conformers(mol, data, smi, a, model, device, tg, a.arm, 0, n)
        if mols is not None:
            state['confs'][smi] = mols
            state['seeds'][smi] = info.pop('seeds')
        state['meta'][smi] = info
        if (row_i + 1) % a.checkpoint_every == 0:
            with open(part, 'wb') as f:
                pickle.dump(state, f)
    meta = state['meta']
    with open(out / 'confs.pkl', 'wb') as f:
        pickle.dump(state['confs'], f)
    with open(out / 'seeds.pkl', 'wb') as f:
        pickle.dump(state['seeds'], f)
    (out / 'meta.json').write_text(json.dumps(meta, default=float) + '\n')
    fails = {k: v['failure'] for k, v in meta.items() if 'failure' in v}
    budget = dict(nfe=int(sum(v.get('nfe', 0) for v in meta.values())),
                  extra_nfe=int(sum(v.get('extra_nfe', 0) for v in meta.values())),
                  n_conformers=int(sum(len(v) for v in state['confs'].values())),
                  pg_wall_s=float(sum(v.get('pg_wall_s', 0.0) for v in meta.values())),
                  wall_s=time.time() - t0, device=str(device), arm=a.arm, n_mode=a.n, seed=a.seed)
    (out / 'budget.json').write_text(json.dumps(budget, indent=1) + '\n')
    counts = {k: int(sum(sum(v['counts'][k]) for v in meta.values() if 'counts' in v))
              for k in ('ez_violations', 'planarity_violations', 'fold_corrections')}
    (out / 'population.txt').write_text(
        f'rows {len(df)}; generated {len(state["confs"])}; failures {len(fails)} {sorted(set(fails.values()))}\n'
        f'counters {counts}\nreplaced {sum(v.get("replaced", 0) for v in meta.values())} '
        f'filled {sum(v.get("filled", 0) for v in meta.values())} '
        f'molecules_with_fill {sum(1 for v in meta.values() if v.get("filled", 0))}\n'
        + ''.join(f'failure {k} {v}\n' for k, v in sorted(fails.items())))
    write_provenance(out)
    if part.exists():
        part.unlink()
    print('SUMMARY', json.dumps({**budget, **counts, 'failures': len(fails), 'commit': provenance()['commit']}))


if __name__ == '__main__':
    main()
