"""Torsional Diffusion backend for explicit-state generation (frozen TD checkout imported read-only).

Reuses TD's own featurize_mol / get_transformation_mask / embed_seeds / perturb_seeds / sample / pyg_to_mol. The only
differences from TD's generate_confs.py are:
  - the input is the state Mol built by states.build_state (TD re-parses a SMILES; same graph, but we keep atom order
    and stereo exactly as specified);
  - ETKDG seed embedding is seeded (TD calls EmbedMultipleConfs without a seed -> not reproducible): conformer k is
    embedded on its own with RDKit randomSeed = conformer_seeds(seed, n)[k], drawn from numpy's generator. Passing
    randomSeed=seed to EmbedMultipleConfs is NOT used: with seed 0 RDKit returned n identical conformers (smoke 1,
    2026-10-28). Same default embedding parameters as TD otherwise;
  - optional `lock_stereo_double_bonds`: drop double bonds with a specified E/Z (STEREOE/Z/CIS/TRANS; never
    STEREONONE or STEREOANY) from the torsion set, so their E/Z cannot flip. True = every such bond; a collection of
    atom-index pairs = only those bonds. Only the locked bonds leave the torsion set: neighbouring single bonds stay
    rotatable and the model's per-bond readout for the remaining bonds is unchanged. Off by default here (TD
    behaviour, used by the wrapper-fidelity check); the run-level default is set in run.py (--ez_policy).
"""
import hashlib
import os
import random
import sys
from argparse import Namespace

import numpy as np
import torch
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

TD_REPO = os.environ.get('TD_REPO', '/scratch/nishanth.r/tordiff/torsional-diffusion')
if TD_REPO not in sys.path:
    sys.path.insert(0, TD_REPO)

# TD's diffusion/torus.py builds its score-norm table at import by Monte Carlo with the global, unseeded numpy RNG, so
# each process gets a slightly different table (relative differences ~1e-5) that scales every model output; the SDE
# amplified this to 0.03-1.7 A between processes for the same seed (drift diagnosis 2026-10-28). Importing TD under a
# fixed numpy seed makes the table identical in every process; the caller's numpy RNG state is restored afterwards.
TORUS_IMPORT_SEED = 20221024
if 'diffusion.torus' in sys.modules:
    raise ImportError('diffusion.torus was imported before jmcg.condgen.td_backend: its score-norm table is unseeded')
_np_state = np.random.get_state()
np.random.seed(TORUS_IMPORT_SEED)
from utils.featurization import featurize_mol, qm9_types, drugs_types  # noqa: E402  (TD, read-only)
from utils.torsion import get_transformation_mask  # noqa: E402
from utils.utils import get_model  # noqa: E402
from diffusion import sampling as td_sampling  # noqa: E402
from diffusion import torus as td_torus  # noqa: E402
np.random.set_state(_np_state)


def torus_table_sha():
    """Hash of TD's score-norm table actually in use (recorded in run metadata)."""
    return hashlib.sha256(np.ascontiguousarray(td_torus.score_norm_).tobytes()).hexdigest()[:16]

TYPES = {'qm9': qm9_types, 'drugs': drugs_types}
EZ_SPECIFIED = (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS,
                Chem.BondStereo.STEREOTRANS)
CHARGE_BUCKETS = (-1, 0, 1)  # TD featurize_mol: one-hot over these + "other"


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


class TDGenerator:
    def __init__(self, model_dir, ckpt='best_model.pt', steps=20, device='cpu'):
        with open(os.path.join(model_dir, 'model_parameters.yml')) as f:
            self.params = yaml.full_load(f)
        self.args = Namespace(**self.params)
        self.dataset = self.args.dataset
        self.types = TYPES[self.dataset]
        self.ckpt_path = os.path.join(model_dir, ckpt)
        self.ckpt_sha256 = sha256(self.ckpt_path)
        self.model = get_model(self.args)
        self.model.load_state_dict(torch.load(self.ckpt_path, map_location='cpu'), strict=True)
        self.model = self.model.to(device).eval()
        self.steps = steps
        self.device = device
        td_sampling.device = torch.device(device)

    def describe(self):
        return dict(backend='torsional_diffusion', td_repo=TD_REPO, checkpoint=self.ckpt_path,
                    checkpoint_sha256=self.ckpt_sha256, dataset=self.dataset,
                    in_node_features=self.args.in_node_features, atom_types=sorted(self.types, key=self.types.get),
                    sigma_min=float(self.args.sigma_min), sigma_max=float(self.args.sigma_max),
                    inference_steps=self.steps, sde='reverse SDE (TD default, ode=False)', device=self.device,
                    torus_import_seed=TORUS_IMPORT_SEED, torus_score_norm_sha256=torus_table_sha())

    def compatibility(self, mol):
        """Precise reasons this checkpoint cannot represent the state graph (hard), plus out-of-distribution notes."""
        hard, ood = [], []
        bad = sorted({a.GetSymbol() for a in mol.GetAtoms() if a.GetSymbol() not in self.types})
        if bad:
            hard.append(f'atom types {bad} not in the {self.dataset} checkpoint vocabulary {sorted(self.types)}')
        q = sorted({a.GetFormalCharge() for a in mol.GetAtoms()} - set(CHARGE_BUCKETS))
        if q:
            ood.append(f'formal charges {q} fall in the featuriser "other" bucket')
        if not hard:
            n = featurize_mol(mol, self.types).x.shape[1]
            if n != self.args.in_node_features:
                hard.append(f'featuriser gives {n} node features, checkpoint expects {self.args.in_node_features}')
        if self.dataset == 'qm9' and mol.GetNumHeavyAtoms() > 9:
            ood.append(f'{mol.GetNumHeavyAtoms()} heavy atoms > 9 (QM9 training size)')
        if any(a.GetFormalCharge() for a in mol.GetAtoms()):
            ood.append('charged atoms: representable by the featuriser, but whether the checkpoint saw such graphs '
                       'is a training-data question (inventory.json), not shown by successful inference')
        return hard, ood

    def size_or_bucket_ood(self, mol):
        q = {a.GetFormalCharge() for a in mol.GetAtoms()} - set(CHARGE_BUCKETS)
        return bool(q) or (self.dataset == 'qm9' and mol.GetNumHeavyAtoms() > 9)

    def torsion_data(self, mol, lock_stereo_double_bonds=False):
        data = featurize_mol(mol, self.types)
        edge_mask, mask_rotate = get_transformation_mask(data)
        locked = []
        if lock_stereo_double_bonds:
            want = None if lock_stereo_double_bonds is True else {frozenset(map(int, p)) for p in lock_stereo_double_bonds}
            ei = data.edge_index.T.numpy()
            rows = np.flatnonzero(edge_mask)
            keep = np.ones(len(rows), dtype=bool)
            for k, e in enumerate(rows):
                i, j = int(ei[e, 0]), int(ei[e, 1])
                b = mol.GetBondBetweenAtoms(i, j)
                if (b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() in EZ_SPECIFIED
                        and (want is None or frozenset((i, j)) in want)):
                    keep[k] = False
                    locked.append((i, j))
            edge_mask = edge_mask.copy()
            edge_mask[rows[~keep]] = False
            mask_rotate = mask_rotate[keep]
        data.edge_mask = torch.tensor(edge_mask)
        data.mask_rotate = mask_rotate
        return data, locked

    def generate(self, mol, n_confs, seed, lock_stereo_double_bonds=False):
        """Returns (list of RDKit Mols, one conformer each; info). Raises RuntimeError on embedding failure."""
        random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
        mol = Chem.Mol(mol)
        data, locked = self.torsion_data(mol, lock_stereo_double_bonds)
        n_rot = int(data.edge_mask.sum())
        n_gen = max(n_confs, 2)  # TD's edge_pred.squeeze() breaks for a 1-conformer batch with 1 torsion

        def embed_func(m, num):
            return embed_n(m, num, seed)

        conformers, _ = td_sampling.embed_seeds(mol, data, n_gen, embed_func=embed_func)
        if not conformers:
            raise RuntimeError('ETKDG seed embedding failed')
        if n_rot > 0:
            conformers = td_sampling.perturb_seeds(conformers)
            conformers = td_sampling.sample(conformers, self.model, self.args.sigma_max, self.args.sigma_min,
                                            self.steps, batch_size=n_gen, ode=False)
        mols = [td_sampling.pyg_to_mol(mol, c, mmff=False, rmsd=False) for c in conformers][:n_confs]
        return mols, dict(n_rotatable_td=n_rot, model_applied=n_rot > 0, locked_double_bonds=locked,
                          edge_index=data.edge_index, edge_mask=data.edge_mask.numpy())


def conformer_seeds(seed, n):
    """One independent RDKit randomSeed per conformer (all > 0), reproducible from the user seed."""
    return [int(x) for x in np.random.default_rng(seed).integers(1, 2 ** 31 - 1, n)]


def embed_n(mol, n, seed, params=None):
    """Add n conformers to mol in place, conformer k embedded alone with conformer_seeds(seed, n)[k].
    params=None: RDKit EmbedMolecule defaults (= TD's EmbedMultipleConfs defaults); else an EmbedParameters factory."""
    mol.RemoveAllConformers()
    for s in conformer_seeds(seed, n):
        t = Chem.Mol(mol)
        if params is None:
            cid = AllChem.EmbedMolecule(t, randomSeed=s)
        else:
            p = params()
            p.randomSeed = s
            cid = AllChem.EmbedMolecule(t, p)
        if cid < 0:
            break
        mol.AddConformer(Chem.Conformer(t.GetConformer(cid)), assignId=True)
    return mol


def etkdg_baseline(mol, n_confs, seed):
    """Clearly labelled baseline: RDKit ETKDGv3, same Mol, same per-conformer seeds as the TD seed embedding, no
    force-field refinement (TD seeds have none)."""
    m = embed_n(Chem.Mol(mol), n_confs, seed, AllChem.ETKDGv3)
    if m.GetNumConformers() != n_confs:
        raise RuntimeError(f'ETKDGv3 embedded {m.GetNumConformers()}/{n_confs}')
    out = []
    for conf in m.GetConformers():
        c = Chem.Mol(m)
        c.RemoveAllConformers()
        c.AddConformer(Chem.Conformer(conf), assignId=True)
        out.append(c)
    return out
