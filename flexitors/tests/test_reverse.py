"""Golden and pairing tests of the reverse loop and generation driver against TD's own code (gnode118, CPU, released
qm9_default checkpoint). plan_P2 G1, V-E3-2/3/4; EXP-3 §6 test 5; EXP-5 §9 (prefix property)."""
import copy
import os
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

MODEL = Path(os.environ.get('FLEXITORS_MODEL_DIR', '/scratch/nishanth.r/tordiff/workdir/pretrained/qm9_default'))
pytestmark = pytest.mark.cluster
SMILES = ['CC/C=C/C(C)O', 'CCC(=O)NC', 'OCC1CC1C=O']


@pytest.fixture(scope='module')
def td():
    if not MODEL.exists():
        pytest.skip(f'{MODEL} missing')
    import torch
    torch.set_num_threads(1)
    from flexitors.sampling.generate import load_model
    from flexitors.td import import_baseline
    model, margs = load_model(MODEL, 'best_model.pt', torch.device('cpu'))
    return Namespace(model=model, margs=margs, sampling=import_baseline('diffusion.sampling'),
                     tors=import_baseline('utils.torsion'), torch=torch)


def seeded(td, smi, n=6, seed=3):
    """TD's own seeds and perturb_seeds under one global seed (generate_confs.py:94-95, 125-151)"""
    import random
    random.seed(seed)
    np.random.seed(seed)
    td.torch.manual_seed(seed)
    mol, data = td.sampling.get_seed(smi, dataset='qm9')
    from rdkit.Chem import AllChem

    def embed_func(m, k):
        AllChem.EmbedMultipleConfs(m, numConfs=k, randomSeed=seed + 1)
        return m
    confs, _ = td.sampling.embed_seeds(mol, data, n, embed_func=embed_func)
    return mol, td.sampling.perturb_seeds(confs)


def test_reverse_without_hooks_equals_td_sample(td):
    from flexitors.sampling.reverse import reverse_diffusion
    for smi in SMILES:
        _, confs = seeded(td, smi)
        a, b = copy.deepcopy(confs), copy.deepcopy(confs)
        td.torch.manual_seed(11)
        td.sampling.sample(a, td.model, td.margs.sigma_max, td.margs.sigma_min, 20, 4)
        td.torch.manual_seed(11)
        nfe = reverse_diffusion(b, td.model, td.margs.sigma_max, td.margs.sigma_min, 20, 4, td.torch.device('cpu'))
        assert nfe == 20 * len(confs)
        for x, y in zip(a, b):
            assert td.torch.equal(x.pos, y.pos), smi


def test_o_hooks_change_nothing_and_k1_freezes(td):
    from flexitors.geometry.internal import dihedral
    from flexitors.sampling.reverse import reverse_diffusion
    from flexitors.sampling.stereo_guard import StereoHooks, constrainable_bonds
    mol, confs = seeded(td, 'CC/C=C/C(C)O')
    bonds, m = constrainable_bonds(mol, {b.GetIdx(): np.pi for b in mol.GetBonds()
                                         if b.GetBondTypeAsDouble() == 2.0})
    runs = {}
    for arm in (None, 'O', 'K1'):
        c = copy.deepcopy(confs)
        td.torch.manual_seed(5)
        hooks = StereoHooks(arm, bonds, m, 20) if arm else None
        reverse_diffusion(c, td.model, td.margs.sigma_max, td.margs.sigma_min, 20, 4, td.torch.device('cpu'),
                          hooks=hooks)
        runs[arm] = c
    for x, y in zip(runs[None], runs['O']):
        assert td.torch.equal(x.pos, y.pos)
    spec = [b for b in bonds if b.kind == 'specified'][0]
    for x, y in zip(confs, runs['K1']):
        q = np.array([spec.quad])
        assert dihedral(x.pos.numpy().astype(float), q)[0] == pytest.approx(dihedral(y.pos.numpy().astype(float), q)[0],
                                                                         abs=1e-5)


def test_update_sign_against_modify_conformer(td):
    from flexitors.geometry.internal import dihedral, wrap_angle
    from flexitors.sampling.stereo_guard import UPDATE_SIGN, constrainable_bonds
    mol, confs = seeded(td, 'CC/C=C/C(C)O')
    d = confs[0]
    bonds, m = constrainable_bonds(mol, {})
    for b in bonds:
        upd = np.zeros(m)
        upd[b.slot] = 0.2
        p = td.tors.modify_conformer(d.pos.clone(), d.edge_index.T[d.edge_mask], d.mask_rotate, upd)  # in place
        q = np.array([b.quad])
        delta = wrap_angle(dihedral(p.numpy().astype(float), q) - dihedral(d.pos.numpy().astype(float), q))[0]
        assert delta == pytest.approx(UPDATE_SIGN * 0.2, abs=1e-5)


def test_update_mask_equals_edge_mask_removal_for_other_scores(td):
    """EXP-3 §3.2 removes the bond from edge_mask; the update mask keeps the graph. The other bonds' scores agree."""
    from flexitors.sampling.stereo_guard import constrainable_bonds
    mol, confs = seeded(td, 'CC/C=C/C(C)O')
    bonds, m = constrainable_bonds(mol, {})
    spec = [b for b in bonds if b.kind != 'triple'][0]
    batch = td.sampling.DataLoader(td.sampling.InferenceDataset(copy.deepcopy(confs)), batch_size=1).__iter__().__next__()
    full, cut = copy.deepcopy(batch), copy.deepcopy(batch)
    masked = np.flatnonzero(batch.edge_mask.numpy())
    cut.edge_mask[int(masked[spec.slot])] = False
    for d in (full, cut):
        d.node_sigma = 0.5 * td.torch.ones(d.num_nodes)
    with td.torch.no_grad():
        a = td.model(full).edge_pred.numpy()
        b = td.model(cut).edge_pred.numpy()
    assert np.allclose(np.delete(a, spec.slot), b, atol=1e-6)


def test_conformer_trajectory_independent_of_pool_size_and_subset(td, tmp_path):
    """EXP-5 §9 / EXP-3 §9.1: a 2K run is the prefix of a larger run. Identical streams; the score model's float32
    reductions differ slightly with batch composition, so equality is to 1e-4 Å (CPU), and exact for a repeat run"""
    from flexitors.sampling.generate import run_conformers
    a = Namespace(seed=0, steps=20, batch_size=4, sigma_max=td.margs.sigma_max, sigma_min=td.margs.sigma_min)
    mol, data = td.sampling.get_seed('CC/C=C/C(C)O', dataset='qm9')
    small, _ = run_conformers(mol, data, 'k', a, td.model, td.torch.device('cpu'), {}, 'O', 0, 3)
    big, _ = run_conformers(mol, data, 'k', a, td.model, td.torch.device('cpu'), {}, 'O', 0, 8)
    for x, y in zip(small, big[:3]):
        assert np.abs(x.GetConformer().GetPositions() - y.GetConformer().GetPositions()).max() < 1e-4
    other, _ = run_conformers(mol, data, 'k', a, td.model, td.torch.device('cpu'), {}, 'O', 0, 3)
    assert np.array_equal(other[0].GetConformer().GetPositions(), small[0].GetConformer().GetPositions())
