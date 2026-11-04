"""EXP-1 V2 building blocks (§5.3 test_projector torch part, test_parity, test_golden_td, test_transform, loss limits).
All need torch + the TD environment (torch_geometric, e3nn, torus tables): marked ``cluster``."""
import math

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

pytestmark = pytest.mark.cluster

SMILES = ('CC(O)CC1CCNC1C(=O)N', 'CCOC(=O)C=CC')


def featurised(smiles, seed=3):
    import torch
    from flexitors.td import import_baseline
    feat = import_baseline('utils.featurization')
    tors = import_baseline('utils.torsion')
    mol, data = feat.featurize_mol_from_smiles(smiles, dataset='qm9')
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    data.pos = torch.tensor(mol.GetConformer().GetPositions(), dtype=torch.float32)
    em, mr = tors.get_transformation_mask(data)
    data.edge_mask, data.mask_rotate = torch.tensor(em), mr
    return mol, data


def model_args():
    from flexitors.training.loop import parse_args
    return parse_args(['--cache_train', 'x', '--cache_val', 'x', '--log_dir', 'x'])


def batch_of(datas, t=0.3, sigma=0.5):
    import torch
    from torch_geometric.data import Batch, Data
    ds = []
    for d in datas:
        e = Data(x=d.x, edge_index=d.edge_index, edge_attr=d.edge_attr, edge_mask=d.edge_mask, pos=d.pos.clone())
        e.node_sigma = sigma * torch.ones(len(d.pos)); e.node_lambda = t * torch.ones(len(d.pos))
        ds.append(e)
    return Batch.from_data_list(ds)


def test_batched_projector_matches_numpy():
    import torch
    from flexitors.geometry.l_step import w_quadruples
    from flexitors.geometry.projector import dihedral_preserving_projector
    from flexitors.geometry.projector_torch import BatchProjector
    mols = [featurised(s) for s in SMILES]
    pos, batch, quads, per = [], [], [], []
    off = 0
    rng = np.random.default_rng(0)
    for g, (mol, d) in enumerate(mols):
        x = d.pos.double().numpy()
        z = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
        ei = d.edge_index.numpy()
        nb = [[] for _ in z]
        for u, v in ei[:, ::2].T:
            nb[u].append(int(v)); nb[v].append(int(u))
        q, _, _, _ = w_quadruples(z, ei.T[d.edge_mask.numpy()], nb, x)
        v = rng.normal(size=x.shape)
        P, r = dihedral_preserving_projector(x, q)
        per.append(((P @ v.reshape(-1)).reshape(x.shape), r, v))
        pos.append(x); batch += [g] * len(x); quads.append(q + off); off += len(x)
    proj = BatchProjector(torch.tensor(np.vstack(pos)), torch.tensor(batch), torch.tensor(np.vstack(quads)))
    out = proj.project(torch.tensor(np.vstack([p[2] for p in per]))).numpy()
    assert np.allclose(out, np.vstack([p[0] for p in per]), atol=1e-8)
    assert proj.rank.tolist() == [p[1] for p in per]


def test_golden_torsion_output_equals_td_and_head_is_1o():
    import torch
    from flexitors.models.flexitors_score import _td_score_model, build_model
    args = model_args()
    torch.manual_seed(0)
    m = build_model(args).eval()
    td = _td_score_model()(in_node_features=44, in_edge_features=4, ns=32, nv=8, sigma_embed_dim=32,
                           sigma_min=args.sigma_min, sigma_max=args.sigma_max, num_conv_layers=4, max_radius=5.0,
                           radius_embed_dim=50, scale_by_sigma=True, use_second_order_repr=True,
                           lambda_embed_dim=32).eval()
    missing = td.load_state_dict({k: v for k, v in m.state_dict().items() if not k.startswith('l_head.')}, strict=True)
    assert not missing.missing_keys
    datas = [featurised(s)[1] for s in SMILES]
    with torch.no_grad():
        a = m(batch_of(datas)); b = td(batch_of(datas))
        assert torch.equal(a.edge_pred, b.edge_pred)  # geometry head present, torsion output unchanged (A3 golden)
        l0 = a.l_pred.clone()
        for d in datas:
            d.pos = -d.pos
        assert torch.allclose(m(batch_of(datas)).l_pred, -l0, atol=1e-5)  # parity: 1o negates under x -> -x
        R = torch.tensor([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]) @ torch.tensor(
            [[1.0, 0.0, 0.0], [0.0, math.cos(0.7), -math.sin(0.7)], [0.0, math.sin(0.7), math.cos(0.7)]])
        for d in datas:
            d.pos = -d.pos @ R.T
        assert torch.allclose(m(batch_of(datas)).l_pred, l0 @ R.T, atol=1e-4)  # rotation equivariance


def paired_datapoint(smiles='CC1CCC(O)CC1C=CC'):
    import torch
    mol, d = featurised(smiles)
    m2 = Chem.Mol(mol)
    AllChem.MMFFOptimizeMolecule(m2, maxIters=300)
    d.pos = [d.pos]
    d.gt_pos = [torch.tensor(m2.GetConformer().GetPositions(), dtype=torch.float32)]
    d.pair_ok, d.weights, d.mol = [True], [1.0], mol
    return d


def test_transform_same_dtau_on_both_endpoints():
    import torch
    from flexitors.geometry.internal import dihedral, wrap_angle
    from flexitors.td import import_baseline
    from flexitors.training.transforms import PairedJointTransform
    np.random.seed(1)
    d = paired_datapoint()
    xR = d.pos[0].numpy().astype(np.float64)
    out = PairedJointTransform(sigma_b=0.15)(d)
    q = out.w_quad_index.T.numpy()
    assert out.w_quad_index.shape[0] == 4 and len(q) > 0 and float(out.geo_mask) == 1.0
    tors = import_baseline('utils.torsion')
    em = np.asarray(out.edge_mask, dtype=bool)
    xRp = tors.modify_conformer(xR.copy(), torch.as_tensor(out.edge_index.numpy().T[em]), out.mask_rotate,
                                out.edge_rotate.numpy(), as_numpy=True)
    # x_G*' has x_R''s W dihedrals: the same Δτ moved both endpoints
    assert np.abs(wrap_angle(dihedral(out.geo_target.double().numpy(), q) - dihedral(xRp, q))).max() < 1e-3


def test_loss_limits():
    import torch
    from flexitors.training.losses import projected_l_loss
    from flexitors.training.transforms import PairedJointTransform
    from torch_geometric.data import Batch
    np.random.seed(2)
    outs = [PairedJointTransform(sigma_b=0.15)(paired_datapoint(s)) for s in ('CC1CCC(O)CC1C=CC', 'OC1CCCCC1CN')]
    b = Batch.from_data_list(outs)
    b.l_pred = (b.geo_target - b.pos).float()
    loss, base, n = projected_l_loss(b)
    assert n == 2 and float(loss) < 1e-10 < float(base)
    b.l_pred = torch.zeros_like(b.pos)
    loss0, base0, _ = projected_l_loss(b)
    assert abs(float(loss0) - float(base0)) < 1e-12


def test_geo_norm_base_objective():
    import torch
    from flexitors.training.loop import geo_objective
    lg = torch.tensor(0.02, dtype=torch.float64, requires_grad=True)
    bg = torch.tensor(0.01, dtype=torch.float64, requires_grad=True)
    go = geo_objective(lg, bg, 'base')
    assert abs(float(go) - 2.0) < 1e-12
    go.backward()
    assert abs(float(lg.grad) - 100.0) < 1e-9 and bg.grad is None  # base detached
    assert geo_objective(lg, bg, 'none') is lg
    z = torch.tensor(0.0, dtype=torch.float64)
    assert float(geo_objective(z, z, 'base')) == 0.0  # no active molecule: raw (0), no division


def test_zero_init_head_starts_at_zero_and_learns():
    import torch
    from flexitors.models.flexitors_score import build_model
    args = model_args()
    args.l_head_init = 'zero'
    torch.manual_seed(0)
    m = build_model(args)
    out = m(batch_of([featurised(s)[1] for s in SMILES]))
    assert torch.count_nonzero(out.l_pred) == 0
    (out.l_pred * torch.randn_like(out.l_pred)).sum().backward()
    assert float(m.l_head.weight.grad.abs().sum()) > 0  # weights get a gradient although Δ̂ = 0


def test_loss_without_projector_is_plain_mse():
    import torch
    from flexitors.training.losses import projected_l_loss
    from flexitors.training.transforms import PairedJointTransform
    from torch_geometric.data import Batch
    np.random.seed(4)
    b = Batch.from_data_list([PairedJointTransform(sigma_b=0.15)(paired_datapoint(s))
                              for s in ('CC1CCC(O)CC1C=CC', 'OC1CCCCC1CN')])
    b.l_pred = torch.randn_like(b.pos)
    loss, base, n = projected_l_loss(b, project=False)
    r = (b.pos - b.geo_target).double() + b.l_pred.double()
    per = [float((r[b.batch == g] ** 2).sum() / (3 * int((b.batch == g).sum()))) for g in range(2)]
    assert n == 2 and abs(float(loss) - np.mean(per)) < 1e-9
    b.l_pred = (b.geo_target - b.pos).float()
    assert float(projected_l_loss(b, project=False)[0]) < 1e-10


def test_rigid_molecule_geometry_head():
    """No rotatable bond: TD's bond head would fail; l_pred must equal the one computed with any bond marked rotatable
    (the atom layers do not read edge_mask), edge_pred is empty and edge_mask is restored."""
    import torch
    from flexitors.models.flexitors_score import build_model
    args = model_args()
    torch.manual_seed(0)
    m = build_model(args).eval()
    _, d = featurised('O=C1N[C@H]2C[C@H]2O1')
    assert not bool(d.edge_mask.any())
    with torch.no_grad():
        a = m(batch_of([d]))
        e = d.clone(); e.edge_mask = torch.zeros_like(d.edge_mask); e.edge_mask[0] = True
        b = m(batch_of([e]))
    assert torch.equal(a.l_pred, b.l_pred) and a.edge_pred.numel() == 0
    assert not bool(a.edge_mask.any())


def test_td_loss_all_rigid_batch_is_zero_not_nan():
    import torch
    from torch_geometric.data import Data
    from flexitors.training.losses import td_loss
    d = Data(edge_pred=torch.zeros(0, requires_grad=True), edge_rotate=torch.zeros(0), edge_sigma=torch.ones(1))
    lt, bt = td_loss(d)
    assert float(lt) == 0.0 and float(bt) == 0.0
    lt.backward()


def test_rigid_molecule_in_mixed_batch():
    """A rigid molecule batched with a flexible one: TD's own path runs (some bond is rotatable); the flexible molecule's
    torsion output and both molecules' geometry outputs equal those computed alone."""
    import torch
    from flexitors.models.flexitors_score import build_model
    args = model_args()
    torch.manual_seed(0)
    m = build_model(args).eval()
    _, rigid = featurised('O=C1N[C@H]2C[C@H]2O1')
    _, flex = featurised(SMILES[1])
    with torch.no_grad():
        both = m(batch_of([flex, rigid]))
        a, b = m(batch_of([flex])), m(batch_of([rigid]))
    n = len(flex.pos)
    assert torch.allclose(both.edge_pred, a.edge_pred, atol=1e-5)
    assert torch.allclose(both.l_pred[:n], a.l_pred, atol=1e-5) and torch.allclose(both.l_pred[n:], b.l_pred, atol=1e-5)
