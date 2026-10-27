"""Tests for the TD backend, SDF round trip and geometry validation (needs the TD venv + a checkpoint).

Run: python -m jmcg.condgen.test_backend --model_dir <TD qm9 model dir> --states_csv jmcg/condgen/smoke_states.csv
"""
import argparse
import os
import tempfile

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from jmcg.condgen import states as st
from jmcg.condgen import validate as va
from jmcg.condgen.run import write_sdf
from jmcg.condgen.td_backend import TDGenerator, embed_n, etkdg_baseline, featurize_mol

RESULTS = []


def check(name, fn):
    try:
        info = fn()
        RESULTS.append((name, 'PASS', info or ''))
    except Exception as e:  # noqa: BLE001
        RESULTS.append((name, 'FAIL', f'{type(e).__name__}: {e}'[:400]))


def xyz(mols):
    return [m.GetConformer().GetPositions() for m in mols]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model_dir', required=True)
    ap.add_argument('--states_csv', required=True)
    a = ap.parse_args()
    gen = TDGenerator(a.model_dir)
    table = pd.read_csv(a.states_csv, comment='#').fillna('')
    built = {}
    for r in table.itertuples():
        try:
            built[r.id] = st.build_state(r.state_smiles, r.parent_smiles or None)
        except st.StateError:
            pass
    T = len(gen.types)
    q0 = T + 2 + 8 + 6 + 8  # featurize_mol layout: types | Z, aromatic | degree(8) | hybrid(6) | implicit valence(8) | charge(4)

    def dims():
        bad = {k: featurize_mol(m, gen.types).x.shape for k, (m, _) in built.items() if not gen.compatibility(m)[0]
               and featurize_mol(m, gen.types).x.shape != (m.GetNumAtoms(), gen.args.in_node_features)}
        assert not bad, bad
        return f'{sum(1 for m, _ in built.values() if not gen.compatibility(m)[0])} states: x = (n_atoms incl. H, {gen.args.in_node_features})'
    check('featuriser width == checkpoint in_node_features; one node per atom incl. explicit H', dims)

    def charges():
        m, _ = built['glycine_zwitterion']
        x = featurize_mol(m, gen.types).x[:, q0:q0 + 4].argmax(1).tolist()
        want = [{-1: 0, 0: 1, 1: 2}.get(at.GetFormalCharge(), 3) for at in m.GetAtoms()]
        assert x == want, (x, want)
        return 'N+ -> bucket +1, O- -> bucket -1, all others 0'
    check('formal-charge one-hot matches the prepared state (glycine zwitterion)', charges)

    def compat():
        hard, _ = gen.compatibility(Chem.AddHs(Chem.MolFromSmiles('[NH3+][C@@H](CS)C(=O)[O-]')))
        assert hard and 'S' in hard[0], hard
        _, ood = gen.compatibility(built['ibuprofen_S_anion'][0])
        assert any('heavy atoms' in o for o in ood), ood
        return f'S -> hard: {hard[0][:60]}...; ibuprofen -> OOD size'
    check('compatibility: unsupported element is hard failure, oversize is OOD', compat)

    def torsions():
        got = {k: int(gen.torsion_data(built[k][0])[0].edge_mask.sum())
               for k in ['phenol', 'phenolate', 'hydroxypyridine', 'pyridone', 'pyridinium', 'glycine_neutral',
                         'glycine_zwitterion']}
        assert (got['phenol'], got['phenolate'], got['hydroxypyridine'], got['pyridone'], got['pyridinium']) == (1, 0, 1, 0, 0), got
        return str(got)
    check('state-dependent TD torsion count (C-OH rotatable, C-O(-) not)', torsions)

    def repro():
        m = built['glycine_zwitterion'][0]
        A, _ = gen.generate(m, 4, 0); B, _ = gen.generate(m, 4, 0); C, _ = gen.generate(m, 4, 1)
        same = max(float(np.abs(x - y).max()) for x, y in zip(xyz(A), xyz(B)))
        diff = max(float(np.abs(x - y).max()) for x, y in zip(xyz(A), xyz(C)))
        assert same < 1e-4 and diff > 1e-2, (same, diff)
        return f'seed 0 twice: max |dx| = {same:.1e} A; seed 0 vs 1: {diff:.2f} A'
    check('TD generation reproducible for a fixed seed (same process), differs across seeds', repro)

    def distinct_seeds():
        m = built['ibuprofen_S_neutral'][0]
        out = []
        for seed in (0, 1):
            e = xyz(etkdg_baseline(m, 4, seed))
            t = [c.GetPositions() for c in embed_n(Chem.Mol(m), 4, seed).GetConformers()]
            for name, X in (('etkdg', e), ('td_seed', t)):
                dmin = min(float(np.abs(X[i] - X[j]).max()) for i in range(4) for j in range(i))
                assert dmin > 1e-2, (seed, name, dmin)
            out.append(f'seed {seed}: 4 distinct ETKDGv3 and 4 distinct TD seed conformers')
        return '; '.join(out)
    check('conformers within one call are distinct for seed 0 and seed 1 (ETKDG baseline and TD seeds)', distinct_seeds)

    def lock():
        m, info = built['acac_Z_enol']
        L, gi = gen.generate(m, 8, 0, lock_stereo_double_bonds=True)
        U, _ = gen.generate(m, 8, 0)
        ok_l = [va.canon_from_3d(x) == info['state_smiles_canonical'] for x in L]
        ok_u = [va.canon_from_3d(x) == info['state_smiles_canonical'] for x in U]
        assert gi['locked_double_bonds'] and all(ok_l), (gi['locked_double_bonds'], ok_l)
        return f'locked: {sum(ok_l)}/8 correct E/Z; unlocked TD: {sum(ok_u)}/8 correct E/Z'
    check('E/Z lock keeps the specified double-bond configuration', lock)

    def sdf():
        n = 0
        with tempfile.TemporaryDirectory() as d:
            for k, (m, info) in built.items():
                confs = etkdg_baseline(m, 1, 0)
                p = os.path.join(d, f'{k}.sdf')
                write_sdf(p, confs, [dict(state_id=k, method='etkdg_v3', state_smiles=info['state_smiles_canonical'])])
                rt = next(iter(Chem.SDMolSupplier(p, removeHs=False)))
                r = va.roundtrip(rt, info['state_smiles_canonical'])
                assert r['rt_atoms'] and r['rt_bonds'] and r['rt_isomeric'], (k, r)
                assert [x.GetAtomicNum() for x in rt.GetAtoms()] == [x.GetAtomicNum() for x in m.GetAtoms()], k
                assert Chem.GetFormalCharge(rt) == info['net_charge'] and rt.GetProp('state_id') == k, k
                n += 1
        return f'{n} states: atom order, charges, H, bond types, isomeric SMILES and SD properties preserved'
    check('SDF write/read round trip for every accepted smoke state', sdf)

    def geometry_controls():
        out = []
        for k in ['ibuprofen_S_neutral', 'acac_Z_enol']:
            m, info = built[k]
            c = etkdg_baseline(m, 1, 0)[0]
            tgt = info['state_smiles_canonical']
            base = va.check_conformer(c, tgt)
            assert base['valid'], (k, 'unperturbed ETKDG conformer flagged', base)
            pos = c.GetConformer().GetPositions()
            # 1 NaN
            x = Chem.Mol(c); x.GetConformer().SetAtomPosition(0, (float('nan'), 0., 0.))
            assert not va.check_conformer(x, tgt)['valid']
            # 2 stretched bond: move a terminal heavy atom 0.5 A further along its bond
            term = next(at for at in m.GetAtoms() if at.GetAtomicNum() > 1 and at.GetDegree() == 4 and
                        sum(nb.GetAtomicNum() > 1 for nb in at.GetNeighbors()) == 1)
            nb = next(n_.GetIdx() for n_ in term.GetNeighbors() if n_.GetAtomicNum() > 1)
            v = pos[term.GetIdx()] - pos[nb]; v /= np.linalg.norm(v)
            x = Chem.Mol(c); conf = x.GetConformer()
            moved = [term.GetIdx()] + [h.GetIdx() for h in term.GetNeighbors() if h.GetAtomicNum() == 1]
            for i in moved:
                conf.SetAtomPosition(i, (pos[i] + 0.5 * v).tolist())
            r2 = va.check_conformer(x, tgt)
            assert not r2['bonds_angles'], r2['max_12_dev']
            # 3 broken bond: 1.5 A -> connectivity must change
            x = Chem.Mol(c); conf = x.GetConformer()
            for i in moved:
                conf.SetAtomPosition(i, (pos[i] + 1.5 * v).tolist())
            r3 = va.check_conformer(x, tgt)
            assert r3['connectivity'] is False, r3.get('connectivity')
            # 4 clash: put a far atom 0.6 A from another
            topo = Chem.GetDistanceMatrix(m)
            i, j = map(int, np.unravel_index(np.argmax(topo), topo.shape))
            x = Chem.Mol(c); x.GetConformer().SetAtomPosition(j, (pos[i] + np.array([0.6, 0, 0])).tolist())
            assert not va.check_conformer(x, tgt)['clash']
            # 5 mirror image -> stereo must fail for the chiral molecule / still pass for achiral E/Z
            x = Chem.Mol(c); conf = x.GetConformer()
            for i in range(x.GetNumAtoms()):
                conf.SetAtomPosition(i, (pos[i] * np.array([-1, 1, 1])).tolist())
            rm = va.check_conformer(x, tgt)
            assert rm['stereo'] is (k != 'ibuprofen_S_neutral'), (k, rm['stereo_smiles_3d'])
            out.append(f'{k}: base valid, NaN/stretch(+0.5 A, dev {r2["max_12_dev"]:.2f})/break(+1.5 A)/clash/mirror detected')
        # 6 E/Z flip of the enol C=C
        m, info = built['acac_Z_enol']
        c = etkdg_baseline(m, 1, 0)[0]
        db = next(b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE and not b.IsInRing()
                  and b.GetStereo() != Chem.BondStereo.STEREONONE)
        i, j = db.GetBeginAtomIdx(), db.GetEndAtomIdx()
        a_ = next(n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() != j and n.GetAtomicNum() > 1)
        d_ = next(n.GetIdx() for n in m.GetAtomWithIdx(j).GetNeighbors() if n.GetIdx() != i and n.GetAtomicNum() > 1)
        conf = c.GetConformer()
        rdMolTransforms.SetDihedralDeg(conf, a_, i, j, d_, rdMolTransforms.GetDihedralDeg(conf, a_, i, j, d_) + 180)
        rf = va.check_conformer(c, info['state_smiles_canonical'])
        assert not rf['stereo'], rf['stereo_smiles_3d']
        out.append('Z-enol C=C rotated 180 deg -> stereo failure detected')
        return '; '.join(out)
    check('geometry validation detects NaN, stretched/broken bonds, clashes, mirror image, E/Z flip', geometry_controls)

    def hbond():
        m, info = built['acac_Z_enol']
        c = etkdg_baseline(m, 1, 0)[0]
        tgt = info['state_smiles_canonical']
        o_h = next(a for a in m.GetAtoms() if a.GetSymbol() == 'O' and a.GetTotalNumHs(includeNeighbors=True) == 1)
        h = next(n.GetIdx() for n in o_h.GetNeighbors() if n.GetAtomicNum() == 1)
        c2 = next(n for n in o_h.GetNeighbors() if n.GetAtomicNum() == 6)
        c3 = next(b.GetOtherAtom(c2) for b in c2.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
        o_acc = next(a for a in m.GetAtoms() if a.GetSymbol() == 'O' and a.GetTotalNumHs(includeNeighbors=True) == 0)
        c4 = next(n for n in o_acc.GetNeighbors())
        o_acc = o_acc.GetIdx()
        conf = c.GetConformer()
        # H-bonded geometry: C=O turned towards the enol O (C2=C3-C4=O cis), O-H turned towards it (H-O-C2=C3 cis)
        rdMolTransforms.SetDihedralDeg(conf, c2.GetIdx(), c3.GetIdx(), c4.GetIdx(), o_acc, 0.0)
        rdMolTransforms.SetDihedralDeg(conf, h, o_h.GetIdx(), c2.GetIdx(), c3.GetIdx(), 0.0)
        p = conf.GetPositions()
        dho = float(np.linalg.norm(p[h] - p[o_acc]))
        out = []
        for d_target in (None, 1.7, 1.2):  # as built; H moved to 1.7 A and 1.2 A from the acceptor (clash check only)
            x = Chem.Mol(c)
            if d_target:
                v = p[h] - p[o_acc]
                x.GetConformer().SetAtomPosition(h, (p[o_acc] + d_target * v / np.linalg.norm(v)).tolist())
            r = va.check_conformer(x, tgt)
            if d_target == 1.2:
                assert not r['clash'], (d_target, r['clash_pairs'])
            else:
                assert r['clash'], (d_target or dho, r['clash_pairs'])
                assert r['n_hbond_contacts'] >= (1 if (d_target or dho) < 1.9 else 0), (d_target or dho, r)
            out.append(f"H...O {d_target or dho:.2f} A: clash-free={r['clash']}, hbond contacts={r['n_hbond_contacts']}")
        return '; '.join(out)
    check('clash screen: intramolecular H-bond accepted, H...O < 1.5 A rejected', hbond)

    n_fail = sum(r[1] == 'FAIL' for r in RESULTS)
    for r in RESULTS:
        print(f'{r[1]}  {r[0]}' + (f'  -> {r[2]}' if r[2] else ''))
    print(f'test_backend: {len(RESULTS) - n_fail} passed, {n_fail} failed')
    raise SystemExit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
