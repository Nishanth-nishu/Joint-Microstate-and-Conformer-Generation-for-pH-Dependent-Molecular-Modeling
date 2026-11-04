"""Tests for the joint-prototype feature paths (jmcg venv). Run: python -m jmcg.joint.test_joint"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from jmcg.joint import features as ft
from jmcg.joint import p1_unmin as p1

RESULTS = []


def check(name, fn):
    try:
        RESULTS.append((name, 'PASS', fn() or ''))
    except Exception as e:  # noqa: BLE001
        RESULTS.append((name, 'FAIL', f'{type(e).__name__}: {e}'[:300]))


def main():
    def default_is_mmff():
        a, b = ft.ensemble_features('CC(=O)CC(C)=O', n_confs=6), ft.ensemble_features('CC(=O)CC(C)=O', n_confs=6, path='mmff')
        assert a['n_valid'] == b['n_valid'] and all(a[k] == b[k] for k in ('phi1', 'phi2', 'phi3')), (a, b)
        return f"acac n_valid {a['n_valid']}, phi3 {a['phi3']:.4f}"
    check('default feature path is the MMFF control and is deterministic', default_is_mmff)

    def bad_path():
        try:
            ft.ensemble_features('CCO', path='xtb')
        except ValueError:
            return 'ValueError'
        raise AssertionError('no error for unknown path')
    check('unknown feature path is rejected', bad_path)

    def unmin_uniform():
        r = ft.ensemble_features('CCCO', n_confs=8, path='unmin', return_ensemble=True)
        m, keep, w = r['ensemble']
        assert r['n_valid'] == len(keep) > 0 and np.allclose(w, 1 / len(keep))
        assert abs(r['phi1'] + np.log(len(keep))) < 1e-12, r['phi1']
        # the unminimised path must not touch coordinates: the kept conformers are the ETKDG embeddings themselves
        seed = ft.conformer_seeds(0, 8)[0]
        t = Chem.AddHs(Chem.MolFromSmiles('CCCO'))
        p = AllChem.ETKDGv3(); p.randomSeed = seed
        AllChem.EmbedMolecule(t, p)
        assert np.allclose(t.GetConformer().GetPositions(), m.GetConformer(0).GetPositions())
        return f"n_valid {r['n_valid']}/8, phi1 = -ln n_valid"
    check('unmin path: uniform weights, phi1 = -ln n_valid, coordinates = ETKDG embeddings', unmin_uniform)

    def classify():
        c = {s: p1.classify(s)['category'] for s in ('NCC(=O)O', '[NH3+]CC(=O)[O-]', 'CC(=O)[O-]', '[O-]c1ccccc1',
                                                     'CC[NH3+]', 'C[S-]')}
        assert c == {'NCC(=O)O': 'neutral', '[NH3+]CC(=O)[O-]': 'zwitterion', 'CC(=O)[O-]': 'carboxylate',
                     '[O-]c1ccccc1': 'phenolate', 'CC[NH3+]': 'cation', 'C[S-]': 'other_anion'}, c
        return 'neutral/zwitterion/carboxylate/phenolate/cation/other_anion'
    check('state categories', classify)

    n_fail = sum(r[1] == 'FAIL' for r in RESULTS)
    for r in RESULTS:
        print(f'{r[1]}  {r[0]}' + (f'  -> {r[2]}' if r[2] else ''))
    print(f'test_joint: {len(RESULTS) - n_fail} passed, {n_fail} failed')
    raise SystemExit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
