"""experiments/exp5/analysis.py on synthetic pools: arm prefixes, decision rule, pre-declared families (EXP-5 §12)."""
import importlib.util
import json
import pickle

import numpy as np
import pandas as pd
import pytest

from flexitors.td import REPO_ROOT


@pytest.fixture(scope='module')
def ana():
    spec = importlib.util.spec_from_file_location('exp5_analysis', REPO_ROOT / 'experiments' / 'exp5' / 'analysis.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_pool(d, K, N, rng, missing=()):
    d.mkdir(parents=True)
    res, feats, inval = {}, {}, {}
    for k, n_k in K.items():
        if k in missing:
            continue
        n = N(n_k)
        res[(k, k)] = {'n_true': n_k, 'n_model': n, 'rmsd': rng.uniform(0.05, 1.0, (n_k, n))}
        feats[k] = dict(tors=rng.uniform(-np.pi, np.pi, (n, 2)), energy=rng.normal(0, 3, n), mmff_cpu_s=0.0)
        inval[k] = rng.random(n) < 0.05
    pickle.dump({'results': res}, open(d / 'eval.pkl', 'wb'))
    pickle.dump(feats, open(d / 'features.pkl', 'wb'))
    pickle.dump(inval, open(d / 'validity.pkl', 'wb'))


def test_validation_and_test_pipeline(ana, tmp_path):
    rng = np.random.default_rng(0)
    K = {f'm{i}': int(k) for i, k in enumerate(rng.integers(1, 6, 30))}
    csv = pd.DataFrame(dict(smiles=list(K), n_conformers=list(K.values()), corrected_smiles=list(K)))
    root = tmp_path / 'exp5'
    write_pool(root / 'val' / 'td_s0', K, lambda k: 4 * k, rng)
    for s in ('noninv_recall', 'inv_recall'):
        write_pool(root / 'val' / f'pg_{s}_s0', K, lambda k: 2 * k, rng)
    for s in (0, 1, 2):
        write_pool(root / 'test' / f'td_s{s}', K, (lambda k: max(128, 4 * k)) if s == 0 else (lambda k: 4 * k), rng,
                   missing=('m0',))
        write_pool(root / 'test' / f'pg2k_s{s}', K, lambda k: 2 * k, rng)
        write_pool(root / 'test' / f'pg4k_s{s}', K, lambda k: 4 * k, rng)
    out = tmp_path / 'ana'
    ana.validation(root, csv, out)
    dec = json.loads((out / 'exp5_validation_decision.json').read_text())
    assert dec['selector']['choice'] in ('random',) + tuple(ana.CANDIDATES) and dec['pg']['choice'] in (
        'noninv_recall', 'inv_recall')
    dec['selector']['choice'] = 'fps_lowE'  # force the full test path (random data rarely passes the futility rule)
    (out / 'exp5_validation_decision.json').write_text(json.dumps(dec))
    ana.test(root, csv, out)
    C = pd.read_csv(out / 'exp5_contrasts.csv')
    prim = C[C.family == 'primary']
    assert len(prim) == 4 and prim.p_holm.notna().all()
    assert set(prim.contrast) == {'P3-P1', 'P3-P2b'} and set(prim.metric) == {'AMR-R', 'COV-R@0.5'}
    assert len(C[C.family == 'key secondary']) == 2
    T = pd.read_csv(out / 'exp5_per_molecule.csv.gz')
    f = T[(T.corrected_smiles == 'm0') & T.arm.isin(['P0', 'P1', 'P3'])]
    assert f.failure.all() and (f['COV-R@0.5'] == 0).all()  # TD convention for a failed molecule
    assert 'K3 (selector is a method)' in (out / 'exp5_summary.txt').read_text()


def test_p0_is_the_pool_prefix_and_arms_need_4k(ana):
    rng = np.random.default_rng(1)
    pool = dict(rmsd=rng.uniform(size=(3, 12)), feat=dict(tors=rng.uniform(size=(12, 1)), energy=rng.normal(size=12)),
                invalid=None)
    assert ana.arm_columns('P0', 3, pool, 0, 'k', None).tolist() == list(range(6))
    assert set(ana.arm_columns('P1', 3, pool, 0, 'k', None).tolist()) <= set(range(12))
    with pytest.raises(AssertionError):
        ana.arm_columns('P3', 4, pool, 0, 'k', 'fps_lowE')
