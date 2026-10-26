"""experiments/exp3/analysis.py on synthetic arm outputs: the pre-declared family and non-inferiority (EXP-3 §9.3)."""
import importlib.util
import json
import pickle

import numpy as np
import pandas as pd
import pytest

from flexitors.td import REPO_ROOT


@pytest.fixture(scope='module')
def ana():
    spec = importlib.util.spec_from_file_location('exp3_analysis', REPO_ROOT / 'experiments' / 'exp3' / 'analysis.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_family_and_non_inferiority(ana, tmp_path):
    rng = np.random.default_rng(0)
    mols = [f'm{i}' for i in range(40)]
    for arm, ez_shift, amr_shift in (('O', 0, 0), ('K1', 0.1, 0.0), ('K2_inf', 0.05, 0.01), ('F_flip', 0.02, 0),
                                     ('F_filter', 0.03, 0), ('K1_wrong', -0.2, 0)):
        for s in (0, 1, 2):
            d = tmp_path / 'qm9' / f'{arm}_s{s}'
            d.mkdir(parents=True)
            st = pd.DataFrame(dict(corrected_smiles=mols, n_double=[1] * 30 + [0] * 10,
                                   ez_gt_consistency=np.r_[np.clip(rng.uniform(0.6, 0.9, 30) + ez_shift, 0, 1),
                                                           [np.nan] * 10],
                                   planarity_share=rng.uniform(0.7, 1.0, 40), gt_mixed=[False] * 40))
            st.to_csv(d / 'stereo_eval.csv', index=False)
            res = {(m, m): {'rmsd': rng.uniform(0.05, 0.6, (3, 6)) + amr_shift} for m in mols}
            pickle.dump({'results': res}, open(d / 'eval.pkl', 'wb'))
            (d / 'budget.json').write_text(json.dumps({'nfe': 100, 'extra_nfe': 5 if arm == 'F_filter' else 0}))
    C = ana.analyse(tmp_path, 'qm9', tmp_path / 'out')
    prim = C[C.family == 'primary']
    assert len(prim) == 8 and prim.p_holm.notna().all()
    ni = C[C.family == 'non-inferiority']
    assert len(ni) == 8 and ni.n.max() == 30
    assert set(C[C.family == 'secondary'].contrast) == {'K1_wrong-K1'}
    k1 = prim[(prim.contrast == 'K1-O') & (prim.metric == 'ez_gt_consistency')].iloc[0]
    assert k1['diff'] > 0.05
    assert 'K1_wrong missing' not in (tmp_path / 'out' / 'exp3_qm9_summary.txt').read_text()
