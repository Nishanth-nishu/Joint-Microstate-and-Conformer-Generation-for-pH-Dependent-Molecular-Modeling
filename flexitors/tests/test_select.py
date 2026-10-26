"""Exp 5 selectors and pool-subset metrics (EXP-5 §9, §12.2-§12.4; EVAL §3.1; plan_P2 V-E5-5/7)."""
import inspect

import numpy as np
import pytest

from flexitors.eval.pool_metrics import subset_metrics
from flexitors.sampling import select as S


def data(n=40, m=3, seed=0):
    rng = np.random.default_rng(seed)
    return rng.uniform(-np.pi, np.pi, (n, m)), rng.normal(0, 5, n)


@pytest.mark.parametrize('method', S.SELECTORS)
def test_exact_size_subset_unique_deterministic(method):
    tors, e = data()
    a = S.select(method, tors, e, 20, np.random.default_rng(1))
    b = S.select(method, tors, e, 20, np.random.default_rng(1))
    assert len(a) == 20 and len(set(a.tolist())) == 20 and set(a.tolist()) <= set(range(40))
    assert np.array_equal(a, b)


def test_selector_semantics():
    tors, e = data()
    assert set(S.select('lowest_E', tors, e, 5, None).tolist()) == set(np.argsort(e)[:5].tolist())
    f = S.select('fps_lowE', tors, e, 5, np.random.default_rng(0))
    assert f[0] == int(np.argmin(e))
    keep = set(np.argsort(e)[:30].tolist())  # q = 0.25 of 40 proposals
    assert set(S.select('filter25_fps_lowE', tors, e, 20, np.random.default_rng(0)).tolist()) <= keep
    e2 = e.copy()
    e2[:10] = np.nan  # failed MMFF setups are +inf, i.e. filtered first
    assert not set(S.select('lowest_E', tors, e2, 10, None).tolist()) & set(range(10))
    assert 'filter50' not in ''.join(S.SELECTORS)  # degenerate at c = 2 (Q-P2-009)
    with pytest.raises(ValueError):
        S.select('random', tors, e, 41, np.random.default_rng(0))


def test_fps_spreads_points_on_the_circle():
    tors = np.linspace(-np.pi, np.pi, 36, endpoint=False)[:, None]
    sel = S.fps(tors, 4, 0)
    angles = np.sort(tors[sel, 0])
    assert np.allclose(np.diff(angles), np.pi / 2, atol=0.2)


def test_selectors_are_gt_free():
    sig = list(inspect.signature(S.select).parameters)
    assert sig == ['method', 'tors', 'energy', 'n_out', 'rng']
    src = inspect.getsource(S)
    assert 'rmsd' not in src.lower().replace('no rmsd', '') and 'true_mols' not in src and 'evaluate' not in src


def test_subset_metrics_match_evaluator_formulas_and_validity_secondary():
    rng = np.random.default_rng(3)
    R = rng.uniform(0, 1, (5, 12))
    cols = np.arange(6)
    m = subset_metrics(R, cols)
    assert m['AMR-R'] == pytest.approx(R[:, :6].min(1).mean())
    assert m['COV-R@0.5'] == pytest.approx((R[:, :6].min(1) < 0.5).mean())
    assert m['AMR-P'] == pytest.approx(R[:, :6].min(0).mean())
    inv = np.zeros(12, bool)
    inv[int(np.argmin(R[0, :6]))] = True
    v = subset_metrics(R, cols, inv)
    assert v['n_invalid'] == 1 and v['AMR-R_valid'] >= m['AMR-R']
    assert v['AMR-P_valid'] == pytest.approx(R[:, :6][:, ~inv[:6]].min(0).mean())
    assert v['COV-P@0.5_valid'] <= m['COV-P@0.5']
