"""Exp 5 Stage A analysis (selection at matched compute), PRE-DECLARED before any pool exists (EXP-5 §12; D-101, D-102;
Q-P1-008/009, Q-P2-008/009, Q-P3-005/006 answers). Commit time = pre-registration (tools/analyze_round2.py pattern).

Pools (per ``generate.py``, per-conformer streams): test seed 0 TD = max(128, 4K_i) columns; test seeds 1, 2 TD = 4K_i;
test PG2K (2K_i, batch 128) and PG4K (4K_i, batch 128) per seed with the validation-chosen PG setting; validation seed
0 TD = 4K_i and PG2K for each of the four README settings. Each pool is scored once by the unchanged evaluator; arms are
column subsets (``flexitors.eval.pool_metrics``) and read prefixes, never the pool length.

Arms (§12.2), all output exactly 2K_i:
  P0  TD columns [0, 2K)                      P1  random 2K of TD [0, 4K)          (NFE control)
  P2  PG2K, all                               P2b random 2K of PG4K               (NFE-matched PG comparator)
  P3  selector on TD [0, 4K)                  P6  selector on PG4K                (composition)
  ORACLE  greedy min-AMR on P3's proposals (uses GT: analysis bound only)
Selector choice on validation (§12.4), rule fixed here: among the five candidates keep those with AMR-P ≤ AMR-P(random)
+ 0.004 Å and COV-P@0.5 ≥ COV-P(random) − 0.5 pt; take the lowest AMR-R, ties within 0.001 Å to the higher COV-R@0.1.
If none qualifies, or the winner's AMR-R gain over random is < 0.002 Å, or its paired 95 % CI includes 0, the
selector is "random" and P3/P6 are not run on test (K2-equivalent).
PG setting (§12.5): the same rule over the four README settings at c = 1 against TD P0 on validation (no futility
clause: P2b needs one setting; if none passes the precision margins, the lowest AMR-R is used and FLAGGED).
Test (§12.3): primary Holm family {P3 − P1, P3 − P2b} × {AMR-R, COV-R@0.5}, TD convention, paired molecule bootstrap
over the per-molecule mean of 3 sampling seeds (seed SD printed). Key secondary family, own Holm: the same contrasts on
COV-R@0.1, success intersection. Descriptive: P6 − P2b, P2 − P0, P2b − P1, P1 − P0 (sanity, expected ≈ 0), precision,
validity-adjusted metrics (EVAL §3.1 b). K3: the selector is a method only if P3 beats P2b on AMR-R (Holm) AND is
non-inferior on AMR-P (upper 95 % bound ≤ +0.004 Å) AND on COV-P@0.5 (lower bound ≥ −0.5 pt); else a baseline.
"""
from __future__ import annotations

import hashlib
import json
import pickle
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd

from flexitors.eval.pool_metrics import subset_metrics
from flexitors.provenance import write_provenance
from flexitors.sampling.seeding import numpy_rng
from flexitors.sampling.select import CANDIDATES, select
from flexitors.td import import_baseline

N_BOOT = 5000
AMRP_MARGIN = 0.004      # Å, round-2 margin [E-SAMPLE-057]
COVP_MARGIN = 0.005      # 0.5 pt, declared in EXP-5 §12.3
TIE = 0.001              # Å
FUTILITY = 0.002         # Å
METRICS = ['AMR-R', 'COV-R@0.5', 'COV-R@0.1', 'AMR-P', 'COV-P@0.5', 'AMR-R_valid', 'COV-R@0.5_valid', 'AMR-P_valid',
           'COV-P@0.5_valid']
PRIMARY = [('P3', 'P1'), ('P3', 'P2b')]
DESCRIPTIVE = [('P6', 'P2b'), ('P2', 'P0'), ('P2b', 'P1'), ('P1', 'P0')]


def load_pool(d: Path) -> dict:
    """{corrected: dict(rmsd, tors, energy, invalid)} of one pool directory (eval.pkl, features.pkl, validity.pkl)."""
    with open(d / 'eval.pkl', 'rb') as f:
        R = pickle.load(f)['results']
    feats = pickle.load(open(d / 'features.pkl', 'rb')) if (d / 'features.pkl').exists() else {}
    inval = pickle.load(open(d / 'validity.pkl', 'rb')) if (d / 'validity.pkl').exists() else {}
    return {corr: dict(rmsd=r['rmsd'], feat=feats.get(corr), invalid=inval.get(corr)) for (s, corr), r in R.items()}


def arm_columns(arm: str, K: int, pool: dict, seed: int, key: str, selector: str | None) -> np.ndarray:
    """Columns of the arm in its pool (prefix reads only; asserts the pool holds the columns it needs)."""
    need = {'P0': 2 * K, 'P2': 2 * K}.get(arm, 4 * K)
    if pool['rmsd'].shape[1] < need:
        raise AssertionError(f'{key}: pool has {pool["rmsd"].shape[1]} < {need} columns for {arm}')
    if arm in ('P0', 'P2'):
        return np.arange(2 * K)
    rng = numpy_rng(seed, key, 0, arm)
    if arm in ('P1', 'P2b') or selector == 'random':
        return np.sort(rng.choice(4 * K, 2 * K, replace=False))
    f = pool['feat']
    return select(selector, f['tors'][:4 * K], f['energy'][:4 * K], 2 * K, rng)


def oracle_columns(rmsd: np.ndarray, K: int) -> np.ndarray:
    """ORACLE greedy min-AMR selection of 2K of the first 4K columns (prototype select.py:20-24); analysis only."""
    P = rmsd[:, :4 * K]
    sel, cur = [], np.full(P.shape[0], np.inf)
    for _ in range(2 * K):
        cand = np.minimum(cur[:, None], P).mean(0)
        cand[sel] = np.inf
        j = int(np.argmin(cand))
        sel.append(j)
        cur = np.minimum(cur, P[:, j])
    return np.array(sel)


def arm_table(csv: pd.DataFrame, pools: dict[str, dict], arm: str, seed: int, selector: str | None) -> pd.DataFrame:
    """Per-molecule metrics of one arm and seed over every csv row (TD convention: a molecule absent from its pool is
    a failure: COV = 0, AMR = NaN)."""
    src = 'TD' if arm in ('P0', 'P1', 'P3', 'ORACLE') else ('PG2K' if arm == 'P2' else 'PG4K')
    rows = []
    for r in csv.itertuples(index=False):
        p = pools[src].get(r.corrected_smiles)
        row = dict(corrected_smiles=r.corrected_smiles, K=int(r.n_conformers), arm=arm, seed=seed)
        if p is None:
            rows.append({**row, 'failure': True, **{m: (np.nan if m.startswith('AMR') else 0.0) for m in METRICS}})
            continue
        K = int(r.n_conformers)
        cols = oracle_columns(p['rmsd'], K) if arm == 'ORACLE' else arm_columns(arm, K, p, seed, r.corrected_smiles,
                                                                                 selector)
        rows.append({**row, 'failure': False, **{m: np.nan for m in METRICS},
                     **subset_metrics(p['rmsd'], cols, p['invalid'])})
    return pd.DataFrame(rows)


def contrast(T: pd.DataFrame, a: str, b: str, metric: str, rng, universe=None) -> dict:
    """Paired molecule bootstrap of a − b on the per-molecule seed mean (TD convention)."""
    alf = import_baseline('analyze_leakfree')
    piv = T[T.arm.isin([a, b])].groupby(['corrected_smiles', 'arm'])[metric].mean().unstack()
    if universe is not None:
        piv = piv.loc[piv.index.intersection(universe)]
    d = (piv[a] - piv[b]).dropna()
    mean, lo, hi, p = alf.boot(d.values, rng, N_BOOT)
    seed_sd = T[T.arm.isin([a, b])].groupby(['seed', 'arm'])[metric].mean().unstack()
    sd = float((seed_sd[a] - seed_sd[b]).std(ddof=1)) if len(seed_sd) > 1 else np.nan
    return dict(contrast=f'{a}-{b}', metric=metric, n=len(d), diff=mean, ci_lo=lo, ci_hi=hi, p_boot=p, seed_sd=sd)


def choose(table: pd.DataFrame, ref: str, candidates: list[str], futility: bool, rng) -> dict:
    """The §12.4 rule on validation per-molecule metrics ``table`` (columns: arm, metrics)."""
    mean = table.groupby('arm')[['AMR-R', 'AMR-P', 'COV-P@0.5', 'COV-R@0.1']].mean()
    r = mean.loc[ref]
    ok = [c for c in candidates if mean.loc[c, 'AMR-P'] <= r['AMR-P'] + AMRP_MARGIN
          and mean.loc[c, 'COV-P@0.5'] >= r['COV-P@0.5'] - COVP_MARGIN]
    out = dict(reference=ref, means=mean.round(5).to_dict(orient='index'), eligible=ok, flagged=False)
    pool = ok or (candidates if not futility else [])
    if not ok and not futility:
        out['flagged'] = True
    if not pool:
        return {**out, 'choice': 'random', 'reason': 'no candidate passes the precision margins'}
    best = min(pool, key=lambda c: mean.loc[c, 'AMR-R'])
    ties = [c for c in pool if mean.loc[c, 'AMR-R'] - mean.loc[best, 'AMR-R'] <= TIE]
    best = max(ties, key=lambda c: mean.loc[c, 'COV-R@0.1'])
    if futility:
        c = contrast(table, best, ref, 'AMR-R', rng)
        out['winner_vs_ref'] = c
        if -c['diff'] < FUTILITY or c['ci_hi'] >= 0:
            return {**out, 'choice': 'random', 'reason': f'futility: gain {-c["diff"]:.4f} A, CI [{c["ci_lo"]:.4f}, '
                                                         f'{c["ci_hi"]:.4f}]'}
    return {**out, 'choice': best, 'reason': 'rule'}


def sha(paths: list[Path]) -> str:
    """Hash of the decision inputs (recorded with every decision file)."""
    h = hashlib.sha256()
    for p in sorted(paths):
        h.update(p.read_bytes())
    return h.hexdigest()


def validation(root: Path, csv: pd.DataFrame, out: Path) -> None:
    """Selector and PG-setting choices on the QM9 validation subset, written before any test pool is scored."""
    rng = np.random.default_rng(0)
    td = load_pool(root / 'val' / 'td_s0')
    tab = pd.concat([arm_table(csv, {'TD': td}, 'P1', 0, None).assign(arm='random')] +
                    [arm_table(csv, {'TD': td}, 'P3', 0, c).assign(arm=c) for c in CANDIDATES])
    sel = choose(tab, 'random', list(CANDIDATES), True, rng)
    p0 = arm_table(csv, {'TD': td}, 'P0', 0, None)
    pg_tabs = [p0]
    settings = sorted(d.name[len('pg_'):-len('_s0')] for d in (root / 'val').glob('pg_*_s0'))
    for s in settings:
        pg_tabs.append(arm_table(csv, {'PG2K': load_pool(root / 'val' / f'pg_{s}_s0')}, 'P2', 0, None).assign(arm=s))
    pg = choose(pd.concat(pg_tabs), 'P0', settings, False, rng) if settings else {'choice': None}
    inputs = list((root / 'val').rglob('eval.pkl')) + list((root / 'val').rglob('features.pkl'))
    dec = dict(selector=sel, pg=pg, inputs_sha256=sha(inputs))
    out.mkdir(parents=True, exist_ok=True)
    (out / 'exp5_validation_decision.json').write_text(json.dumps(dec, indent=1, default=float) + '\n')
    write_provenance(out)
    print('SUMMARY selector', sel['choice'], '| PG', pg['choice'])


def test(root: Path, csv: pd.DataFrame, out: Path, seeds=(0, 1, 2)) -> None:
    """The pre-declared test analysis; reads the hashed validation decision."""
    dec = json.loads((out / 'exp5_validation_decision.json').read_text())
    selector = dec['selector']['choice']
    tabs = []
    for s in seeds:
        pools = {'TD': load_pool(root / 'test' / f'td_s{s}'), 'PG2K': load_pool(root / 'test' / f'pg2k_s{s}'),
                 'PG4K': load_pool(root / 'test' / f'pg4k_s{s}')}
        arms = ['P0', 'P1', 'P2', 'P2b'] + (['P3', 'P6'] if selector != 'random' else []) + ['ORACLE']
        tabs += [arm_table(csv, pools, a, s, selector) for a in arms]
    T = pd.concat(tabs, ignore_index=True)
    T.to_csv(out / 'exp5_per_molecule.csv.gz', index=False, compression='gzip')
    rng = np.random.default_rng(0)
    alf = import_baseline('analyze_leakfree')
    rows = []
    if selector != 'random':
        prim = [contrast(T, a, b, m, rng) for a, b in PRIMARY for m in ('AMR-R', 'COV-R@0.5')]
        for c, p in zip(prim, alf.holm([c['p_boot'] for c in prim])):
            rows.append({**c, 'family': 'primary', 'p_holm': p})
        ok = T.groupby('corrected_smiles').failure.apply(lambda f: not f.any())
        inter = ok.index[ok]
        sec = [contrast(T, a, b, 'COV-R@0.1', rng, inter) for a, b in PRIMARY]
        for c, p in zip(sec, alf.holm([c['p_boot'] for c in sec])):
            rows.append({**c, 'family': 'key secondary', 'p_holm': p})
        for m in ('AMR-P', 'COV-P@0.5'):
            rows.append({**contrast(T, 'P3', 'P2b', m, rng), 'family': 'K3 non-inferiority'})
    pairs = [x for x in DESCRIPTIVE if selector != 'random' or 'P6' not in x]
    rows += [{**contrast(T, a, b, m, rng), 'family': 'descriptive'} for a, b in pairs for m in METRICS]
    C = pd.DataFrame(rows)
    C.to_csv(out / 'exp5_contrasts.csv', index=False)
    k3 = None
    if selector != 'random':
        g = C.set_index(['family', 'contrast', 'metric'])
        win = g.loc[('primary', 'P3-P2b', 'AMR-R')]
        k3 = bool(win['p_holm'] < 0.05 and win['diff'] < 0
                  and g.loc[('K3 non-inferiority', 'P3-P2b', 'AMR-P'), 'ci_hi'] <= AMRP_MARGIN
                  and g.loc[('K3 non-inferiority', 'P3-P2b', 'COV-P@0.5'), 'ci_lo'] >= -COVP_MARGIN)
    levels = T.groupby('arm')[METRICS + ['failure']].mean()
    txt = [f'# Exp 5 Stage A test analysis (selector {selector}; PG {dec["pg"]["choice"]})',
           f'K3 (selector is a method): {k3}', '', levels.round(4).to_string(), '', C.round(5).to_string(index=False)]
    (out / 'exp5_summary.txt').write_text('\n'.join(txt) + '\n')
    write_provenance(out)
    print('\n'.join(txt[:2]))


def main() -> None:
    ap = ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('mode', choices=['validation', 'test'])
    ap.add_argument('--root', required=True, help='results/round5/exp5 (val/, test/)')
    ap.add_argument('--csv', required=True, help='val_smiles.csv or test_smiles.csv')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    csv = pd.read_csv(a.csv)
    (validation if a.mode == 'validation' else test)(Path(a.root), csv, Path(a.out))


if __name__ == '__main__':
    main()
