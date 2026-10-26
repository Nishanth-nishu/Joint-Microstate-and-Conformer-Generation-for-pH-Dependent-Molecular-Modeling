"""Exp 3 Stage A analysis, PRE-DECLARED before any arm runs (EXP-3 §9.1, §9.3, §9.4; Q-P1-006, Q-P2-007, Q-P3-004
answers). Commit time = pre-registration.

Layout: ``<root>/<dataset>/<arm>_s<seed>/`` holds the spliced ``eval.pkl`` (unchanged evaluator), ``stereo_eval.csv``
(``flexitors.eval.stereo_eval``) and ``budget.json`` / ``meta.json``. Arms: O, K1, K1_wrong, K1_canonical, K2_inf
(sampler), F_flip, F_filter (post hoc on O); QM9 3 sampling seeds, DRUGS seed 0 first (§9.4).

Endpoints per dataset (§9.3):
  E/Z endpoint  QM9: share of conformers consistent with the GT ensemble (``ez_gt_consistency``) on the structural
                molecules (37); DRUGS: E/Z fidelity to the SMILES target (``ez_fidelity``) on the 126 specified.
  planarity     share of (conformer, bond) pairs with |δ_b| ≤ 20° on every TD-rotatable double bond (QM9: the 102).
Primary family (Holm over 8, per dataset): {K1, K2_inf, F_flip, F_filter} − O × {E/Z endpoint, planarity}; paired
molecule bootstrap (5000) on the per-molecule seed mean.
Non-inferiority (co-primary): the same four contrasts on the AFFECTED subset (molecules with a TD-rotatable double
bond): AMR-R upper 95 % bound ≤ +0.004 Å and COV-R@0.5 lower bound ≥ −0.5 pt, hierarchical bootstrap (seeds jointly,
then molecules). The all-molecule difference is reported as derived: (n_affected / N) × subset difference (§9.1).
Secondary: K1_wrong − K1 and K1_canonical − K1 on both endpoints; F_filter extra NFE, replacements, fills. Sensitivity:
molecules whose GT ensemble mixes E/Z states excluded (QM9). If K1_wrong is missing, the report states that K1's effect
cannot be called "about stereo" (Q-P3-004).
"""
from __future__ import annotations

import json
import pickle
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd

from flexitors.provenance import write_provenance
from flexitors.td import import_baseline

N_BOOT = 5000
AMR_MARGIN = 0.004       # Å (EXP-3 §9.3)
COV_MARGIN = 0.5         # points (EXP-3 §9.3)
PRIMARY_ARMS = ['K1', 'K2_inf', 'F_flip', 'F_filter']
SECONDARY = [('K1_wrong', 'K1'), ('K1_canonical', 'K1')]
EZ_ENDPOINT = {'qm9': 'ez_gt_consistency', 'drugs': 'ez_fidelity'}


def recall(path: Path) -> pd.DataFrame:
    """Per-molecule AMR-R and COV-R@0.5 of one spliced eval.pkl (TD convention handled by the universe below)."""
    with open(path, 'rb') as f:
        R = pickle.load(f)['results']
    rows = [dict(corrected_smiles=c, **{'AMR-R': float(np.min(r['rmsd'], 1).mean()),
                                        'COV-R@0.5': 100 * float((np.min(r['rmsd'], 1) < 0.5).mean())})
            for (s, c), r in R.items()]
    return pd.DataFrame(rows).set_index('corrected_smiles')


def load(root: Path, dataset: str) -> tuple[dict, dict, dict]:
    """{(arm, seed): stereo frame}, {(arm, seed): recall frame}, {(arm, seed): budget}."""
    stereo, rec, bud = {}, {}, {}
    for d in sorted((root / dataset).glob('*_s*')):
        arm, seed = d.name.rsplit('_s', 1)
        f = d / 'stereo_eval.csv'
        if f.exists() and f.stat().st_size > 1:  # empty when the arm's subset has no double bond
            stereo[(arm, int(seed))] = pd.read_csv(f).set_index('corrected_smiles')
        if (d / 'eval.pkl').exists():
            rec[(arm, int(seed))] = recall(d / 'eval.pkl')
        if (d / 'budget.json').exists():
            bud[(arm, int(seed))] = json.loads((d / 'budget.json').read_text())
    return stereo, rec, bud


def per_seed(frames: dict, arm: str, col: str, mols) -> np.ndarray:
    """(seeds, molecules) array of one arm and column on ``mols``; NaN where missing."""
    seeds = sorted(s for a, s in frames if a == arm)
    return np.stack([frames[(arm, s)][col].reindex(mols).values.astype(float) for s in seeds])


def paired(frames: dict, a: str, b: str, col: str, mols, rng) -> dict:
    """Paired molecule bootstrap of a − b on the per-molecule seed mean."""
    alf = import_baseline('analyze_leakfree')
    d = np.nanmean(per_seed(frames, a, col, mols), 0) - np.nanmean(per_seed(frames, b, col, mols), 0)
    d = d[np.isfinite(d)]
    mean, lo, hi, p = alf.boot(d, rng, N_BOOT)
    return dict(contrast=f'{a}-{b}', metric=col, n=len(d), diff=mean, ci_lo=lo, ci_hi=hi, p_boot=p)


def hierarchical(frames: dict, a: str, b: str, col: str, mols, rng) -> tuple[float, float]:
    """Seeds jointly (arms share seeds), then molecules; 95 % percentile interval of a − b."""
    alf = import_baseline('analyze_leakfree')
    D = per_seed(frames, a, col, mols) - per_seed(frames, b, col, mols)
    D = D[:, np.isfinite(D).all(0)]
    return alf.hier_ci(D, None, rng, N_BOOT)


def analyse(root: Path, dataset: str, out: Path) -> pd.DataFrame:
    """The pre-declared family, non-inferiority and secondary contrasts of one dataset."""
    rng = np.random.default_rng(0)
    alf = import_baseline('analyze_leakfree')
    stereo, rec, bud = load(root, dataset)
    arms = {a for a, _ in stereo}
    O = stereo[('O', 0)]
    ez = EZ_ENDPOINT[dataset]
    ez_mols = O.index[O[ez].notna()] if ez in O else O.index[:0]
    plan_mols = O.index[O.n_double > 0]
    rows = []
    prim = [paired(stereo, a, 'O', col, mols, rng) for a in PRIMARY_ARMS if a in arms
            for col, mols in ((ez, ez_mols), ('planarity_share', plan_mols))]
    for c, p in zip(prim, alf.holm([c['p_boot'] for c in prim]) if prim else []):
        rows.append({**c, 'family': 'primary', 'p_holm': p})
    affected = plan_mols
    n_all = len(next(iter(rec.values()))) if rec else 0
    for a in [x for x in PRIMARY_ARMS if x in {k for k, _ in rec}]:
        for col, margin, side in (('AMR-R', AMR_MARGIN, 'upper'), ('COV-R@0.5', -COV_MARGIN, 'lower')):
            c = paired(rec, a, 'O', col, affected, rng)
            lo, hi = hierarchical(rec, a, 'O', col, affected, rng)
            ok = hi <= margin if side == 'upper' else lo >= margin
            rows.append({**c, 'family': 'non-inferiority', 'hier_lo': lo, 'hier_hi': hi, 'margin': margin,
                         'non_inferior': bool(ok), 'derived_all_molecule_diff': c['diff'] * len(affected) / max(n_all, 1)})
    for a, b in SECONDARY:
        if a in arms and b in arms:
            rows += [{**paired(stereo, a, b, col, mols, rng), 'family': 'secondary'}
                     for col, mols in ((ez, ez_mols), ('planarity_share', plan_mols))]
    if 'gt_mixed' in O:
        clean = O.index[(O.gt_mixed != True) & O[ez].notna()] if ez in O else O.index[:0]  # noqa: E712
        rows += [{**paired(stereo, a, 'O', ez, clean, rng), 'family': 'sensitivity (no mixed-GT molecules)'}
                 for a in PRIMARY_ARMS if a in arms]
    C = pd.DataFrame(rows)
    out.mkdir(parents=True, exist_ok=True)
    C.to_csv(out / f'exp3_{dataset}_contrasts.csv', index=False)
    extra = {f'{a}_s{s}': {k: b.get(k) for k in ('nfe', 'extra_nfe')} for (a, s), b in bud.items()}
    notes = [] if 'K1_wrong' in arms else ['K1_wrong missing: no claim that K1\'s effect is "about stereo" (Q-P3-004)']
    txt = [f'# Exp 3 {dataset}: arms {sorted(arms)}; E/Z endpoint {ez} on {len(ez_mols)} molecules; planarity on '
           f'{len(plan_mols)}', *notes, json.dumps(extra), '', C.round(5).to_string(index=False)]
    (out / f'exp3_{dataset}_summary.txt').write_text('\n'.join(txt) + '\n')
    write_provenance(out)
    print(txt[0])
    return C


def main() -> None:
    ap = ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--root', required=True, help='results/round5/exp3')
    ap.add_argument('--dataset', required=True, choices=sorted(EZ_ENDPOINT))
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    analyse(Path(a.root), a.dataset, Path(a.out))


if __name__ == '__main__':
    main()
