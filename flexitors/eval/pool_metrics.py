"""Metrics of a column subset of an evaluated pool (EVAL §2-§3.1; EXP-5 §12.1-§12.2, §12.7).

The pool is scored once with the unchanged ``evaluate_confs.py`` (rows = cleaned GT conformers, columns = proposals);
every arm is a column subset, and the metrics are the evaluator's (``evaluate_confs.py:39-45``): AMR-R = mean over
rows of the row minimum, COV-R@t = share of rows with minimum < t (strict), AMR-P = mean over columns of the column
minimum, COV-P@t likewise. Arms read a prefix or a subset of the pool and never its length (D-101).

Validity-adjusted secondary (EVAL §3.1 (b)): an invalid conformer keeps its slot; its recall RMSD is +∞, it counts as
uncovered in COV-P, and AMR-P is over valid conformers only. The primary leaves the evaluator unchanged.
"""
from __future__ import annotations

import numpy as np

THRESHOLDS = (0.5, 0.1)   # QM9 primary δ = 0.5 Å; COV-R@0.1 key secondary (EVAL §2)


def subset_metrics(rmsd: np.ndarray, cols: np.ndarray, invalid: np.ndarray | None = None) -> dict:
    """Metrics of columns ``cols`` of one molecule's RMSD matrix (n_true, N); ``invalid`` (N,) bool for the secondary.

    NaN rows (evaluator 'additional failures') propagate as in ``evaluate_confs.py``."""
    S = rmsd[:, cols]
    out = {'AMR-R': float(S.min(axis=1).mean()), 'AMR-P': float(S.min(axis=0).mean())}
    for t in THRESHOLDS:
        out[f'COV-R@{t}'] = float((S.min(axis=1) < t).mean())
        out[f'COV-P@{t}'] = float((S.min(axis=0) < t).mean())
    if invalid is not None:
        bad = np.asarray(invalid, dtype=bool)[cols]
        R = np.where(bad[None, :], np.inf, S)
        out['AMR-R_valid'] = float(R.min(axis=1).mean())
        out['COV-R@0.5_valid'] = float((R.min(axis=1) < 0.5).mean())
        out['COV-P@0.5_valid'] = float(np.mean((S.min(axis=0) < 0.5) & ~bad))
        out['AMR-P_valid'] = float(S[:, ~bad].min(axis=0).mean()) if (~bad).any() else np.nan
        out['n_invalid'] = int(bad.sum())
    return out
