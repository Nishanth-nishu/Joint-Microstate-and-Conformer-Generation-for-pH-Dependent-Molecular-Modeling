"""Figures for the top-level README from committed results (research/experiments/N1_pilot/run2)."""
import csv
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'   # reference palette slots 1-3 (validated instance)
STRATA = [('sensitive', BLUE), ('control', ORANGE), ('other', AQUA)]

plt.rcParams.update({'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
                     'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'text.color': INK, 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8, 'axes.axisbelow': True})


def fig_gate(an, out):
    names = ['novartis_acid', 'novartis_base', 'novartis', 'sampl6', 'sampl7', 'sampl8']
    labels = ['Novartis\nacid', 'Novartis\nbase', 'Novartis\ntotal', 'SAMPL6', 'SAMPL7', 'SAMPL8']
    g = an['sanity_gate']
    rep = [np.mean(g[n]['reported']) if isinstance(g[n]['reported'], list) else g[n]['reported'] for n in names]
    ours = [g[n]['ours'] for n in names]
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(7, 3.4))
    for xi, r in zip(x, rep):
        ax.fill_between([xi - 0.3, xi + 0.3], r - 0.1, r + 0.1, color=GRID, lw=0)
    ax.scatter(x, rep, s=70, marker='s', color=BLUE, label='Uni-pKa paper, Table 1', zorder=3, edgecolor=SURFACE, lw=2)
    ax.scatter(x, ours, s=70, marker='o', color=ORANGE, label='Our re-run (5-fold mean)', zorder=4, edgecolor=SURFACE, lw=2)
    ax.set_xticks(x, labels)
    ax.set_ylabel('Macro-pKa MAE (pKa units)')
    ax.set_ylim(0.3, 1.0)
    ax.set_title('Pipeline reproduces published Uni-pKa accuracy (grey band = ±0.10 gate)', loc='left', fontsize=11)
    ax.legend(frameon=False, loc='upper right')
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def fig_forest(an, out):
    rows = []
    for s, _ in STRATA:
        st = an['strata'][s]
        rows.append((f'{s} (n={st["n"]}): 2D − 3D', st['delta_2D_minus_3Donly']))
    rows.append(('DiD: sensitive − control', an['did_sensitive_minus_control']))
    st = an['strata']['sensitive']
    rows += [('sensitive: Boltzmann(LSE) − mean', st['delta_lse_minus_mean']),
             ('sensitive: min − mean', st['delta_min_minus_mean']),
             ('sensitive: 1 conformer − 10', st['delta_single_minus_mean'])]
    fig, ax = plt.subplots(figsize=(8.6, 3.8))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, (d, lo, hi)) in zip(y, rows):
        ax.plot([lo, hi], [yi, yi], color=BLUE, lw=2, solid_capstyle='round')
        ax.scatter([d], [yi], s=64, color=BLUE, zorder=3, edgecolor=SURFACE, lw=2)
        ax.text(1.02, yi, f'{d:+.3f} [{lo:+.3f}, {hi:+.3f}]', va='center', fontsize=8.5, color=INK2,
                transform=ax.get_yaxis_transform(), clip_on=False)
    ax.axvline(0, color=INK2, lw=1)
    ax.axvspan(0.05, 0.2, color=GRID, alpha=0.6, lw=0)
    ax.text(0.052, y[0] + 0.55, 'effect size expected\nbefore running', fontsize=8, color=INK2, va='bottom')
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlim(-0.07, 0.2)
    ax.set_ylim(-0.7, len(rows) + 0.4)
    ax.set_xlabel('Δ MAE in pKa units, 95 % paired bootstrap CI (> 0: 3D helps)')
    ax.set_title('N1 falsified: 3D conformers do not change Uni-pKa error', loc='left', fontsize=11, x=-0.45)
    ax.grid(axis='y', visible=False)
    fig.tight_layout(rect=(0, 0, 0.8, 1)); fig.savefig(out, dpi=160); plt.close(fig)


def fig_scatter(transitions, out):
    rows = list(csv.DictReader(open(transitions)))
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.8))
    ax = axes[0]
    for s, c in STRATA:
        x = [float(r['A3D-only-mean']) for r in rows if r['stratum_v2'] == s]
        y = [float(r['A2D']) for r in rows if r['stratum_v2'] == s]
        ax.scatter(x, y, s=18, color=c, alpha=0.8, label=f'{s} (n={len(x)})', edgecolor=SURFACE, lw=0.6)
    lim = [-1, 14]
    ax.plot(lim, lim, color=INK2, lw=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('Predicted pKa, 10 3D conformers'); ax.set_ylabel('Predicted pKa, 2D coordinates only')
    ax.set_title('2D vs 3D input: same prediction', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    ax = axes[1]
    data = [[float(r['conf_pka_std']) for r in rows if r['stratum_v2'] == s] for s, _ in STRATA]
    bins = np.linspace(0, 0.2, 21)
    for (s, c), d in zip(STRATA, data):
        ax.hist(d, bins=bins, histtype='step', lw=2, color=c, label=f'{s} (median {np.median(d):.3f})')
    ax.set_xlabel('SD of predicted pKa across 10 conformers'); ax.set_ylabel('Transitions')
    ax.set_title('Conformer-to-conformer spread is tiny', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else '.'
    run2 = os.path.join(root, 'research/experiments/N1_pilot/run2')
    figs = os.path.join(root, 'research/figures')
    os.makedirs(figs, exist_ok=True)
    an = json.load(open(os.path.join(run2, 'analysis_gate_table1.json')))
    fig_gate(an, os.path.join(figs, 'n1_sanity_gate.png'))
    fig_forest(an, os.path.join(figs, 'n1_effects_forest.png'))
    fig_scatter(os.path.join(run2, 'transitions.csv'), os.path.join(figs, 'n1_conformer_invariance.png'))
    print('figures written to', figs)


if __name__ == '__main__':
    main()
