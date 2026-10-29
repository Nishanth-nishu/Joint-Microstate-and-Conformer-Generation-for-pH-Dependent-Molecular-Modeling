"""State scoring for the sequential baseline: released Uni-pKa microstate free energies and the populations and
macro-pKa they imply. Thermodynamic convention (stated, not assumed physical):

  G_z      Uni-pKa free-energy head output for microstate z, dimensionless (units of kT), the "standard free energy
           before adding the pH-dependent m ln(10) pH term" (unipkainfer docstring). Per fold: mean over the released
           inference protocol's 11 coordinate sets (10 ETKDG+MMFF + 1 2D; N1 patches applied by importing n1_pilot).
  q_z      net formal charge; within one enumerated ensemble states differ only by protons, so q_z - m_z is constant
           and using q in place of the proton count m only shifts every logit by the same amount.
  logit_z(pH) = -G_z - q_z ln(10) (pH - PKA_MEAN),   p(z | G, pH) = softmax over the ENUMERATED ensemble.
  macro-pKa between charge levels q and q-1: (LSE_{z:q}(-G) - LSE_{z:q-1}(-G)) / ln 10 + PKA_MEAN, which is exactly
  the pH where the two levels have equal population under the logits above (and Uni-pKa's own pKa formula).
  PKA_MEAN = 6.504894871171601 is Uni-pKa's training-label mean (added back at inference).
These are MODEL populations over an enumerated set at the implicit T and solvent of the training data (aqueous pKa,
~25 C); they are not validated equilibrium populations. Fold handling: populations are computed per fold and then
averaged (mixture of the 5 fold models); macro-pKa per fold then averaged (N1/T-pilot protocol). Fold spread is kept.
"""
import math
import time

import numpy as np
import pandas as pd
from rdkit import Chem

LN10 = math.log(10)
PKA_MEAN = 6.504894871171601


def lse(x):
    x = np.asarray(x, dtype=float)
    m = x.max()
    return float(m + np.log(np.exp(x - m).sum()))


def free_energies(smiles, model_dir, folds=(0, 1, 2, 3, 4), threads=8, task='seq'):
    """DataFrame(smiles, fold, G, n_coords) for every (unique) SMILES and fold. Requires the jmcg venv."""
    from jmcg.experiments import n1_pilot as n1  # applies the faithful-conformer and capture patches
    from unipkainfer.pka_predictor.free_energy import UnipkaFreeEnergyConfig, predict_standard_free_energies
    uniq = sorted(set(smiles))
    rows = []
    for fold in folds:
        t0 = time.time()
        n1._CONF.clear()
        cfg = UnipkaFreeEnergyConfig(model_dir=model_dir, folds=(fold,), nthreads=threads, gpu=False)
        predict_standard_free_energies([Chem.MolFromSmiles(s) for s in uniq], task_name=f'{task}_fold{fold}', config=cfg)
        conf = pd.concat(n1._CONF)
        conf['smiles_c'] = conf['smiles'].map(n1.canon)
        for smi, g in conf.groupby('smiles_c'):
            rows.append(dict(smiles=smi, fold=fold, G=n1.state_energies(g, 'all_mean'), n_coords=len(g)))
        print(f'fold {fold}: {len(uniq)} states in {time.time() - t0:.0f}s', flush=True)
    out = pd.DataFrame(rows)
    missing = set(map(lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)), uniq)) - set(out.smiles)
    if missing:
        raise RuntimeError(f'{len(missing)} states without a free energy, e.g. {sorted(missing)[:3]}')
    return out


def populations(G, q, ph):
    """p over one ensemble at one pH (arrays G, q in the same order)."""
    logit = -np.asarray(G, float) - np.asarray(q, float) * LN10 * (ph - PKA_MEAN)
    logit -= logit.max()
    p = np.exp(logit)
    return p / p.sum()


def macro_pkas(G, q):
    """{(q_acid, q_base): pKa} for every adjacent pair of charge levels present in the ensemble."""
    G, q = np.asarray(G, float), np.asarray(q)
    out = {}
    for qa in sorted(set(q.tolist())):
        if qa - 1 in set(q.tolist()):
            out[(qa, qa - 1)] = (lse(-G[q == qa]) - lse(-G[q == qa - 1])) / LN10 + PKA_MEAN
    return out


def log10_k(G1, G2):
    """Tautomer ratio log10([T2]/[T1]) from two same-charge states: (G_T1 - G_T2) / ln 10."""
    return (G1 - G2) / LN10
