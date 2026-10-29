"""Joint-generation prototype sampler: parent + pH -> (state, conformer) pairs from p_theta(z, x | G, h).

  level (net charge) free energies F_L = -LSE_{z in L}(-G_U(z))  (released Uni-pKa; keeps its macro-pKa exactly)
  within-level   p(z | L) = softmax_{z in L}(-(a G_U(z) + theta . phi(X_z)))   (E2 model; a, theta from e2 fit)
  pH             p(L | h) = softmax_L(-F_L - L ln10 (pH - PKA_MEAN))
  conformers     X_z ~ ETKDGv3 x30 + MMFF94s (validated), Boltzmann weights within z at 298.15 K
p(z, x_k | h) = p(L | h) p(z | L) w_z(x_k); sampled (z, x) pairs are written to sampled_pairs.sdf. With a = 1, theta = 0 this is the sequential baseline exactly (S0); with
theta = 0 it is the calibrated sequential baseline (S1). A charge level where any state lacks a valid conformer
ensemble (no phi) is scored with theta = 0 and listed in run_info.fallback_levels_theta0. This is a prototype: the within-level model was fitted on
neutral tautomer ratios only and is applied to every charge level; nothing here validates conformer populations.

Run (jmcg venv): python -m jmcg.joint.sample --state_dir <state_stage out> --coef <e2_result.json> --arm J1
                    --ph 7.4 --out <dir> [--n_pairs 20 --seed 0]
"""
import argparse
import json
import math
import os

import numpy as np
import pandas as pd
from rdkit import Chem

from jmcg.joint import features as ft
from jmcg.seq import score as sc

LN10 = math.log(10)
COLMAP = {'dGU': 'G_U', 'dphi1': 'phi1', 'dphi2': 'phi2', 'dphi3': 'phi3', 'dpsi1': 'psi1', 'dpsi2': 'psi2', 'dpsi3': 'psi3'}


def state_scores(ens, G_U, feats, coef):
    """Within-level energy a*G_U + theta.phi (kT) per state; coef keys as in e2 ARMS columns (log10 K units)."""
    e = np.zeros(len(ens))
    for col, b in coef.items():
        v = G_U if col == 'dGU' else feats[COLMAP[col]].to_numpy()
        # log10 K = sum b * d(col); G difference in kT = ln10 * log10 K, and dGU was already divided by ln10
        e += b * (v if col == 'dGU' else v * LN10)
    return e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--state_dir', required=True)
    ap.add_argument('--coef', required=True)
    ap.add_argument('--arm', default='J1')
    ap.add_argument('--ph', type=float, required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--n_pairs', type=int, default=20)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ens = pd.read_csv(os.path.join(a.state_dir, 'ensembles.csv'))
    G = pd.read_csv(os.path.join(a.state_dir, 'energies.csv')).groupby('smiles').G.mean()
    coef = json.load(open(a.coef))['final_coefficients_all_eligible'][a.arm] if a.arm != 'S0' else {'dGU': 1.0}
    rng = np.random.default_rng(a.seed)
    rows, pairs, fallback_levels = [], [], {}
    sdf = Chem.SDWriter(os.path.join(a.out, 'sampled_pairs.sdf'))
    for pid, g in ens.groupby('parent'):
        g = g.reset_index(drop=True)
        GU = G[g.smiles].to_numpy()
        feats, X = [], {}
        for s in g.smiles:
            f = ft.ensemble_features(s, return_ensemble=True)
            f.update(ft.graph_features(s))
            X[s] = f.pop('ensemble', None)
            feats.append(f)
        F = pd.DataFrame(feats).reindex(columns=sorted(set(COLMAP.values()) - {'G_U'} | {'failures', 'n_valid'} |
                                                       {c for f in feats for c in f}))
        ok = F.get('phi1', pd.Series([np.nan] * len(F))).notna().to_numpy()
        e = state_scores(g, GU, F.fillna(0.0), coef)
        levels = sorted(set(g.charge))
        F_L = {L: -sc.lse(-GU[g.charge == L]) for L in levels}
        pL = sc.populations([F_L[L] for L in levels], levels, a.ph)
        p, coupled = np.zeros(len(g)), np.zeros(len(g), dtype=bool)
        e_fb = coef.get('dGU', 1.0) * GU   # theta = 0 fallback: the calibrated sequential score (S1 form)
        for L, pl in zip(levels, pL):
            m = (g.charge == L).to_numpy()
            # The coupled score needs phi for every state of the level; masking out states without a valid ensemble
            # would move their mass onto the scorable ones (job 45581: glycine's zwitterion has no valid MMFF ensemble,
            # and masking gave 100 % neutral glycine at pH 7.4). Such a level uses theta = 0 instead, flagged.
            use_coupled = not any(c.startswith('dphi') for c in coef) or ok[m].all()
            if not use_coupled:
                fallback_levels[f'{pid}|charge {L}'] = dict(level_mass=float(pl), n_states=int(m.sum()),
                                                            n_without_phi=int((m & ~ok).sum()))
            ee = e[m] if use_coupled else e_fb[m]
            w = np.exp(-(ee - ee.min()))
            p[m] = pl * w / w.sum()
            coupled[m] = use_coupled and any(c.startswith('dphi') for c in coef)
        for k, r in g.iterrows():
            rows.append(dict(parent=pid, state_id=r.state_id, smiles=r.smiles, charge=r.charge, G_U=GU[k],
                             within_level_energy=e[k] if coupled[k] else e_fb[k], p_joint=p[k], phi_ok=bool(ok[k]),
                             coupled_score_used=bool(coupled[k]),
                             failures=F.loc[k, 'failures'] if 'failures' in F else ''))
        for z in rng.choice(len(g), size=a.n_pairs, p=p):
            ens_z = X[g.smiles[z]]
            k = None
            if ens_z is not None:
                mol, cids, w = ens_z
                k = int(rng.choice(len(cids), p=w))
                x = Chem.Mol(mol, confId=cids[k])
                for key, val in (('parent', pid), ('state_id', g.state_id[z]), ('state_smiles', g.smiles[z]),
                                 ('p_state', float(p[z])), ('w_conf', float(w[k])), ('arm', a.arm), ('pH', a.ph)):
                    x.SetProp(key, str(val))
                sdf.write(x)
            pairs.append(dict(parent=pid, state_id=g.state_id[z], smiles=g.smiles[z], conf_index=k,
                              has_conformer=k is not None))
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'joint_state_populations.csv'), index=False)
    pd.DataFrame(pairs).to_csv(os.path.join(a.out, 'sampled_pairs.csv'), index=False)
    sdf.close()
    json.dump(dict(arm=a.arm, coefficients=coef, ph=a.ph, seed=a.seed, fallback_levels_theta0=fallback_levels,
                   n_pairs_without_conformer=int(sum(not q['has_conformer'] for q in pairs))),
              open(os.path.join(a.out, 'run_info.json'), 'w'), indent=1)
    print(pd.DataFrame(rows).groupby('parent').apply(lambda d: d.sort_values('p_joint', ascending=False).head(3)
                                                     [['smiles', 'charge', 'p_joint']]).to_string())


if __name__ == '__main__':
    main()
