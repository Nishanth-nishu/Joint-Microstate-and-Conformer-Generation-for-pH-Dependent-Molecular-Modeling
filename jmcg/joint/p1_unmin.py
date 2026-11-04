"""P1 (EXPLORATORY, not pre-registered; added 2026-10-29): does an unminimised conformer feature path remove the M5
feature-extraction blocker? The MMFF94s-minimised path (features.ensemble_features(path='mmff'), the pre-registered
E2 definition) is the control and is not changed; validator thresholds are not changed.

Arms (same state panel, same ETKDGv3 per-conformer seeds, n = 30, same validator):
  mmff   ETKDGv3 -> MMFF94s minimisation -> validate -> Boltzmann weights on MMFF94s energies (control)
  unmin  ETKDGv3 -> validate -> uniform weights
  td     Torsional Diffusion qm9_default (E/Z lock default, jmcg.condgen.run, TD venv) -> re-validated here -> uniform
         weights. TD keeps the ETKDG seed's bonds/angles, so its bonds_angles verdict is the seed's.

Run (jmcg venv unless noted):
  python -m jmcg.joint.p1_unmin panel   --data data_manifests/canonical_v1 --e1 <E1 full dir> --smoke <smoke states dir> --out <dir>
  python -m jmcg.joint.p1_unmin extract --panel <dir>/panel.csv --path mmff|unmin --out <csv> [--procs 16]
  (TD venv) python -m jmcg.condgen.run --model_dir <qm9_default> --states_csv <dir>/td_states_<k>.csv --n_confs 30 --out <td dir>
  python -m jmcg.joint.p1_unmin td_features --panel <dir>/panel.csv --td_dir <td dir> --out <csv>
  python -m jmcg.joint.p1_unmin compare --panel ... --mmff <csv> --unmin <csv> [--td <csv>] --data ... --energies ...
         --e2_features <E2 features_v2/features.csv> --out <dir>
"""
import argparse
import glob
import json
import math
import os
from argparse import Namespace
from multiprocessing import Pool

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

from jmcg.condgen import validate as va
from jmcg.joint import e2_tautomer as e2
from jmcg.joint import features as ft

CARBOXYLATE = Chem.MolFromSmarts('[CX3](=O)[O-]')
PHENOLATE = Chem.MolFromSmarts('[O-]c')
CATS = ['neutral', 'zwitterion', 'carboxylate', 'phenolate', 'other_anion', 'cation']


def classify(smiles):
    """Primary category + flags from formal charges. zwitterion: both signs present (any net charge)."""
    m = Chem.MolFromSmiles(smiles)
    q = [a.GetFormalCharge() for a in m.GetAtoms()]
    pos, neg = any(x > 0 for x in q), any(x < 0 for x in q)
    carb, phen = m.HasSubstructMatch(CARBOXYLATE), m.HasSubstructMatch(PHENOLATE)
    if pos and neg:
        cat = 'zwitterion'
    elif neg:
        cat = 'carboxylate' if carb else 'phenolate' if phen else 'other_anion'
    elif pos:
        cat = 'cation'
    else:
        cat = 'neutral'
    return dict(category=cat, net_charge=int(sum(q)), has_carboxylate=carb, has_phenolate=phen,
                n_heavy=m.GetNumHeavyAtoms(), elements=''.join(sorted({a.GetSymbol() for a in m.GetAtoms()})))


def cmd_panel(a):
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for s in sorted(set(e2.measured_pairs(a.data)[['t1', 't2']].to_numpy().ravel())):
        rows.append(dict(smiles=s, source='tautobase_E2'))
    sm = pd.read_csv(os.path.join(a.smoke, 'ensembles.csv'))
    rows += [dict(smiles=s, source='smoke_M3') for s in sorted(set(sm.smiles))]
    # E1 external-pKa states: only the states Uni-pKa scored (energies.csv); sampled per category (seeded)
    E = pd.read_csv(os.path.join(a.e1, 'energies.csv'))
    e1 = pd.DataFrame(dict(smiles=sorted(set(E.smiles))))
    e1 = e1.join(pd.DataFrame([classify(s) for s in e1.smiles]))
    for cat, g in e1.groupby('category'):
        take = g.sample(n=min(a.per_cat, len(g)), random_state=a.seed) if len(g) > a.per_cat else g
        rows += [dict(smiles=s, source='e1_external') for s in sorted(take.smiles)]
    P = pd.DataFrame(rows)
    P = P.groupby('smiles').source.agg(lambda x: '+'.join(sorted(set(x)))).reset_index()
    P = P.join(pd.DataFrame([classify(s) for s in P.smiles]))
    P['state_key'] = [f'P{i:05d}' for i in range(len(P))]
    P.to_csv(os.path.join(a.out, 'panel.csv'), index=False)
    pool = pd.DataFrame(dict(category=e1.category)).category.value_counts().to_dict()
    json.dump(dict(per_cat=a.per_cat, seed=a.seed, e1_pool_by_category=pool,
                   panel_by_category=P.category.value_counts().to_dict(),
                   panel_by_source=P.source.value_counts().to_dict()), open(os.path.join(a.out, 'panel_info.json'), 'w'),
              indent=1)
    # TD input tables (parent = the state itself: every state here is already explicit), split for parallel runs
    T = pd.DataFrame(dict(id=P.state_key, parent_smiles=P.smiles, state_smiles=P.smiles, note=P.category))
    for k, chunk in enumerate(np.array_split(T, a.td_chunks)):
        chunk.to_csv(os.path.join(a.out, f'td_states_{k}.csv'), index=False)
    print(P.groupby(['source', 'category']).size().to_string())


def worst_pair(m, conf_id, bounds):
    """Atom labels of the worst 1-2 and 1-3 bound violations of one conformer (diagnostic for the MMFF failures)."""
    pos = m.GetConformer(conf_id).GetPositions()
    topo = Chem.GetDistanceMatrix(m)
    d = np.linalg.norm(pos[:, None] - pos[None], axis=-1)
    lb, ub = np.tril(bounds).T, np.triu(bounds)
    dev = np.maximum(lb - d, d - ub)
    out = {}
    for name, t in (('12', 1), ('13', 2)):
        mask = np.triu(topo == t, 1)
        if not mask.any():
            continue
        i, j = np.unravel_index(np.argmax(np.where(mask, dev, -np.inf)), dev.shape)
        ai, aj = m.GetAtomWithIdx(int(i)), m.GetAtomWithIdx(int(j))
        lab = lambda x: x.GetSymbol() + ('-' if x.GetFormalCharge() < 0 else '+' if x.GetFormalCharge() > 0 else '')
        mid = ''
        if t == 2:
            k = [n.GetIdx() for n in ai.GetNeighbors() if topo[n.GetIdx(), j] == 1][0]
            mid = m.GetAtomWithIdx(k).GetSymbol() + ('(ar)' if m.GetAtomWithIdx(k).GetIsAromatic() else '')
        out[f'worst_{name}'] = f'{lab(ai)}-{mid + "-" if mid else ""}{lab(aj)}'
        out[f'worst_{name}_dev'] = round(float(dev[i, j]), 3)
        out[f'worst_{name}_dist'] = round(float(d[i, j]), 3)
        out[f'worst_{name}_bounds'] = f'{lb[i, j]:.3f}-{ub[i, j]:.3f}'
    return out


def _extract(args):
    smiles, path = args
    try:
        r = ft.ensemble_features(smiles, path=path)
    except Exception as e:  # noqa: BLE001 - recorded as a failure
        r = dict(n_generated=0, n_valid=0, failures=[f'{type(e).__name__}: {e}'[:200]])
    r.pop('ensemble', None)
    r.update(ft.graph_features(smiles))
    if path == 'mmff' and r.get('n_generated', 0) and not r.get('n_valid'):
        # diagnostic on the first minimised conformer of a state with no valid MMFF conformer: which bound is violated
        try:
            m2 = Chem.AddHs(Chem.MolFromSmiles(smiles))
            p = AllChem.ETKDGv3(); p.randomSeed = ft.conformer_seeds(0, 1)[0]
            if AllChem.EmbedMolecule(m2, p) == 0:
                b = rdDistGeom.GetMoleculeBoundsMatrix(m2)
                mf = Chem.Mol(m2)
                ff = AllChem.MMFFGetMoleculeForceField(mf, AllChem.MMFFGetMoleculeProperties(mf, mmffVariant='MMFF94s'))
                ff.Minimize(maxIts=2000)
                r.update(worst_pair(_with_coords(m2, mf), 0, b))
                r.update({'unmin_' + k: v for k, v in worst_pair(m2, 0, b).items()})
        except Exception as e:  # noqa: BLE001
            r['diag_error'] = str(e)[:120]
    fails = pd.Series(r.get('failures', []), dtype=str)
    r['failure_counts'] = json.dumps(fails.value_counts().to_dict())
    r['failures'] = '; '.join(sorted(set(r.get('failures', []))))
    r['smiles'], r['path'] = smiles, path
    return r


def _with_coords(m, mf):
    out = Chem.Mol(m)
    c = out.GetConformer(0)
    for i, xyz in enumerate(mf.GetConformer(0).GetPositions()):
        c.SetAtomPosition(i, xyz.tolist())
    return out


def cmd_extract(a):
    P = pd.read_csv(a.panel)
    with Pool(a.procs) as p:
        rows = p.map(_extract, [(s, a.path) for s in P.smiles], chunksize=2)
    F = pd.DataFrame(rows)
    F.to_csv(a.out, index=False)
    print(a.path, 'states', len(F), 'with phi', int(F.phi1.notna().sum()) if 'phi1' in F else 0)


def cmd_td_features(a):
    P = pd.read_csv(a.panel)
    rows = []
    for _, r in P.iterrows():
        sdf = os.path.join(a.td_dir, f'{r.state_key}__td.sdf')
        rec = dict(smiles=r.smiles, path='td', **ft.graph_features(r.smiles))
        js = os.path.join(a.td_dir, f'{r.state_key}.json')
        info = json.load(open(js)) if os.path.exists(js) else {}
        rec['td_category'] = info.get('category')
        rec['td_failures'] = '; '.join(f"{x.get('method', x['stage'])}: {x['reason']}" for x in info.get('failures', []))[:300]
        if not os.path.exists(sdf):
            rec.update(n_generated=0, n_valid=0, failures='no TD output: ' + (rec['td_failures'] or 'missing'))
            rows.append(rec)
            continue
        ms = [x for x in Chem.SDMolSupplier(sdf, removeHs=False, sanitize=True)]
        rec['td_valid_in_run'] = int(sum(x is not None and x.GetProp('valid') == 'True' for x in ms))
        ms = [x for x in ms if x is not None]
        if not ms:
            rec.update(n_generated=0, n_valid=0, failures='TD SDF unreadable')
            rows.append(rec)
            continue
        m = Chem.Mol(ms[0]); m.RemoveAllConformers()
        for x in ms:
            m.AddConformer(x.GetConformer(), assignId=True)
        target = Chem.MolToSmiles(Chem.MolFromSmiles(r.smiles))
        f = ft.features_from_conformers(m, target)
        fails = pd.Series(f['failures'], dtype=str)
        f['failure_counts'] = json.dumps(fails.value_counts().to_dict())
        f['failures'] = '; '.join(sorted(set(f['failures'])))
        rec.update(f)
        rows.append(rec)
    F = pd.DataFrame(rows)
    F.to_csv(a.out, index=False)
    print('td states', len(F), 'with phi', int(F.phi1.notna().sum()) if 'phi1' in F else 0)


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [round(c - h, 4), round(c + h, 4)]


def _oof(d, cols, y, cv):
    p = np.zeros(len(d))
    for f in range(5):
        p[cv == f] = e2.fit_predict(d[cols].to_numpy(), y, cv != f, cv == f)[0]
    return p


def cmd_compare(a):
    os.makedirs(a.out, exist_ok=True)
    P = pd.read_csv(a.panel).set_index('smiles')
    arms = {'mmff': a.mmff, 'unmin': a.unmin}
    if a.td:
        arms['td'] = a.td
    F = {}
    for k, v in arms.items():
        G = pd.read_csv(v).set_index('smiles')
        if 'psi1' not in G:  # files written before graph features were added to extract/td_features (job 46561)
            G = G.join(pd.DataFrame([ft.graph_features(s) for s in G.index], index=G.index))
            arms[k] = os.path.join(a.out, f'features_{k}_with_psi.csv')
            G.reset_index().to_csv(arms[k], index=False)
        F[k] = G
    R = dict(note='EXPLORATORY (not pre-registered). Validator thresholds unchanged. mmff = control (E2 definition).')
    # 1. control reproduction: mmff path here vs the E2 features_v2 file on the shared states
    V2 = pd.read_csv(a.e2_features).set_index('smiles')
    common = V2.index.intersection(F['mmff'].index)
    x, y = V2.loc[common], F['mmff'].loc[common]
    R['control_reproduces_E2_features_v2'] = dict(
        n_states=len(common), n_valid_equal=int((x.n_valid == y.n_valid).sum()),
        phi_present_equal=int((x.phi1.notna() == y.phi1.notna()).sum()),
        max_abs_diff={k: float(np.nanmax(np.abs(x[k] - y[k]))) for k in ('phi1', 'phi2', 'phi3')})
    # 2. extraction success per category and arm
    rows = []
    for arm, G in F.items():
        G = G.join(P[['category', 'source', 'has_carboxylate', 'has_phenolate', 'state_key']], how='inner')
        G['ok'] = G.phi1.notna() if 'phi1' in G else False
        for scope, H in (('all', G), ('e1_external', G[G.source.str.contains('e1_external')]),
                         ('tautobase_E2', G[G.source.str.contains('tautobase')])):
            for cat, g in list(H.groupby('category')) + [('zwitterion_with_carboxylate',
                                                          H[(H.category == 'zwitterion') & H.has_carboxylate]),
                                                         ('any_phenolate', H[H.has_phenolate])]:
                if not len(g):
                    continue
                fc = {}
                for s in g.failure_counts.dropna():
                    for kk, vv in json.loads(s).items():
                        key = kk.split(': ', 1)[-1] if 'invalid' in kk else kk[:60]
                        fc[key] = fc.get(key, 0) + vv
                rows.append(dict(arm=arm, scope=scope, category=cat, n_states=len(g), n_with_phi=int(g.ok.sum()),
                                 frac_with_phi=round(float(g.ok.mean()), 4), wilson95=wilson(int(g.ok.sum()), len(g)),
                                 conformers_valid=int(g.n_valid.sum()), conformers_generated=int(g.n_generated.sum()),
                                 frac_conformers_valid=round(float(g.n_valid.sum() / max(1, g.n_generated.sum())), 4),
                                 conformer_failure_reasons=json.dumps(dict(sorted(fc.items(), key=lambda t: -t[1])[:6]))))
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(a.out, 'extraction_rates.csv'), index=False)
    R['extraction_rates_all'] = S[S.scope == 'all'][['arm', 'category', 'n_states', 'n_with_phi', 'frac_with_phi',
                                                     'wilson95', 'frac_conformers_valid']].to_dict('records')
    # paired success: states without phi under mmff that gain phi under unmin (and the reverse)
    ok = {k: G.phi1.notna() if 'phi1' in G else pd.Series(False, index=G.index) for k, G in F.items()}
    idx = ok['mmff'].index.intersection(ok['unmin'].index)
    tr = pd.DataFrame(dict(mmff=ok['mmff'][idx], unmin=ok['unmin'][idx], category=P.category[idx]))
    R['paired_success_mmff_vs_unmin'] = {cat: dict(both=int((g.mmff & g.unmin).sum()), mmff_only=int((g.mmff & ~g.unmin).sum()),
                                                   unmin_only=int((~g.mmff & g.unmin).sum()), neither=int((~g.mmff & ~g.unmin).sum()))
                                         for cat, g in tr.groupby('category')}
    # MMFF failure diagnostics (worst bound violation after minimisation, same conformer before minimisation)
    M = F['mmff'].join(P[['category']], how='inner')
    if 'worst_13' in M:
        D = M[M.worst_13.notna()]
        R['mmff_failure_diagnostics'] = {cat: dict(n=len(g), worst_13_top=g.worst_13.value_counts().head(5).to_dict(),
                                                   worst_12_top=g.worst_12.value_counts().head(3).to_dict(),
                                                   median_worst_13_dev_after_mmff=float(g.worst_13_dev.median()),
                                                   median_worst_13_dev_same_conformer_unminimised=float(g.unmin_worst_13_dev.median()),
                                                   median_worst_12_dev_after_mmff=float(g.worst_12_dev.median()))
                                         for cat, g in D.groupby('category')}
        D.reset_index()[['smiles', 'category', 'worst_12', 'worst_12_dev', 'worst_12_dist', 'worst_12_bounds', 'worst_13',
                         'worst_13_dev', 'worst_13_dist', 'worst_13_bounds', 'unmin_worst_13_dev', 'unmin_worst_12_dev']
                        ].to_csv(os.path.join(a.out, 'mmff_failure_diagnostics.csv'), index=False)
    # 3. feature distributions on matched states (phi under both arms)
    dist = {}
    for other in [k for k in F if k != 'mmff']:
        mi = F['mmff'].index[F['mmff'].phi1.notna()].intersection(F[other].index[F[other].phi1.notna()])
        dd = {}
        for scope, sel in (('all', mi), ('neutral', mi[P.category[mi] == 'neutral']),
                           ('charged', mi[P.category[mi] != 'neutral'])):
            if len(sel) < 3:
                dd[scope] = dict(n=len(sel))
                continue
            st = dict(n=len(sel))
            for k in ('phi2', 'phi3'):
                u, v = F['mmff'].loc[sel, k], F[other].loc[sel, k]
                st[k] = dict(mean_mmff=float(u.mean()), mean_other=float(v.mean()), sd_mmff=float(u.std()), sd_other=float(v.std()),
                             mean_diff_other_minus_mmff=float((v - u).mean()), median_abs_diff=float((v - u).abs().median()),
                             pearson=float(np.corrcoef(u, v)[0, 1]), spearman=float(u.rank().corr(v.rank())))
            st['phi1_note'] = 'uniform weights on this arm: phi1 = -ln n_valid (no conformational information)'
            dd[scope] = st
        dist[f'{other}_vs_mmff'] = dd
    R['feature_distributions_matched'] = dist
    # 4. E2 predictions (exploratory): same pre-registered fitting procedure, each arm's features
    t = e2.measured_pairs(a.data)
    E = pd.read_csv(a.energies).groupby('smiles').G.mean()
    designs = {}
    for arm, path in arms.items():
        _, d = e2.design(Namespace(data=a.data, features=path, energies=a.energies))
        designs[arm] = d.set_index('obs_id')
    pred = {}
    ids_all = {k: set(v.index) for k, v in designs.items()}
    R['e2_eligible_pairs'] = {k: len(v) for k, v in ids_all.items()} | dict(n_measured=len(t))

    def run_set(ids, arm_list, tag):
        ids = sorted(ids)
        base = designs['mmff'] if 'mmff' in arm_list and set(ids) <= ids_all['mmff'] else designs[arm_list[0]]
        d0 = base.loc[ids]
        cv = e2.scaffold_folds(d0.group.tolist(), 5)
        y = d0.value.to_numpy()
        err, out = {'S0': d0.dGU.to_numpy() - y}, {}
        err['S1'] = _oof(d0, ['dGU'], y, cv) - y
        err['C1'] = _oof(d0, e2.ARMS['C1'], y, cv) - y
        preds = {}
        for arm in arm_list:
            dd = designs[arm].loc[ids]
            preds[arm] = _oof(dd, e2.ARMS['J1'], y, cv)
            err[f'J1_{arm}'] = preds[arm] - y
            for k in ('J1-phi3',):
                err[f'{k}_{arm}'] = _oof(dd, e2.ARMS[k], y, cv) - y
            err[f'S1+phi3_{arm}'] = _oof(dd, ['dGU', 'dphi3'], y, cv) - y
        g = d0.group
        out = dict(n_pairs=len(ids), n_groups=int(g.nunique()), predict_zero_MAE=float(np.abs(y).mean()),
                   MAE={k: float(np.abs(v).mean()) for k, v in err.items()},
                   delta_S={arm: e2.boot_delta(err['S1'], err[f'J1_{arm}'], g) for arm in arm_list},
                   delta_C={arm: e2.boot_delta(err['C1'], err[f'J1_{arm}'], g) for arm in arm_list})
        if len(arm_list) > 1:
            out['J1_mmff_minus_J1_other'] = {arm: e2.boot_delta(err['J1_mmff'], err[f'J1_{arm}'], g)
                                             for arm in arm_list if arm != 'mmff'}
            out['prediction_agreement'] = {arm: dict(pearson=float(np.corrcoef(preds['mmff'], preds[arm])[0, 1]),
                                                     mean_abs_diff=float(np.abs(preds['mmff'] - preds[arm]).mean()))
                                           for arm in arm_list if arm != 'mmff'}
        pd.DataFrame(dict(obs_id=ids, y=y, **{f'err_{k}': v for k, v in err.items()})).to_csv(
            os.path.join(a.out, f'oof_{tag}.csv'), index=False)
        return out

    matched = set.intersection(*ids_all.values())
    R['predictions_matched_pairs'] = run_set(matched, list(arms), 'matched')
    for arm in arms:
        if arm != 'mmff':
            R[f'predictions_all_eligible_{arm}'] = run_set(ids_all[arm], [arm], f'eligible_{arm}')
    # exploratory final coefficients per arm (for the sampler only; not evaluated)
    final = {}
    for arm, d in designs.items():
        y = d.value.to_numpy()
        for name, cols in (('J1', e2.ARMS['J1']), ('S1', ['dGU'])):
            X = d[cols].to_numpy()
            rms = np.sqrt((X ** 2).mean(0)); rms[rms == 0] = 1.0
            beta = np.linalg.solve((X / rms).T @ (X / rms) + e2.LAM * np.eye(len(cols)), (X / rms).T @ y) / rms
            final.setdefault(arm, {})[name] = dict(zip(cols, beta.tolist()))
    for arm in arms:
        json.dump(dict(note=f'EXPLORATORY P1 refit of the E2 J1/S1 arms on {arm}-path features (all pairs eligible under {arm})',
                       n_pairs=len(designs[arm]), final_coefficients_all_eligible=final[arm]),
                  open(os.path.join(a.out, f'coef_{arm}.json'), 'w'), indent=1)
    R['final_coefficients_exploratory'] = final
    json.dump(R, open(os.path.join(a.out, 'p1_result.json'), 'w'), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float))


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('panel')
    for k in ('data', 'e1', 'smoke', 'out'):
        s.add_argument(f'--{k}', required=True)
    s.add_argument('--per_cat', type=int, default=150); s.add_argument('--seed', type=int, default=0)
    s.add_argument('--td_chunks', type=int, default=12)
    s = sp.add_parser('extract')
    s.add_argument('--panel', required=True); s.add_argument('--path', choices=['mmff', 'unmin'], required=True)
    s.add_argument('--out', required=True); s.add_argument('--procs', type=int, default=16)
    s = sp.add_parser('td_features')
    for k in ('panel', 'td_dir', 'out'):
        s.add_argument(f'--{k}', required=True)
    s = sp.add_parser('compare')
    for k in ('panel', 'mmff', 'unmin', 'data', 'energies', 'e2_features', 'out'):
        s.add_argument(f'--{k}', required=True)
    s.add_argument('--td', default=None)
    a = ap.parse_args()
    dict(panel=cmd_panel, extract=cmd_extract, td_features=cmd_td_features, compare=cmd_compare)[a.cmd](a)


if __name__ == '__main__':
    main()
