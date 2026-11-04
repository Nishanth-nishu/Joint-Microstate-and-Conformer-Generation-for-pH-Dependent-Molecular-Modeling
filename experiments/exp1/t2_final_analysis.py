"""EXP1-V2-T2 part D: the pre-declared final analysis of the six 20-epoch runs (DESIGN.md "Endpoints and rules").

Read only. Inputs are the outputs of ``exp1_v2_t2_e1.sbatch``: per TD sampling seed s, the evaluate_confs.py results of
TD, ETKDGv3, the T1 arms (constrained xTB, MMFF, T1 L(proj)fix / L(noproj)fix) and the six T2 read-outs
``td{s}_T2_{proj,noproj}_s{0,1,2}``; ``geom_validate_seed{s}.json`` written with every arm; the read-out logs; the
training logs. Metrics are the E1 compare's (per molecule, then the unweighted mean over the paired molecules); CIs are
the same percentile bootstrap over molecules (5,000 draws, ``default_rng(0)``). Seed-averaged = per-molecule metric
averaged over the three training seeds, then the paired bootstrap. Endpoint verdicts follow the DESIGN.md rules
mechanically; nothing here selects a seed or a setting.

Run (TD venv): python experiments/exp1/t2_final_analysis.py --n200 <n200> --t1 <T1_e1> --t2 <T2_e1> --geom_dir <dir>
                   --train_dir <T2 train dir> --out final_analysis.json
"""
from __future__ import annotations

import json
import os
from argparse import ArgumentParser

import numpy as np

from jmcg.condgen.e1_confs import per_mol

ARMS = ('proj', 'noproj')
SEEDS = (0, 1, 2)
MET = ('MAT_R', 'MAT_P', 'COV_R@0.25', 'COV_R@0.5', 'COV_P@0.25', 'COV_P@0.5')
GEO = ('bond_mae', 'angle_mae', 'ring_dih_mae', 'acyclic_dih_mae')
LOWER_BETTER = {'MAT_R', 'MAT_P', *GEO}


def boot(d, rng):
    d = np.asarray(d, float)
    d = d[~np.isnan(d)]
    bs = d[rng.integers(0, len(d), (5000, len(d)))].mean(1)
    return dict(mean=float(d.mean()), ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], n=int(len(d)))


def favourable(r, metric):
    lo, hi = r['ci95']
    return hi < 0 if metric in LOWER_BETTER else lo > 0


def excludes0(r):
    return r['ci95'][0] > 0 or r['ci95'][1] < 0


NEAR_DUP = 'N#CC12C3C4CN1C2C43'  # the one E1 graph found in the training cache with stereo removed (t1_leak_check.py)


def flat(smi):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)


def strata(common):
    """Addendum C (exploratory, fixed before any T2 number): molecule groups from the graph only."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    out = {}
    for cs in common:
        m = Chem.MolFromSmiles(cs)
        nrot = rdMolDescriptors.CalcNumRotatableBonds(m)
        ri = m.GetRingInfo()
        small = any(len(r) <= 4 for r in ri.AtomRings())
        out[cs] = {'rotatable': '0' if nrot == 0 else '1' if nrot == 1 else '>=2',
                   'rings': str(min(ri.NumRings(), 2)).replace('2', '>=2'),
                   'ring_3or4': 'yes' if small else 'no'}
    return out


def e1(a, s):
    p = {'td': f'{a.n200}/td{s}.eval.pkl', 'etkdg': f'{a.n200}/etkdg{s}.eval.pkl',
         **{k: f'{a.t1}/td{s}_{k}.eval.pkl' for k in ('xtb', 'mmff', 'Lprojfix', 'Lnoprojfix')},
         **{f'{m}_s{k}': f'{a.t2}/td{s}_T2_{m}_s{k}.eval.pkl' for m in ARMS for k in SEEDS}}
    dfs = {k: per_mol(v) for k, v in p.items()}
    common = sorted(set.intersection(*[set(df.index) for df, _ in dfs.values()]))
    t = {k: df.loc[common] for k, (df, _) in dfs.items()}
    for m in ARMS:
        t[f'{m}_avg'] = sum(t[f'{m}_s{k}'] for k in SEEDS) / len(SEEDS)
    return t, common, {k: f for k, (_, f) in dfs.items()}


def readout_logs(a, s):
    out = {}
    for k, f in [(f'{m}_s{j}', f'{a.t2}/td{s}_T2_{m}_s{j}.pkl.log.json') for m in ARMS for j in SEEDS] + \
            [(k, f'{a.t1}/td{s}_{k}.pkl.log.json') for k in ('xtb', 'mmff', 'Lprojfix', 'Lnoprojfix')]:
        if os.path.exists(f):
            out[k] = json.load(open(f))['summary']
    return out


def training(a):
    out = {}
    for m in ARMS:
        for k in SEEDS:
            d = f'{a.train_dir}/T2_{m}_s{k}'
            ep = [json.loads(x) for x in open(f'{d}/epochs.jsonl')]
            it = list({(r['epoch'], r['it']): r for r in map(json.loads, open(f'{d}/train_log.jsonl'))}.values())
            best = min(ep, key=lambda e: e['val_objective'])
            g = {}
            for key in ('grad_l_head', 'grad_conv_layers', 'grad_embeddings'):
                by_ep = {}
                for r in it:
                    by_ep.setdefault(r['epoch'], []).append(r.get(key, np.nan))
                g[key] = {str(e): dict(median=float(np.nanmedian(v)), p95=float(np.nanpercentile(v, 95)),
                                       max=float(np.nanmax(v))) for e, v in sorted(by_ep.items())}
                g[key]['nonfinite'] = int(sum(not np.isfinite(r.get(key, 0.0)) for r in it))
            out[f'{m}_s{k}'] = dict(
                epochs=len(ep), last_epoch=ep[-1]['epoch'], best_epoch=best['epoch'],
                best_val_objective=best['val_objective'], best_val_geo=best['val_geo'],
                best_val_geo_base=best['val_geo_base'], best_val_td=best['val_td'],
                train_hours=float(sum(e['train_seconds'] for e in ep) / 3600),
                curve=[{k2: e[k2] for k2 in ('epoch', 'train_loss', 'train_geo', 'train_geo_base', 'train_td',
                                              'val_objective', 'val_geo', 'val_geo_base', 'val_td', 'train_seconds')}
                       for e in ep],
                grad=g, nonfinite_loss=int(sum(not np.isfinite(r['loss']) for r in it)))
    return out


def main():
    ap = ArgumentParser()
    for x in ('--n200', '--t1', '--t2', '--geom_dir', '--train_dir', '--out'):
        ap.add_argument(x, required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(0)
    rng_c = np.random.default_rng(1)  # addendum C draws kept apart, so the endpoint CIs do not depend on them
    res = dict(training=training(a), e1={}, geometry={}, readout={}, endpoints={})
    for s in (0, 1):
        t, common, fails = e1(a, s)
        r = dict(n_attempted=200, n_common=len(common), failures=fails,
                 means={k: {c: float(v[c].mean()) for c in MET} for k, v in t.items()},
                 between_seed_sd={m: {c: float(np.std([t[f'{m}_s{k}'][c].mean() for k in SEEDS], ddof=1)) for c in MET}
                                  for m in ARMS},
                 paired={})
        pairs = [(f'{m}_s{k}', 'td') for m in ARMS for k in SEEDS] + \
                [(f'{m}_avg', ref) for m in ARMS for ref in ('td', 'xtb', 'mmff', 'etkdg')] + \
                [('proj_avg', 'noproj_avg'), ('Lprojfix', 'td'), ('xtb', 'td'), ('proj_s0', 'Lprojfix'),
                 ('noproj_s0', 'Lnoprojfix')] + [(f'proj_s{k}', f'noproj_s{k}') for k in SEEDS] + \
                [(f'proj_s{k}', 'xtb') for k in SEEDS]
        for x, y in pairs:
            r['paired'][f'{x}-{y}'] = {c: boot(t[x][c] - t[y][c], rng) for c in MET}
        st = strata(common)
        r['strata'] = {}
        for key in ('rotatable', 'rings', 'ring_3or4'):
            for g in sorted({v[key] for v in st.values()}):
                ids = [cs for cs in common if st[cs][key] == g]
                r['strata'][f'{key}={g}'] = dict(n=len(ids), **{
                    f'{x}-{y}': boot((t[x].loc[ids, 'MAT_R'] - t[y].loc[ids, 'MAT_R']).to_numpy(), rng_c)
                    for x, y in (('proj_avg', 'td'), ('proj_avg', 'xtb'), ('proj_avg', 'noproj_avg'))})
        keep = [cs for cs in common if flat(cs) != flat(NEAR_DUP)]
        r['without_near_duplicate'] = dict(n=len(keep), removed=len(common) - len(keep), **{
            f'{x}-{y}': boot((t[x].loc[keep, 'MAT_R'] - t[y].loc[keep, 'MAT_R']).to_numpy(), rng_c)
            for x, y in [(f'proj_s{k}', 'td') for k in SEEDS] + [('proj_avg', 'noproj_avg'), ('proj_avg', 'xtb')]})
        r['variability_MAT_R'] = {m: dict(per_training_seed=[float(t[f'{m}_s{k}']['MAT_R'].mean()) for k in SEEDS],
                                          between_training_seed_sd=r['between_seed_sd'][m]['MAT_R'],
                                          bootstrap_ci_width_minus_td=float(np.diff(
                                              r['paired'][f'{m}_avg-td']['MAT_R']['ci95'])[0]))
                                  for m in ARMS}
        res['e1'][f'seed{s}'] = r
        g = json.load(open(f'{a.geom_dir}/geom_validate_seed{s}.json'))
        per, summ = g['per_molecule'], g['summary']
        ok = sorted(set.intersection(*[set(v) for v in per.values()]))
        for m in ARMS:
            per[f'{m}_avg'] = {cs: {q: float(np.nanmean([per[f'{m}_s{k}'][cs][q] for k in SEEDS])) for q in GEO}
                               for cs in ok}
        gp = {}
        for x, y in [(f'{m}_avg', ref) for m in ARMS for ref in ('td', 'xtb')] + [('proj_avg', 'noproj_avg')] + \
                [(f'{m}_s{k}', 'td') for m in ARMS for k in SEEDS]:
            gp[f'{x}-{y}'] = {q: boot([per[x][cs][q] - per[y][cs][q] for cs in ok], rng) for q in GEO}
        res['geometry'][f'seed{s}'] = dict(molecules=len(ok), conformer_checks=summ['conformer_checks'],
                                           means={k: {q: float(np.nanmean([per[k][cs][q] for cs in ok])) for q in GEO}
                                                  for k in per},
                                           molecule_errors=summ['molecule_errors'], paired=gp)
        res['readout'][f'seed{s}'] = readout_logs(a, s)

    P0, P1 = res['e1']['seed0']['paired'], res['e1']['seed1']['paired']
    ep = res['endpoints']
    per_seed = {k: P0[f'proj_s{k}-td']['MAT_R'] for k in SEEDS}
    ep['1_improves_on_td'] = dict(
        per_training_seed=per_seed, td_seed1_sign={k: P1[f'proj_s{k}-td']['MAT_R']['mean'] for k in SEEDS},
        verdict='yes' if all(favourable(v, 'MAT_R') for v in per_seed.values())
        and all(P1[f'proj_s{k}-td']['MAT_R']['mean'] < 0 for k in SEEDS) else 'no')
    ep['2_persists_across_seeds'] = dict(
        between_seed_sd_MAT_R={m: res['e1']['seed0']['between_seed_sd'][m]['MAT_R'] for m in ARMS},
        seed_averaged_proj_minus_td=P0['proj_avg-td']['MAT_R'],
        seed_averaged_proj_minus_td_tdseed1=P1['proj_avg-td']['MAT_R'])
    sign = sum(P0[f'proj_s{k}-noproj_s{k}']['MAT_R']['mean'] < 0 for k in SEEDS)
    pr = P0['proj_avg-noproj_avg']['MAT_R']
    ep['3_projector'] = dict(seed_averaged=pr, seeds_favouring_proj=int(sign),
                             tdseed1=P1['proj_avg-noproj_avg']['MAT_R'],
                             verdict='improves' if favourable(pr, 'MAT_R') and sign >= 2 else 'not shown')
    vx = P0['proj_avg-xtb']['MAT_R']
    ep['4_vs_xtb'] = dict(seed_averaged=vx, tdseed1=P1['proj_avg-xtb']['MAT_R'],
                          verdict='beats' if favourable(vx, 'MAT_R') else 'indistinguishable' if not excludes0(vx)
                          else 'worse')
    lt = P0['proj_s0-Lprojfix']['MAT_R']
    ep['5_longer_training'] = dict(T2s0_minus_T1s0=lt, tdseed1=P1['proj_s0-Lprojfix']['MAT_R'],
                                   verdict='helped' if favourable(lt, 'MAT_R') else 'no change' if not excludes0(lt)
                                   else 'hurt')
    ep['6_where'] = dict(geometry_seed0={k: res['geometry']['seed0']['paired'][k]
                                         for k in ('proj_avg-td', 'proj_avg-xtb', 'proj_avg-noproj_avg')},
                         cov_proj_avg_minus_td={c: P0['proj_avg-td'][c] for c in MET if c.startswith('COV')})
    ep['7_cost'] = {s: {k: v.get('cpu_seconds_per_conformer') for k, v in res['readout'][s].items()}
                    for s in res['readout']}
    json.dump(res, open(a.out, 'w'), indent=1, default=float)

    def f(r):
        return f"{r['mean']:+.4f} [{r['ci95'][0]:+.4f}, {r['ci95'][1]:+.4f}]"
    for s in (0, 1):
        e = res['e1'][f'seed{s}']
        print(f"\n== TD seed {s}: {e['n_common']} paired molecules")
        for k, v in e['means'].items():
            print(f"  {k:12s} " + ' '.join(f"{c}={v[c]:.3f}" for c in MET))
        for k in ('proj_avg-td', 'noproj_avg-td', 'proj_avg-noproj_avg', 'proj_avg-xtb', 'noproj_avg-xtb',
                  'proj_s0-Lprojfix', *[f'proj_s{j}-td' for j in SEEDS], *[f'proj_s{j}-noproj_s{j}' for j in SEEDS]):
            print(f"  {k:22s} MAT_R {f(e['paired'][k]['MAT_R'])}  COV_R@0.25 {f(e['paired'][k]['COV_R@0.25'])}")
        g = res['geometry'][f'seed{s}']
        for k, v in g['means'].items():
            c = g['conformer_checks'].get(k, {})
            print(f"  geo {k:12s} " + ' '.join(f"{q}={v[q]:.4f}" for q in GEO) +
                  f"  idfail={c.get('identity_fail')} stereoΔ={c.get('stereo_changed_vs_td_input')}"
                  f" clash={c.get('clash_conformers')}")
    for k, v in ep.items():
        print(k, v.get('verdict', ''), json.dumps({x: y for x, y in v.items() if x != 'verdict'}, default=float)[:400])


if __name__ == '__main__':
    main()
