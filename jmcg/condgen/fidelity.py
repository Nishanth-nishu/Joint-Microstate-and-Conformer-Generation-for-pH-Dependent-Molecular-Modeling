"""Wrapper-fidelity check: does our TD wrapper (td_backend.TDGenerator) reproduce TD's own generate_confs.py?

Same molecules (a fixed random subset of the GEOM-QM9 test split), same checkpoint (qm9_default, pinned by path and
sha256), same conformer budget (2 x n_true, TD's default), both scored by TD's own evaluate_confs.py. Arms:
  R0, R1 = TD generate_confs.py --seed 0 / 1   (reference, and its seed-to-seed noise)
  W0, W1 = TDGenerator.generate(seed 0 / 1)    (our wrapper, on the same Mol TD builds: MolFromSmiles + AddHs)
The state layer (states.build_state) is NOT used for generation, so this isolates the generator; its verdict on each
molecule is only recorded (accepted / rejection reason).

Pre-registered rule (written before the run, 2026-10-28), per-molecule paired differences over the subset:
  noise  = R1 - R0                      wrapper = W0 - R0 (and W1 - R1)
  FAITHFUL  iff for both wrapper pairs |mean dMAT-R| <= |mean dMAT-R(noise)| + 0.02 A
            and |dCOV-R@0.5| <= |dCOV-R@0.5(noise)| + 5 pp, and the wrapper has no more model failures than TD;
  else DEGRADED if the wrapper is worse on MAT-R in both pairs, else DIFFERENT.
COV-R@0.25 and the bootstrap 95% CIs are reported but not part of the verdict (QM9 at 0.5 A is near saturation).

Run: python -m jmcg.condgen.fidelity subset  --test_csv <TD data/QM9/test_smiles.csv> --n 50 --seed 0 --out subset.csv
     python -m jmcg.condgen.fidelity wrapper --model_dir <qm9_default> --subset subset.csv --seed 0 --out W0.pkl
     python -m jmcg.condgen.fidelity compare --R0 r0.pkl --R1 r1.pkl --W0 w0.pkl --W1 w1.pkl --out compare.json
(r*/w*.pkl are evaluate_confs.py --out_results dumps.)
"""
import argparse
import json
import pickle
import time

import numpy as np
import pandas as pd
from rdkit import Chem

THRS = (0.25, 0.5)


def cmd_subset(a):
    df = pd.read_csv(a.test_csv)
    idx = np.sort(np.random.default_rng(a.seed).choice(len(df), a.n, replace=False))
    df.iloc[idx].to_csv(a.out, index=False)
    print(f'subset: {a.n} of {len(df)} rows (rng seed {a.seed}), row indices {idx.tolist()}')


def cmd_wrapper(a):
    from jmcg.condgen import states as st
    from jmcg.condgen.td_backend import TDGenerator
    gen = TDGenerator(a.model_dir, steps=a.steps, device='cpu')
    print(json.dumps(gen.describe()))
    out, log = {}, []
    for raw_smi, n_true, smi in pd.read_csv(a.subset).values:
        t0 = time.time()
        rec = dict(smiles=raw_smi, corrected_smiles=smi, n_requested=2 * int(n_true))
        try:
            st.build_state(smi)
            rec['build_state'] = 'accepted'
        except st.StateError as e:
            rec['build_state'] = f'rejected: {e}'
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rec['error'] = 'MolFromSmiles failed'
        else:
            try:
                mols, info = gen.generate(Chem.AddHs(mol), 2 * int(n_true), a.seed)
                out[smi] = mols
                rec.update(n_generated=len(mols), n_rotatable_td=info['n_rotatable_td'])
            except Exception as e:  # recorded per molecule, counted as a model failure by evaluate_confs.py
                rec['error'] = f'{type(e).__name__}: {e}'
        rec['seconds'] = round(time.time() - t0, 2)
        log.append(rec)
        print(json.dumps(rec), flush=True)
    pickle.dump(out, open(a.out, 'wb'))
    json.dump(log, open(a.out + '.log.json', 'w'), indent=1)
    print(f'wrapper: generated for {len(out)} of {len(log)} molecules')


def per_mol(path):
    d = pickle.load(open(path, 'rb'))
    rows = {}
    for (smi, csmi), r in d['results'].items():
        m = r['rmsd']
        row = dict(n_true=r['n_true'], n_model=r['n_model'], MAT_R=float(np.min(m, axis=1).mean()))
        for t in THRS:
            row[f'COV_R@{t}'] = float(np.mean(np.min(m, axis=1) < t) * 100)
        rows[csmi] = row
    return pd.DataFrame.from_dict(rows, orient='index'), d['num_failures'], d['summary']


def paired(a, b, rng):
    """Mean of b - a over molecules present in both, with a percentile bootstrap 95% CI."""
    j = a.index.intersection(b.index)
    out = {'n_paired': len(j)}
    for c in ['MAT_R'] + [f'COV_R@{t}' for t in THRS]:
        d = (b.loc[j, c] - a.loc[j, c]).to_numpy()
        bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)]
        out[c] = dict(mean=float(d.mean()), ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])
    return out


def cmd_compare(a):
    arms = {k: per_mol(getattr(a, k)) for k in ('R0', 'R1', 'W0', 'W1')}
    rng = np.random.default_rng(0)
    res = dict(arms={k: dict(summary={x: y for x, y in v[2].items() if x != 'sweep'}, model_failures=v[1])
                     for k, v in arms.items()})
    res['noise_R1_minus_R0'] = noise = paired(arms['R0'][0], arms['R1'][0], rng)
    res['W0_minus_R0'] = p0 = paired(arms['R0'][0], arms['W0'][0], rng)
    res['W1_minus_R1'] = p1 = paired(arms['R1'][0], arms['W1'][0], rng)
    res['W1_minus_W0'] = paired(arms['W0'][0], arms['W1'][0], rng)
    tol_mat = abs(noise['MAT_R']['mean']) + 0.02
    tol_cov = abs(noise['COV_R@0.5']['mean']) + 5
    fails_ok = max(arms['W0'][1], arms['W1'][1]) <= max(arms['R0'][1], arms['R1'][1])
    ok = [abs(p['MAT_R']['mean']) <= tol_mat and abs(p['COV_R@0.5']['mean']) <= tol_cov for p in (p0, p1)]
    if all(ok) and fails_ok:
        verdict = 'FAITHFUL'
    elif p0['MAT_R']['mean'] > tol_mat and p1['MAT_R']['mean'] > tol_mat:
        verdict = 'DEGRADED'
    else:
        verdict = 'DIFFERENT'
    res['rule'] = dict(tol_MAT_R=tol_mat, tol_COV_R_at_0p5=tol_cov, pairs_within_tol=ok, failures_ok=fails_ok)
    res['verdict'] = verdict
    json.dump(res, open(a.out, 'w'), indent=1)
    tab = pd.concat({k: v[0] for k, v in arms.items()}, axis=1)
    tab.to_csv(a.out.replace('.json', '_per_mol.csv'))
    for k, v in arms.items():
        s = v[2]
        print(f"{k}: n_eval={s['n_evaluated']} failures={v[1]} COV-R@0.5={s['COV-R_mean']:.2f} "
              f"MAT-R={s['MAT-R_mean']:.4f} COV-P@0.5={s['COV-P_mean']:.2f} MAT-P={s['MAT-P_mean']:.4f}")
    for k in ('noise_R1_minus_R0', 'W0_minus_R0', 'W1_minus_R1', 'W1_minus_W0'):
        p = res[k]
        print(k, f"n={p['n_paired']}", ' '.join(f"d{c}={p[c]['mean']:+.4f} [{p[c]['ci95'][0]:+.4f},{p[c]['ci95'][1]:+.4f}]"
                                               for c in ['MAT_R'] + [f'COV_R@{t}' for t in THRS]))
    print('RULE', json.dumps(res['rule']))
    print('VERDICT', verdict)


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('subset')
    s.add_argument('--test_csv', required=True); s.add_argument('--n', type=int, default=50)
    s.add_argument('--seed', type=int, default=0); s.add_argument('--out', required=True)
    s = sp.add_parser('wrapper')
    s.add_argument('--model_dir', required=True); s.add_argument('--subset', required=True)
    s.add_argument('--seed', type=int, default=0); s.add_argument('--steps', type=int, default=20)
    s.add_argument('--out', required=True)
    s = sp.add_parser('compare')
    for k in ('R0', 'R1', 'W0', 'W1'):
        s.add_argument(f'--{k}', required=True)
    s.add_argument('--out', required=True)
    a = ap.parse_args()
    dict(subset=cmd_subset, wrapper=cmd_wrapper, compare=cmd_compare)[a.cmd](a)


if __name__ == '__main__':
    main()
