"""Executable audit of the canonical data layer (counts, missingness, duplicates/conflicts, chemistry, splits, leakage).

Run: python -m jmcg.data.audit --data data_manifests/canonical_v1 --out data_manifests/canonical_v1/audit
Writes audit.json (machine-readable) and audit.md (short report). Read-only on the data.
Leakage checks:
  L1 family_key in more than one of train/val/test (must be 0 by construction)
  L2 split_group (scaffold) in more than one of train/val/test (must be 0)
  L3 near duplicates: val/test parents with a train parent at Morgan-2 (2048) Tanimoto >= 0.9 (reported, not removed)
  L4 val/test families present in the released Uni-pKa fine-tuning set (affects only evaluations of released Uni-pKa)
  L5 legacy sensitive_set_v2 split: families whose records were spread over several legacy splits
  L6 tautomer pairs whose two tautomers standardise to different families (canonical-tautomer heuristic failures)
"""
import argparse
import json
import os
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

EVAL = ['train', 'val', 'test']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--tanimoto', type=float, default=0.9)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    obs = pd.read_csv(os.path.join(a.data, 'observations.csv'), keep_default_na=False)
    states = pd.read_csv(os.path.join(a.data, 'states.csv'), keep_default_na=False)
    parents = pd.read_csv(os.path.join(a.data, 'parents.csv'), keep_default_na=False)
    dropped = pd.read_csv(os.path.join(a.data, 'dropped.csv'), keep_default_na=False)
    obs['q'] = obs['quality'].map(json.loads)
    R = {}

    R['counts'] = dict(
        parents=len(parents), states=len(states), observations=len(obs), dropped_records=len(dropped),
        dropped_by_source_reason=Counter(f"{s}: {r.split(':')[0]}" for s, r in zip(dropped.source, dropped.reason)),
        obs_by_task_role_label={f'{t}|{r}|{l}': int(n) for (t, r, l), n in
                                obs.groupby(['task', 'role', 'label_status']).size().items()},
        parents_by_task={t: int(g.parent_id.nunique()) for t, g in obs.groupby('task')},
        tautomer_pairs=int((obs.task == 'tautomer_logk').sum()),
        tautomer_pairs_measured=int(((obs.task == 'tautomer_logk') & (obs.label_status == 'measured')).sum()),
        states_by_origin=Counter(o for x in states.origins for o in x.split(';')),
        conformer_references='none in this layer: no dataset here provides 3D reference ensembles for these states '
                             '(GEOM-QM9 test set is neutral-only and lives in the TD workdir)')

    R['chemistry'] = dict(
        state_formal_charge=Counter(int(x) for x in states.formal_charge),
        state_elements=Counter(e for x in states.elements for e in re.findall('[A-Z][a-z]?', x)),
        state_n_heavy_quantiles={q: float(states.n_heavy.quantile(q)) for q in (0, .25, .5, .75, .95, 1)},
        states_with_unspecified_stereo=int((states.n_unspecified_stereo > 0).sum()),
        td_qm9_vocabulary_states=int(sum(set(re.findall('[A-Z][a-z]?', el)) <= set('HCNOF') for el in states.elements)),
        td_qm9_vocabulary_and_size_le9_states=int(sum(set(re.findall('[A-Z][a-z]?', el)) <= set('HCNOF') and nh <= 9
                                                      for el, nh in zip(states.elements, states.n_heavy))),
        parents_by_method=Counter(parents.parent_method))

    miss = {}
    for c in ['value', 'temperature', 'solvent', 'states_A', 'states_B', 'reported_state', 'scaffold']:
        col = obs[c].astype(str).str.strip()
        bad = col.isin(['', 'not recorded', 'not stated', 'nan', 'None'])
        miss[c] = {t: f'{int(bad[obs.task == t].sum())}/{int((obs.task == t).sum())}' for t in obs.task.unique()}
    R['missingness'] = miss

    # duplicates / conflicts
    pk = obs[(obs.task == 'macro_pka') & (obs.role == 'curated_eval')].copy()
    pk['ptype'] = pk.q.map(lambda q: q.get('pka_type'))
    g = pk.groupby(['family_key', 'ptype']).value
    dup = g.agg(['count', 'min', 'max'])
    dup = dup[dup['count'] > 1]
    ta = obs[obs.task == 'tautomer_logk'].copy()
    ta['pair'] = [tuple(sorted((x, y))) for x, y in zip(ta.states_A, ta.states_B)]
    ta['signed'] = [v if x <= y else -v for v, x, y in zip(ta.value, ta.states_A, ta.states_B)]
    tg = ta.groupby('pair').signed.agg(['count', 'min', 'max'])
    tg = tg[tg['count'] > 1]
    R['duplicates_conflicts'] = dict(
        iupac_same_family_same_index=int(len(dup)), iupac_same_family_same_index_spread_gt1=int(((dup['max'] - dup['min']) > 1).sum()),
        iupac_records_with_internal_spread_gt1=int(sum(q.get('spread', 0) > 1 for q in pk.q)),
        tautomer_pair_records_duplicated=int(len(tg)), tautomer_pairs_conflicting_gt1=int(((tg['max'] - tg['min']) > 1).sum()),
        tautomer_pairs_sign_conflict=int(((tg['max'] > 0) & (tg['min'] < 0)).sum()))

    cur = obs[obs.role == 'curated_eval']
    R['splits'] = dict(
        obs_by_task_split={f'{t}|{s}': int(n) for (t, s), n in cur.groupby(['task', 'split']).size().items()},
        parents_by_split={s: int(g.parent_id.nunique()) for s, g in cur.groupby('split')},
        split_groups_by_split={s: int(g.split_group.nunique()) for s, g in cur.groupby('split')},
        measured_taut_by_split={s: int(n) for s, n in cur[(cur.task == 'tautomer_logk') & (cur.label_status == 'measured')]
                                .groupby('split').size().items()})

    ev = cur[cur.split.isin(EVAL)]
    L = {}
    L['L1_families_in_several_splits'] = int((ev.groupby('family_key').split.nunique() > 1).sum())
    L['L2_split_groups_in_several_splits'] = int((ev.groupby('split_group').split.nunique() > 1).sum())
    fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    psmi = parents.set_index('parent_id').parent_smiles
    fam_split = ev.groupby('parent_id').split.first()
    fps = {p: fpg.GetFingerprint(Chem.MolFromSmiles(psmi[p])) for p in fam_split.index}
    train_ids = [p for p, s in fam_split.items() if s == 'train']
    train_fps = [fps[p] for p in train_ids]
    near = defaultdict(list)
    for p, s in fam_split.items():
        if s == 'train':
            continue
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fps[p], train_fps))
        if sims.max() >= a.tanimoto:
            j = int(sims.argmax())
            near[s].append(dict(parent=psmi[p], nearest_train=psmi[train_ids[j]], tanimoto=round(float(sims[j]), 3)))
    L['L3_near_duplicates_vs_train'] = {s: dict(n_parents=len(v), of=int((fam_split == s).sum()), examples=v[:10])
                                       for s, v in near.items()}
    upk = set(obs[obs.role == 'unipka_train'].family_key)
    L['L4_eval_families_in_unipka_train'] = {s: f"{int(g.family_key.isin(upk).groupby(g.parent_id).first().sum())}/"
                                                f"{g.parent_id.nunique()}" for s, g in ev.groupby('split')}
    leg = cur.copy()
    leg['legacy'] = leg.q.map(lambda q: q.get('legacy_split'))
    L['L5_legacy_v2_families_in_several_legacy_splits'] = int(
        (leg[leg.legacy.isin(EVAL)].groupby('family_key').legacy.nunique() > 1).sum())
    L['L5_legacy_vs_new_split_agreement'] = f"{int((leg.legacy == leg.split).sum())}/{int(leg.legacy.isin(EVAL).sum())}"
    L['L6_tautomer_pairs_two_families'] = int(sum(not q.get('same_family', True) for q in ta.q))
    lp = cur[cur.split.isin(EVAL)].copy()
    fam_of = parents.set_index('parent_id').family_key
    rows = [(f, s) for lk, s in zip(lp.linked_parents, lp.split) for f in {fam_of[x] for x in lk.split(';') if x}]
    L['L1b_linked_families_in_several_splits'] = int(pd.DataFrame(rows, columns=['f', 's']).groupby('f').s.nunique().gt(1).sum())
    L['excluded_external_overlap_obs'] = int((cur.split == 'excluded_external_overlap').sum())
    R['leakage'] = L
    R['inputs_manifest'] = json.load(open(os.path.join(a.data, 'manifest.json')))['outputs']
    json.dump(R, open(os.path.join(a.out, 'audit.json'), 'w'), indent=1, default=str)

    with open(os.path.join(a.out, 'audit.md'), 'w') as f:
        f.write('# Canonical data layer v1: audit\n\n(generated by jmcg.data.audit; numbers are from audit.json)\n\n')
        for k, v in R.items():
            if k == 'inputs_manifest':
                continue
            f.write(f'## {k}\n\n```\n{json.dumps(v, indent=1, default=str)[:6000]}\n```\n\n')
    print(json.dumps({k: R[k] for k in ('counts', 'splits', 'leakage', 'duplicates_conflicts')}, indent=1, default=str)[:12000])


if __name__ == '__main__':
    main()
