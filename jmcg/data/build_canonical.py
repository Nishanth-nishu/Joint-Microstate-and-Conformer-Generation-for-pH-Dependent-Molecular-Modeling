"""Build the canonical data layer v1: parents, microstates and observations with provenance and splits.

Run (TD-free jmcg venv): python -m jmcg.data.build_canonical --ext /scratch/nishanth.r/jmcg/ext \
        --pka_curated data_manifests/sensitive_set_v2/pka_curated.csv --out data_manifests/canonical_v1

Sources (experimental labels; nothing is imputed):
  macro_pka       IUPAC v2.4a as curated in sensitive_set_v2 (one row per molecule x pKa index, median of repeats).
                  Only the reported microstate is known; the acid/base microstate sets are NOT (states_A/B empty).
  macro_pka       Uni-pKa external test sets (SAMPL6/7/8, Novartis acid/base): explicit acid/base microstate sets.
  macro_pka       Dwar-iBonD = the released Uni-pKa model's fine-tuning set (role unipka_train; for leakage flags and
                  for fitting only).
  tautomer_logk   Tautobase, aqueous, numeric log K; label_status measured / estimated / calculated from its tag.
Roles and splits:
  curated_eval  -> train / val / test by canonical.assign_split(split_group) (scaffold, or family if acyclic)
  external_test -> 'external_test' (never used for fitting or selection)
  unipka_train  -> 'unipka_train'
  A curated family that also occurs in an external test set -> 'excluded_external_overlap'.
Every dropped record is written to dropped.csv with its reason.
"""
import argparse
import csv
import hashlib
import json
import os
import platform
from collections import Counter, defaultdict

import rdkit
from rdkit import Chem

from jmcg.data import canonical as cn
from jmcg.data.build_sensitive_set import MEASURED, build_tautobase

UNIPKA_SETS = {'sampl6': 'external_test', 'sampl7': 'external_test', 'sampl8': 'external_test',
               'novartis_acid': 'external_test', 'novartis_base': 'external_test', 'dwar-iBond': 'unipka_train'}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


class Builder:
    def __init__(self):
        self.states, self.parents, self.obs, self.dropped = {}, {}, [], []
        self.state_parent_conflicts = []

    def add_state(self, smiles, origin):
        can, mol = cn.microstate(smiles)
        sid = cn.state_id(can)
        if sid not in self.states:
            unspec = sum(1 for s in Chem.FindPotentialStereo(mol) if s.specified == Chem.StereoSpecified.Unspecified)
            self.states[sid] = dict(state_id=sid, smiles=can, formal_charge=Chem.GetFormalCharge(mol),
                                    n_heavy=mol.GetNumHeavyAtoms(), n_charged_atoms=sum(1 for a in mol.GetAtoms()
                                                                                        if a.GetFormalCharge()),
                                    elements=''.join(sorted({a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms()})),
                                    n_unspecified_stereo=unspec, origins=set(), parent_id=None, _mol=mol)
        self.states[sid]['origins'].add(origin)
        return sid

    def add_parent(self, sid):
        st = self.states[sid]
        ident = cn.parent_identity(st['_mol'])
        pid = ident['parent_id']
        if st['parent_id'] and st['parent_id'] != pid:
            self.state_parent_conflicts.append((sid, st['parent_id'], pid))
        st['parent_id'] = st['parent_id'] or pid
        if pid not in self.parents:
            pm = Chem.MolFromSmiles(ident['parent_smiles'])
            self.parents[pid] = dict(**ident, split_group=cn.split_group(ident), n_heavy=pm.GetNumHeavyAtoms(),
                                     elements=''.join(sorted({a.GetSymbol() for a in Chem.AddHs(pm).GetAtoms()})))
        return pid

    def drop(self, source, record, reason):
        self.dropped.append(dict(source=source, source_record=record, reason=str(reason)[:300]))


def pka_curated(b, path):
    for r in csv.DictReader(open(path)):
        rec = f"{r['inchikey']}|{r['pka_type']}"
        try:
            sid = b.add_state(r['smiles'], 'iupac_v2_4')
            pid = b.add_parent(sid)
        except cn.IdentityError as e:
            b.drop('iupac_v2_4', rec, e); continue
        b.obs.append(dict(task='macro_pka', source='iupac_v2_4(sensitive_set_v2)', source_record=rec, role='curated_eval',
                          parent_id=pid, value=float(r['pka_median']), units='pKa (log10 units)',
                          convention='macroscopic pKa; acid/base microstate sets unknown; index ' + r['pka_type'],
                          temperature='20-25 C (filter)', solvent='water (no cosolvent)', label_status='measured',
                          states_A='', states_B='', reported_state=sid, _linked=[pid],
                          quality=json.dumps(dict(n_measurements=int(r['n_measurements']), spread=float(r['spread']),
                                                  any_reliable=int(r['any_reliable']), pka_type=r['pka_type'],
                                                  legacy_split=r['split'], stratum_v2=r['stratum_v2']))))


def unipka_sets(b, ext):
    for name, role in UNIPKA_SETS.items():
        path = os.path.join(ext, 'unipka_orig', 'dataset', f'{name}.tsv')
        for k, r in enumerate(csv.DictReader(open(path), delimiter='\t')):
            rec = f'{name}:{k}'
            trans = r.get('SMILES') or ''
            try:
                lhs, rhs = trans.split('>>')
                A = [b.add_state(s, name) for s in lhs.split(',') if s]
                B = [b.add_state(s, name) for s in rhs.split(',') if s]
                pids = {b.add_parent(s) for s in A + B}
            except (cn.IdentityError, ValueError) as e:
                b.drop(name, rec, e); continue
            qa = {b.states[s]['formal_charge'] for s in A}
            qb = {b.states[s]['formal_charge'] for s in B}
            pid = b.states[A[0]]['parent_id']
            if len(qa) != 1 or len(qb) != 1 or qa.pop() - qb.pop() != 1:
                b.drop(name, rec, f'acid/base sets not one proton apart in charge: {sorted(qa)} vs {sorted(qb)}'); continue
            b.obs.append(dict(task='macro_pka', source=f'unipka:{name}', source_record=rec, role=role, parent_id=pid,
                              value=float(r['TARGET']), units='pKa (log10 units)',
                              convention='macroscopic pKa of {A} <=> {B} + H+', temperature='not recorded',
                              solvent='water', label_status='measured', states_A=';'.join(A), states_B=';'.join(B),
                              reported_state='', _linked=sorted(pids),
                              quality=json.dumps(dict(ref=r.get('ref.', ''), n_parent_forms=len(pids)))))


def tautobase(b, ext):
    path = os.path.join(ext, 'tautobase', 'Tautobase_SMIRKS.txt')
    rows, keep = build_tautobase(path, kinds=MEASURED + ('e', 'c'))
    status = {'M': 'measured', 'E': 'estimated', 'C': 'calculated'}
    for r in rows:
        rec = f"row {r['row_nr']}"
        try:
            s1, s2 = b.add_state(r['smiles_t1'], 'tautobase'), b.add_state(r['smiles_t2'], 'tautobase')
            p1, p2 = b.add_parent(s1), b.add_parent(s2)
        except cn.IdentityError as e:
            b.drop('tautobase', rec, e); continue
        b.obs.append(dict(task='tautomer_logk', source='tautobase', source_record=rec, role='curated_eval', parent_id=p1,
                          value=float(r['log_k']), units='log10 K (dimensionless)',
                          convention='log10 K = log10([T2]/[T1]); T1 = states_A, T2 = states_B',
                          temperature=r['temperature'] or 'not stated', solvent=r['solvent'],
                          label_status=status[r['label_kind']], states_A=s1, states_B=s2, reported_state='',
                          _linked=sorted({p1, p2}),
                          quality=json.dumps(dict(reference=r['reference'], same_family=p1 == p2,
                                                  parent_id_T2=p2, legacy_split=r['split']))))
    return keep


def assign_splits(b):
    """Split groups are merged (union-find) whenever one observation links several parents (a tautomer pair whose
    tautomers standardise to different families, or a Uni-pKa transition with several parent forms), so no linked
    family can straddle train/val/test. A curated observation is excluded when any family it links also occurs in an
    external test set (family level, not component level: a component can be a whole common scaffold, e.g. benzene)."""
    up = {}

    def find(x):
        up.setdefault(x, x)
        while up[x] != x:
            up[x] = up[up[x]]
            x = up[x]
        return x
    for o in b.obs:
        gs = [b.parents[p]['split_group'] for p in o['_linked']]
        for g in gs[1:]:
            ra, rb = find(gs[0]), find(g)
            if ra != rb:
                up[max(ra, rb)] = min(ra, rb)
    comp = lambda g: find(g)  # noqa: E731
    external = {b.parents[p]['family_key'] for o in b.obs if o['role'] == 'external_test' for p in o['_linked']}
    for p in b.parents.values():
        p['split_component'] = comp(p['split_group'])
        p['split_curated'] = cn.assign_split(p['split_component'])
    for k, o in enumerate(b.obs):
        p = b.parents[o['parent_id']]
        o['obs_id'] = f'O{k:06d}'
        o['family_key'], o['scaffold'], o['split_group'] = p['family_key'], p['scaffold'], p['split_component']
        o['linked_parents'] = ';'.join(o['_linked'])
        if o['role'] == 'curated_eval':
            linked_fams = {b.parents[x]['family_key'] for x in o['_linked']}
            o['split'] = 'excluded_external_overlap' if linked_fams & external else p['split_curated']
        else:
            o['split'] = o['role']
    b.n_merged_components = sum(1 for g in {p['split_group'] for p in b.parents.values()} if find(g) != g)


def write_csv(path, rows, cols):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ext', required=True)
    ap.add_argument('--pka_curated', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    b = Builder()
    pka_curated(b, a.pka_curated)
    unipka_sets(b, a.ext)
    taut_keep = tautobase(b, a.ext)
    assign_splits(b)
    for s in b.states.values():
        s['origins'] = ';'.join(sorted(s['origins']))
    ocols = ['obs_id', 'task', 'source', 'source_record', 'role', 'split', 'parent_id', 'family_key', 'scaffold',
             'split_group', 'linked_parents', 'value', 'units', 'convention', 'temperature', 'solvent', 'label_status', 'states_A',
             'states_B', 'reported_state', 'quality']
    scols = ['state_id', 'smiles', 'formal_charge', 'n_heavy', 'n_charged_atoms', 'elements', 'n_unspecified_stereo',
             'parent_id', 'origins']
    pcols = ['parent_id', 'family_key', 'parent_smiles', 'scaffold', 'split_group', 'split_component', 'split_curated',
             'parent_method',
             'n_heavy', 'elements']
    write_csv(os.path.join(a.out, 'observations.csv'), b.obs, ocols)
    write_csv(os.path.join(a.out, 'states.csv'), sorted(b.states.values(), key=lambda s: s['state_id']), scols)
    write_csv(os.path.join(a.out, 'parents.csv'), sorted(b.parents.values(), key=lambda p: p['parent_id']), pcols)
    write_csv(os.path.join(a.out, 'dropped.csv'), b.dropped, ['source', 'source_record', 'reason'])
    inputs = [a.pka_curated, os.path.join(a.ext, 'tautobase', 'Tautobase_SMIRKS.txt')] + [
        os.path.join(a.ext, 'unipka_orig', 'dataset', f'{n}.tsv') for n in UNIPKA_SETS]
    outputs = {f: sha256(os.path.join(a.out, f)) for f in ('observations.csv', 'states.csv', 'parents.csv', 'dropped.csv')}
    manifest = dict(schema_version=cn.SCHEMA_VERSION, rdkit=rdkit.__version__, python=platform.python_version(),
                    split_policy='sha256(jmcg-split-v1|component) % 100: <70 train, <85 val, else test; '
                                 'split_group = Murcko scaffold of parent form, or acyclic:<family_key>; component = '
                                 'union of split groups linked by any observation (union-find); obs.split_group holds '
                                 'the component',
                    inputs={p: sha256(p) for p in inputs}, outputs=outputs, tautobase_filter=dict(taut_keep),
                    counts=dict(observations=len(b.obs), states=len(b.states), parents=len(b.parents),
                                dropped=len(b.dropped), state_parent_conflicts=len(b.state_parent_conflicts),
                                split_groups_merged_by_links=b.n_merged_components,
                                obs_by_task_role_split={f'{t}|{r}|{s}': n for (t, r, s), n in sorted(
                                    Counter((o['task'], o['role'], o['split']) for o in b.obs).items())}))
    json.dump(manifest, open(os.path.join(a.out, 'manifest.json'), 'w'), indent=1)
    print(json.dumps(manifest['counts'], indent=1))


if __name__ == '__main__':
    main()
