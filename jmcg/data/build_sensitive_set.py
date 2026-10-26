"""Build the conformation-sensitivity-stratified pKa / tautomer-ratio test set (plan step 2).

Inputs (experimental labels only):
  IUPAC Digitized pKa v2.4a  (CC BY-NC 4.0)  ext/iupac_pka/iupac_v2_4.csv
  Tautobase SMIRKS file      (see its README) ext/tautobase/Tautobase_SMIRKS.txt
  SAMPL6 macro-pKa           ext/sampl6_pka/pKa_experimental_values.csv  (blind test only)

Outputs (in --out):
  pka_curated.csv        one row per (molecule, pKa index) with median value, n, spread, flags, stratum, split
  tautobase_curated.csv  one row per aqueous measured tautomer pair with flags, stratum, split
  manifest.json          filters, counts per stratum x split, input sha256s, versions

Stratum definitions are fixed here BEFORE any model is run (pre-registration):
  sensitive = flexible IMHB partner of an ionisable/tautomeric site, OR zwitterion-capable,
              OR two ionisable sites within <= 5 bonds
  control   = exactly one ionisable site and no H-bond partner within 6 bonds of it
  other     = everything else
The ionisable-site SMARTS list below is for stratification only; state enumeration uses Dimorphite-DL.
"""
import argparse
import csv
import hashlib
import json
import os
import statistics
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog('rdApp.*')

ACID_SMARTS = {
    'carboxylic_acid': '[CX3](=O)[OX2H1,OX1-]',
    'phenol': 'c[OX2H1,OX1-]',
    'sulfonamide_NH': '[SX4](=O)(=O)[NX3;H1,H2]',
    'tetrazole': 'c1nn[nH]n1',
    'phosph_acid': '[PX4](=O)[OX2H1,OX1-]',
    'sulfonic_acid': '[SX4](=O)(=O)[OX2H1,OX1-]',
    'thiol': '[SX2H1]',
    'imide': '[CX3](=O)[NX3H1][CX3](=O)',
    'hydroxamic': '[CX3](=O)[NX3H1][OX2H1]',
}
BASE_SMARTS = {
    'aliphatic_amine': '[NX3;H2,H1,H0;!$(N-[C,S,P]=[O,S,N]);!$(N-a);!$(N#*);!$(N=*);!$(N-N);!$(N-O)]',
    'aniline': '[NX3;H2,H1,H0;$(N-a);!$(N-[C,S]=[O,S])]',
    'amidine_guanidine': '[CX3](=[NX2;!$(N-O)])[NX3]',
    'pyridine_like': '[nX2;r6]',
    'imidazole_like': '[nX2;r5]',
}
HBOND_PARTNER = Chem.MolFromSmarts('[O,N;!$([N+](=O)[O-])]')
ALLOWED = set('H C N O S P F Cl Br I'.split())


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def site_atoms(mol):
    acid, base = set(), set()
    for patts, out in ((ACID_SMARTS, acid), (BASE_SMARTS, base)):
        for s in patts.values():
            for match in mol.GetSubstructMatches(Chem.MolFromSmarts(s)):
                # the titratable heavy atom: last O/N/S in the match for acids, first N for bases
                cand = [i for i in match if mol.GetAtomWithIdx(i).GetSymbol() in 'ONS']
                if cand:
                    out.add(cand[-1] if out is acid else cand[0])
    return acid, base


def path_has_rotor(mol, path):
    for a, b in zip(path[:-1], path[1:]):
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond.GetBondType() == Chem.BondType.SINGLE and not bond.IsInRing():
            ends = [mol.GetAtomWithIdx(a), mol.GetAtomWithIdx(b)]
            if all(e.GetDegree() > 1 for e in ends):
                return True
    return False


def flags(mol, extra_sites=()):
    """Conformation-sensitivity flags for one heavy-atom molecule (H implicit)."""
    acid, base = site_atoms(mol)
    sites = acid | base | set(extra_sites)
    dm = Chem.GetDistanceMatrix(mol)
    partners = [m[0] for m in mol.GetSubstructMatches(HBOND_PARTNER)]
    imhb_flex = imhb_rigid = False
    near_partner = False
    for s in sites:
        for p in partners:
            d = dm[s][p]
            if p == s or d < 3:
                continue
            if d <= 6:
                near_partner = True
            if 3 <= d <= 6:  # pseudo-ring of 5-8 atoms including the bridging H
                path = Chem.GetShortestPath(mol, s, p)
                if path_has_rotor(mol, path):
                    imhb_flex = True
                else:
                    imhb_rigid = True
    site_list = sorted(sites)
    close_pair = any(dm[a][b] <= 5 for i, a in enumerate(site_list) for b in site_list[i + 1:])
    zwit = bool(acid) and bool(base)
    sensitive = imhb_flex or zwit or close_pair
    control = len(sites) == 1 and not near_partner
    return {
        'n_acid_sites': len(acid), 'n_base_sites': len(base), 'n_sites': len(sites),
        'imhb_flexible': int(imhb_flex), 'imhb_rigid': int(imhb_rigid), 'zwitterion_capable': int(zwit),
        'close_site_pair': int(close_pair), 'n_rot': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'stratum': 'sensitive' if sensitive else ('control' if control else 'other'),
    }


STRONG_PARTNER = Chem.MolFromSmarts(
    '[$([OX1]=[#6,#16,#15]),$([OX2H1]),$([OX1-]),$([NX3;H1,H2;!$(N-S(=O)=O)]),$([nX2]),'
    '$([NX3H0;!$(N-C=[O,S]);!$(N-a);!$(N-S(=O)=O)])]')


def flags_v2(mol, extra_sites=()):
    """Stricter, through-space proxy (v2; fixed 2026-10-26 before any model output — v1 flagged ~75 % of molecules).

    sensitive_v2 = a 6/7-membered IMHB pseudo-ring (site–partner 4–5 bonds) to a strong donor/acceptor via a rotor,
                   OR acid+base sites 3–8 bonds apart via a rotor (intramolecular salt bridge possible),
                   OR two ionisable sites 4–8 bonds apart via a rotor (through-space site–site coupling).
    """
    acid, base = site_atoms(mol)
    sites = sorted(acid | base | set(extra_sites))
    dm = Chem.GetDistanceMatrix(mol)
    partners = [m[0] for m in mol.GetSubstructMatches(STRONG_PARTNER)]

    def flexible(a, b):
        return path_has_rotor(mol, Chem.GetShortestPath(mol, a, b))

    imhb = any(4 <= dm[s][p] <= 5 and flexible(s, p) for s in sites for p in partners if p != s)
    salt = any(3 <= dm[x][y] <= 8 and flexible(x, y) for x in acid for y in base)
    pair = any(4 <= dm[x][y] <= 8 and flexible(x, y) for i, x in enumerate(sites) for y in sites[i + 1:])
    near = any(dm[s][p] <= 6 for s in sites for p in partners if p != s and dm[s][p] >= 3)
    stratum = 'sensitive' if (imhb or salt or pair) else ('control' if len(sites) == 1 and not near else 'other')
    return {'imhb6_7_flexible_v2': int(imhb), 'salt_bridge_flexible_v2': int(salt),
            'site_pair_flexible_v2': int(pair), 'stratum_v2': stratum}


def scaffold_split(mol, train=70, val=15):
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        scaf = ''
    key = scaf or Chem.MolToSmiles(mol)
    bucket = int(hashlib.sha1(key.encode()).hexdigest(), 16) % 100
    return scaf, 'train' if bucket < train else ('val' if bucket < train + val else 'test')


def ok_mol(mol):
    return (mol is not None and all(a.GetSymbol() in ALLOWED for a in mol.GetAtoms())
            and mol.GetNumHeavyAtoms() <= 50 and len(Chem.GetMolFrags(mol)) == 1)


def build_iupac(path, sampl_keys):
    keep = Counter()
    groups = defaultdict(list)
    meta = {}
    with open(path, encoding='utf-8', errors='replace', newline='') as fh:
        for r in csv.DictReader(fh):
            keep['rows_in'] += 1
            if r['assessment'] not in ('Reliable', 'Approximate'):
                keep['drop_assessment'] += 1; continue
            if r['cosolvent'].strip():
                keep['drop_cosolvent'] += 1; continue
            try:
                t = float(r['T']); v = float(r['pka_value'])
            except ValueError:
                keep['drop_T_or_value_nonnumeric'] += 1; continue
            if not 20 <= t <= 25:
                keep['drop_T_range'] += 1; continue
            if not r['pka_type'].startswith(('pKa', 'pKaH')):
                keep['drop_pka_type'] += 1; continue
            mol = Chem.MolFromSmiles(r['SMILES'])
            if not ok_mol(mol):
                keep['drop_mol'] += 1; continue
            ik = Chem.MolToInchiKey(mol)
            groups[(ik, r['pka_type'])].append((v, r['assessment']))
            meta[ik] = Chem.MolToSmiles(mol)
            keep['rows_kept'] += 1
    rows = []
    for (ik, ptype), vals in sorted(groups.items()):
        mol = Chem.MolFromSmiles(meta[ik])
        values = [v for v, _ in vals]
        f = {**flags(mol), **flags_v2(mol)}
        scaf, split = scaffold_split(mol)
        if ik in sampl_keys:
            split = 'excluded_sampl_overlap'
        rows.append({'inchikey': ik, 'smiles': meta[ik], 'pka_type': ptype,
                     'pka_median': round(statistics.median(values), 3), 'n_measurements': len(values),
                     'spread': round(max(values) - min(values), 3),
                     'any_reliable': int(any(a == 'Reliable' for _, a in vals)),
                     'heavy_atoms': mol.GetNumHeavyAtoms(), 'scaffold': scaf, 'split': split, **f})
    return rows, keep


def build_tautobase(path):
    keep = Counter()
    rows = []
    with open(path, encoding='utf-8', errors='replace', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            keep['rows_in'] += 1
            solv = (r.get('State/Solvent') or '').strip().lower()
            if not any(w in solv for w in ('water', 'aqueous', 'h2o')):
                keep['drop_not_aqueous'] += 1; continue
            if (r.get('Meas/Calc/Est') or '').strip().lower() not in ('', 'meas', 'measured', 'm'):
                keep['drop_not_measured'] += 1; continue
            try:
                logk = float(r['log K'])
            except (ValueError, TypeError):
                keep['drop_logK_nonnumeric'] += 1; continue
            try:
                lhs, rhs = r['Smirks'].split('>>')
                m1, m2 = Chem.MolFromSmiles(lhs), Chem.MolFromSmiles(rhs)
                for m in (m1, m2):
                    for a in m.GetAtoms():
                        a.SetAtomMapNum(0)
                    Chem.SanitizeMol(m)
            except Exception:
                keep['drop_parse'] += 1; continue
            if not (ok_mol(m1) and ok_mol(m2)):
                keep['drop_mol'] += 1; continue
            # tautomeric heavy atoms = those whose H count changes
            moving = [a.GetIdx() for a, b in zip(m1.GetAtoms(), m2.GetAtoms())
                      if a.GetTotalNumHs() != b.GetTotalNumHs() and a.GetSymbol() in 'NOS']
            f = {**flags(m1, extra_sites=moving), **flags_v2(m1, extra_sites=moving)}
            scaf, split = scaffold_split(m1)
            rows.append({'row_nr': r.get('Row Nr'), 'smiles_t1': Chem.MolToSmiles(m1),
                         'smiles_t2': Chem.MolToSmiles(m2), 'log_k': logk,
                         'solvent': r.get('State/Solvent'), 'temperature': r.get('Temperature'),
                         'reference': r.get('Reference'), 'heavy_atoms': m1.GetNumHeavyAtoms(),
                         'scaffold': scaf, 'split': split, **f})
            keep['rows_kept'] += 1
    return rows, keep


def sampl_inchikeys(path):
    keys = set()
    if not os.path.exists(path):
        return keys
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh):
            smi = r.get('canonical isomeric SMILES') or r.get('SMILES') or r.get('smiles')
            m = Chem.MolFromSmiles(smi) if smi else None
            if m is not None:
                keys.add(Chem.MolToInchiKey(m))
    return keys


def unipka_skeleton_keys(path):
    """Protonation-invariant keys (InChIKey first block) of every molecule in a Uni-pKa reaction-SMILES TSV."""
    keys = set()
    if not os.path.exists(path):
        return keys
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            for side in (r.get('SMILES') or '').replace('>>', ',').split(','):
                m = Chem.MolFromSmiles(side)
                if m is not None:
                    keys.add(Chem.MolToInchiKey(m)[:14])
    return keys


def write(rows, path):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ext', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    iupac = os.path.join(a.ext, 'iupac_pka/iupac_v2_4.csv')
    taut = os.path.join(a.ext, 'tautobase/Tautobase_SMIRKS.txt')
    sampl = os.path.join(a.ext, 'sampl6_pka/pKa_experimental_values.csv')
    skeys = sampl_inchikeys(sampl)
    p_rows, p_keep = build_iupac(iupac, skeys)
    t_rows, t_keep = build_tautobase(taut)
    # leakage flags vs Uni-pKa: fine-tuning set (Dwar-iBonD) and its reported test sets
    upk = os.path.join(a.ext, 'unipka_orig/dataset')
    train_keys = unipka_skeleton_keys(os.path.join(upk, 'dwar-iBond.tsv'))
    test_keys = set().union(*(unipka_skeleton_keys(os.path.join(upk, f)) for f in
                              ('sampl6.tsv', 'sampl7.tsv', 'sampl8.tsv', 'novartis_acid.tsv', 'novartis_base.tsv')))
    for r in p_rows:
        r['in_unipka_finetune'] = int(r['inchikey'][:14] in train_keys)
        r['in_unipka_testsets'] = int(r['inchikey'][:14] in test_keys)
    for r in t_rows:
        k = Chem.MolToInchiKey(Chem.MolFromSmiles(r['smiles_t1']))[:14]
        r['in_unipka_finetune'] = int(k in train_keys)
    write(p_rows, os.path.join(a.out, 'pka_curated.csv'))
    write(t_rows, os.path.join(a.out, 'tautobase_curated.csv'))

    def table(rows):
        c = Counter((r['stratum'], r['split']) for r in rows)
        return {f'{s}/{sp}': n for (s, sp), n in sorted(c.items())}
    manifest = {
        'inputs': {p: sha256(p) for p in (iupac, taut) + ((sampl,) if os.path.exists(sampl) else ())},
        'rdkit': rdBase.rdkitVersion,
        'filters_iupac': dict(p_keep), 'filters_tautobase': dict(t_keep),
        'pka_rows': len(p_rows), 'pka_molecules': len({r['inchikey'] for r in p_rows}),
        'pka_stratum_x_split': table(p_rows),
        'pka_stratum_v2_x_split': {f'{a}/{b}': n for (a, b), n in sorted(Counter((r['stratum_v2'], r['split']) for r in p_rows).items())},
        'pka_stratum_v2_clean': dict(Counter(r['stratum_v2'] for r in p_rows if not r['in_unipka_finetune'])),
        'pka_stratum_v2_clean_reliable': dict(Counter(r['stratum_v2'] for r in p_rows if not r['in_unipka_finetune'] and r['any_reliable'])),
        'tautobase_stratum_v2_clean': dict(Counter(r['stratum_v2'] for r in t_rows if not r['in_unipka_finetune'])),
        'pka_reliable_stratum': dict(Counter(r['stratum'] for r in p_rows if r['any_reliable'])),
        'pka_inconsistent_spread_gt1': sum(r['spread'] > 1 for r in p_rows),
        'tautobase_rows': len(t_rows), 'tautobase_stratum_x_split': table(t_rows),
        'sampl6_overlap_excluded': sum(r['split'] == 'excluded_sampl_overlap' for r in p_rows),
        'unipka_keys': {'dwar_ibond': len(train_keys), 'testsets': len(test_keys)},
        'pka_clean_of_unipka_finetune_by_stratum': dict(Counter(r['stratum'] for r in p_rows
                                                                 if not r['in_unipka_finetune'])),
        'pka_clean_reliable_by_stratum': dict(Counter(r['stratum'] for r in p_rows
                                                      if not r['in_unipka_finetune'] and r['any_reliable'])),
        'tautobase_clean_by_stratum': dict(Counter(r['stratum'] for r in t_rows if not r['in_unipka_finetune'])),
        'stratum_flag_counts_pka': {k: sum(int(r[k]) for r in p_rows) for k in
                                    ('imhb_flexible', 'imhb_rigid', 'zwitterion_capable', 'close_site_pair')},
    }
    with open(os.path.join(a.out, 'manifest.json'), 'w') as fh:
        json.dump(manifest, fh, indent=1)
    print(json.dumps(manifest, indent=1), flush=True)


if __name__ == '__main__':
    main()
