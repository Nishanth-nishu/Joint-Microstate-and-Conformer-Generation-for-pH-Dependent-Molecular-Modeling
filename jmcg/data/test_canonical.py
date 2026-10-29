"""Tests for the canonical identity policy (jmcg.data.canonical). Run: python -m jmcg.data.test_canonical"""
from jmcg.data import canonical as cn

RESULTS = []


def check(name, fn):
    try:
        RESULTS.append((name, 'PASS', fn() or ''))
    except Exception as e:  # noqa: BLE001
        RESULTS.append((name, 'FAIL', f'{type(e).__name__}: {e}'[:300]))


def ident(smi):
    can, mol = cn.microstate(smi)
    return cn.state_id(can), cn.parent_identity(mol)


def same_family_distinct_states(pairs):
    out = []
    for a, b in pairs:
        (sa, pa), (sb, pb) = ident(a), ident(b)
        assert sa != sb, (a, b, 'state ids collide')
        assert pa['parent_id'] == pb['parent_id'], (a, b, pa['parent_smiles'], pb['parent_smiles'])
        out.append(pa['parent_smiles'])
    return '; '.join(out)


def main():
    check('protomers share a parent, keep distinct state ids',
          lambda: same_family_distinct_states([('CC(=O)O', 'CC(=O)[O-]'), ('c1ccncc1', 'c1cc[nH+]cc1'),
                                               ('CCCN', 'CCC[NH3+]')]))
    check('zwitterion and neutral form share a parent', lambda: same_family_distinct_states([('NCC(=O)O', '[NH3+]CC(=O)[O-]')]))
    check('tautomers share a parent (canonical tautomer)',
          lambda: same_family_distinct_states([('Oc1ccccn1', 'O=c1cccc[nH]1'), ('CC(=O)CC(C)=O', r'C/C(O)=C/C(C)=O')]))

    def stereo():
        (s1, p1), (s2, p2) = ident('C[C@H](N)C(=O)O'), ident('C[C@@H](N)C(=O)O')
        (s3, p3) = ident('CC(N)C(=O)O')
        assert len({s1, s2, s3}) == 3 and p1['parent_id'] != p2['parent_id'], 'enantiomers must stay distinct'
        assert p1['family_key'] == p2['family_key'] == p3['family_key'], 'family key groups stereoisomers'
        (e, pe), (z, pz) = ident('C/C=C/C'), ident(r'C/C=C\C')
        assert e != z and pe['family_key'] == pz['family_key']
        return 'R/S/unspecified: 3 states, 2+1 parents, 1 family; E/Z: 2 states, 1 family'
    check('stereoisomers: distinct states and parents, one family key (split grouping)', stereo)

    def rejects():
        for bad in ['CC(=O)[O-].[Na+]', 'C1CC1(', '[CH3]']:
            try:
                cn.microstate(bad)
                raise AssertionError(f'accepted {bad}')
            except cn.IdentityError:
                pass
        return 'salt, unparsable, radical rejected (not stripped)'
    check('salts, parse errors and radicals are rejected with a reason', rejects)

    def charges_distinct():
        ids = {ident(s)[0] for s in ['NCCc1c[nH]cn1', 'NCCc1cnc[nH]1', '[NH3+]CCc1c[nH]cn1', '[NH3+]CCc1c[nH]c[nH+]1']}
        assert len(ids) == 4
        return 'histamine N-tau / N-pi / cation / dication: 4 distinct state ids'
    check('chemically distinct microstates never collapse', charges_distinct)

    def splits():
        _, p = ident('CC(=O)O')
        _, q = ident('c1ccccc1C(=O)O')
        g1, g2 = cn.split_group(p), cn.split_group(q)
        assert g1.startswith('acyclic:') and g2 == 'c1ccccc1', (g1, g2)
        assert cn.assign_split(g2) == cn.assign_split(g2) and cn.assign_split(g2) in ('train', 'val', 'test')
        c = [cn.assign_split(f'g{k}') for k in range(20000)]
        frac = {s: round(c.count(s) / len(c), 3) for s in ('train', 'val', 'test')}
        assert abs(frac['train'] - .7) < .02 and abs(frac['val'] - .15) < .02, frac
        return f'acyclic grouped by family; deterministic; fractions on 20k keys {frac}'
    check('split groups and deterministic 70/15/15 assignment', splits)

    n_fail = sum(r[1] == 'FAIL' for r in RESULTS)
    for r in RESULTS:
        print(f'{r[1]}  {r[0]}' + (f'  -> {r[2]}' if r[2] else ''))
    print(f'test_canonical: {len(RESULTS) - n_fail} passed, {n_fail} failed')
    raise SystemExit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
