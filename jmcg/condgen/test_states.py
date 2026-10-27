"""Unit checks for jmcg.condgen.states (RDKit only; no model). Run: python -m jmcg.condgen.test_states"""
from jmcg.condgen import states as st

RESULTS = []


def check(name, fn):
    try:
        fn()
        RESULTS.append((name, 'PASS', ''))
    except Exception as e:  # noqa: BLE001
        RESULTS.append((name, 'FAIL', f'{type(e).__name__}: {e}'[:300]))


def rejects(state, parent, needle):
    def fn():
        try:
            st.build_state(state, parent)
        except st.StateError as e:
            assert needle in str(e), f'wrong reason: {e}'
            return
        raise AssertionError(f'accepted {state} vs {parent}')
    return fn


def accepts(state, parent, cls, **expect):
    def fn():
        m, i = st.build_state(state, parent)
        assert i['state_class'] == cls, i['state_class']
        for k, v in expect.items():
            got = len(i[k[2:]]) if k.startswith('n_') and k[2:] in i else i[k]
            assert got == v, f'{k}: {got} != {v}'
    return fn


def main():
    check('acetate = deprotonation, 1 site, 7 atoms', accepts('CC(=O)[O-]', 'CC(=O)O', 'protonation',
                                                             net_charge=-1, n_changed_sites=1, n_atoms_with_H=7))
    check('2-pyridone = tautomer of 2-hydroxypyridine', accepts('O=c1cccc[nH]1', 'Oc1ccccn1', 'tautomer',
                                                               net_charge=0, n_changed_sites=2))
    check('glycine zwitterion = zwitterionic tautomer', accepts('[NH3+]CC(=O)[O-]', 'NCC(=O)O', 'tautomer_zwitterionic',
                                                               net_charge=0, n_charged_atoms=2))
    check('amidinium: symmetric skeleton -> 1-site mapping', accepts('CC(N)=[NH2+]', 'CC(=N)N', 'protonation',
                                                                     n_changed_sites=1))
    check('histamine cation from N-pi tautomer = tautomer+protonation',
          accepts('[NH3+]CCc1c[nH]cn1', 'NCCc1cnc[nH]1', 'tautomer+protonation', net_charge=1))
    check('(S)-ibuprofen anion keeps parent stereocentre',
          accepts('CC(C)Cc1ccc(cc1)[C@H](C)C(=O)[O-]', 'CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O', 'protonation',
                  parent_stereocentres_checked=[10]))
    check('Z-enol of acetylacetone accepted (stereo specified)', accepts('C/C(O)=C/C(C)=O', 'CC(=O)CC(C)=O', 'tautomer'))
    check('reject redox: ethanol -> acetaldehyde', rejects('CC=O', 'CCO', 'proton transfer'))
    check('reject redox: H2 added', rejects('CCCC', 'C=CCC', 'proton transfer'))
    check('reject charge moved without proton', rejects('C[CH+][O-]', 'CC=O', 'charge moved without a proton'))
    check('reject carbanion/carbocation site', rejects('C[CH2-]', 'CC', 'charge change on C'))
    check('reject different skeleton', rejects('CCN', 'CCO', 'skeleton differs'))
    check('reject salt', rejects('CC(=O)[O-].[Na+]', 'CC(=O)O', 'fragments'))
    check('reject unparsable', rejects('C1CC1(', None, 'does not parse'))
    check('reject valence error', rejects('C[N](C)(C)(C)C', None, 'does not parse'))
    check('reject radical', rejects('C[CH2]', None, 'radical'))
    check('reject unassigned stereocentre + double bond', rejects('CC=CC(C)O', None, 'unassigned stereo'))
    check('reject inverted parent stereocentre', rejects('CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)[O-]',
                                                         'CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O', 'opposite configuration'))

    def enum():
        iso = st.expand_stereo('CC=CC(C)O')
        assert len(iso) == 4 and len(set(iso)) == 4, iso
        for s in iso:
            st.build_state(s, 'CC=CC(C)O')
    check('enumerate policy: 4 explicit stereoisomers, each accepted', enum)

    def imine():
        _, i = st.build_state('CC(=N)N', None)  # =NH E/Z is not a heavy-atom stereo element (documented limitation)
        assert not [x for x in i['stereo_elements'] if not x['specified']], i['stereo_elements']
    check('acetamidine: no heavy-atom stereo element (imine N-H E/Z not representable)', imine)

    n_fail = sum(r[1] == 'FAIL' for r in RESULTS)
    for r in RESULTS:
        print(f'{r[1]}  {r[0]}' + (f'  -> {r[2]}' if r[2] else ''))
    print(f'test_states: {len(RESULTS) - n_fail} passed, {n_fail} failed')
    raise SystemExit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
