"""Tests for the sequential-baseline state stage (math needs no model; enumeration needs the Uni-pKa enumerator).
Run: python -m jmcg.seq.test_seq"""
import math

import numpy as np

from jmcg.seq import enumerate as en
from jmcg.seq import score as sc

RESULTS = []


def check(name, fn):
    try:
        RESULTS.append((name, 'PASS', fn() or ''))
    except Exception as e:  # noqa: BLE001
        RESULTS.append((name, 'FAIL', f'{type(e).__name__}: {e}'[:300]))


def main():
    def two_level():
        g = 5.0  # acid G=0, base G=g (kT): pKa = g/ln10 + PKA_MEAN; equal populations exactly there
        pk = sc.macro_pkas([0.0, g], [0, -1])[(0, -1)]
        assert abs(pk - (g / math.log(10) + sc.PKA_MEAN)) < 1e-12
        p = sc.populations([0.0, g], [0, -1], pk)
        assert abs(p[0] - 0.5) < 1e-12
        lo, hi = sc.populations([0.0, g], [0, -1], pk - 1), sc.populations([0.0, g], [0, -1], pk + 1)
        assert abs(lo[0] - 10 / 11) < 1e-12 and abs(hi[0] - 1 / 11) < 1e-12, (lo, hi)
        return f'pKa {pk:.3f}; acid fraction 10/11 at pKa-1, 1/2 at pKa, 1/11 at pKa+1 (Henderson-Hasselbalch)'
    check('populations and macro-pKa agree with Henderson-Hasselbalch for one transition', two_level)

    def microstates():
        # two equal acid microstates: macro pKa shifts by log10(2) relative to one (statistical factor)
        one = sc.macro_pkas([0.0, 3.0], [1, 0])[(1, 0)]
        two = sc.macro_pkas([0.0, 0.0, 3.0], [1, 1, 0])[(1, 0)]
        assert abs((two - one) - math.log10(2)) < 1e-12
        p = sc.populations([0.0, 0.0, 3.0], [1, 1, 0], 7.0)
        assert abs(p.sum() - 1) < 1e-12 and abs(p[0] - p[1]) < 1e-12
        return 'degenerate microstates add log10(2); populations normalised'
    check('microstate sums (LSE) give the statistical factor', microstates)

    def taut():
        assert abs(sc.log10_k(2.0, 2.0)) < 1e-12 and sc.log10_k(3.0, 0.0) > 0
        return 'log10([T2]/[T1]) = (G_T1 - G_T2)/ln10: lower G_T2 -> positive'
    check('tautomer ratio sign convention', taut)

    def enum_basic():
        st, f = en.enumerate_ensemble('CC(=O)O')
        smis = {s['smiles']: s['charge'] for s in st}
        assert smis.get('CC(=O)[O-]') == -1 and smis.get('CC(=O)O') == 0, smis
        st, f = en.enumerate_ensemble('NCC(=O)O')
        got = {s['smiles'] for s in st}
        need = {'NCC(=O)O', '[NH3+]CC(=O)[O-]', '[NH3+]CC(=O)O', 'NCC(=O)[O-]'}
        assert need <= got, (need - got, got)
        return f'acetic acid {sorted(smis)}; glycine {len(got)} states incl. zwitterion, cation, anion'
    check('enumerator: acetic acid and glycine ensembles', enum_basic)

    def enum_bounds():
        inc, _ = en.enumerate_ensemble('OC(=O)C(=O)O')
        uni, _ = en.enumerate_ensemble('OC(=O)C(=O)O', faithful_bounds=True)
        qi, qu = {s['charge'] for s in inc}, {s['charge'] for s in uni}
        assert -2 in qi and -2 not in qu, (qi, qu)
        return f'oxalic acid charges inclusive {sorted(qi)} vs Uni-pKa bounds {sorted(qu)}'
    check('inclusive charge window reaches -2 (Uni-pKa enum_ensemble stops at -1)', enum_bounds)

    def enum_valid():
        st, f = en.enumerate_ensemble('NCCc1c[nH]cn1')
        assert all(s['charge'] == s['level'] for s in st) and len({s['state_id'] for s in st}) == len(st)
        assert any(s['is_input'] for s in st)
        return f'histamine: {len(st)} states, charges {sorted({s["charge"] for s in st})}, {len(f)} failures'
    check('enumerated states are validated, unique, and include the input', enum_valid)

    n_fail = sum(r[1] == 'FAIL' for r in RESULTS)
    for r in RESULTS:
        print(f'{r[1]}  {r[0]}' + (f'  -> {r[2]}' if r[2] else ''))
    print(f'test_seq: {len(RESULTS) - n_fail} passed, {n_fail} failed')
    raise SystemExit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
