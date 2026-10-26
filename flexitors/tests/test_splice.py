"""flexitors.sampling.splice (EXP-3 §9.1)."""
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from flexitors.sampling.splice import prefix, splice, unaffected_deviation


def mols(smi, seed):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    return [m]


def test_splice_and_identity():
    o = {'A': mols('CCO', 1), 'B': mols('CC=CC', 1), 'C': mols('CCN', 1)}
    arm = {'B': mols('CC=CC', 2)}
    s = splice(o, arm, {'B', 'C'})  # C failed in the arm: it stays a failure
    assert s['A'] is o['A'] and s['B'] is arm['B'] and 'C' not in s
    with pytest.raises(ValueError):
        splice(o, {'A': o['A']}, {'B'})
    assert unaffected_deviation(o, {'A': mols('CCO', 1)}) == 0.0
    assert unaffected_deviation(o, {'B': arm['B']}) > 0
    with pytest.raises(ValueError):
        unaffected_deviation(o, {'Z': arm['B']})


def test_prefix_is_the_first_2k():
    pool = {'A': list(range(10)), 'B': list(range(4))}
    assert prefix(pool, {'A': 2, 'B': 2}) == {'A': [0, 1, 2, 3], 'B': [0, 1, 2, 3]}
    with pytest.raises(ValueError):
        prefix(pool, {'A': 2, 'B': 3})
