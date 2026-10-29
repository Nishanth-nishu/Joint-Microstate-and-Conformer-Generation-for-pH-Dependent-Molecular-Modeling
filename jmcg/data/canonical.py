"""Canonical identifiers for parents (molecular families), microstates and observations. Schema version 1.

Policy (documented in research/experiments/M2_data/README.md):
  microstate  = the exact graph: heavy atoms, bonds/bond orders, formal charges, H counts (implied by SMILES valence),
                aromaticity as perceived by RDKit, and every SPECIFIED stereo element. Canonical isomeric SMILES of the
                sanitised single-fragment molecule; state_id = 'S' + sha1(smiles)[:16]. Distinct protonation states,
                tautomers and stereoisomers therefore always get distinct ids; unspecified stereo stays unspecified.
  parent      = the molecular family shared by all protomers/tautomers of one compound: the RDKit-standardised form
                (largest fragment -> uncharged where possible -> canonical tautomer), stereo kept as given.
                parent_id = 'P' + InChIKey of that form. family_key = first InChIKey block of the same form
                (connectivity only: groups stereoisomers too) and is what splits are grouped by.
                Canonical-tautomer standardisation is a heuristic; when two members of one measured tautomer pair get
                different family keys the audit reports it and the pair keeps the key of its first tautomer.
  scaffold    = Murcko scaffold SMILES of the parent form (charge/tautomer independent); acyclic -> '' and those are
                grouped by family_key.
Salts / several fragments are rejected (not stripped) for microstates; RDKit version is recorded in every manifest.
"""
import hashlib

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog('rdApp.*')
SCHEMA_VERSION = 1

_UNCHARGER = rdMolStandardize.Uncharger()
_TP = rdMolStandardize.CleanupParameters()
_TP.tautomerRemoveSp3Stereo = False   # RDKit's default strips stereo at tautomer-active atoms: (R)/(S)-alanine
_TP.tautomerRemoveBondStereo = False  # collapsed into one parent (test_canonical, job 45290)
_TAUT = rdMolStandardize.TautomerEnumerator(_TP)


class IdentityError(ValueError):
    pass


def microstate(smiles):
    """(canonical isomeric SMILES, Mol) of one explicit microstate; raises IdentityError with the reason."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise IdentityError(f'does not parse/sanitise: {smiles}')
    if len(Chem.GetMolFrags(m)) > 1:
        raise IdentityError(f'several fragments (salt/mixture), not stripped: {smiles}')
    if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
        raise IdentityError(f'radical: {smiles}')
    return Chem.MolToSmiles(m), m


def state_id(canonical_smiles):
    return 'S' + hashlib.sha1(canonical_smiles.encode()).hexdigest()[:16]


def parent_form(mol):
    """Uncharged canonical tautomer of a (single-fragment) microstate, stereo as given, with a recorded fallback:
    each candidate must round-trip through SMILES; else the uncharged form; else the microstate itself."""
    cands = []
    try:
        u = _UNCHARGER.uncharge(Chem.Mol(mol))
    except Exception:  # noqa: BLE001
        u = None
    if u is not None:
        try:
            cands.append((_TAUT.Canonicalize(u), 'canonical_tautomer'))
        except Exception:  # noqa: BLE001
            pass
        cands.append((u, 'uncharged_only'))
    cands.append((Chem.Mol(mol), 'microstate_as_is'))
    for m, how in cands:
        smi = Chem.MolToSmiles(m)
        if Chem.MolFromSmiles(smi) is not None:
            return m, how
    raise IdentityError(f'no parent form round-trips for {Chem.MolToSmiles(mol)}')


def parent_identity(mol):
    """dict(parent_id, family_key, parent_smiles, scaffold, parent_method)."""
    p, how = parent_form(mol)
    ik = Chem.MolToInchiKey(p)
    if not ik:
        raise IdentityError(f'no InChIKey for parent form of {Chem.MolToSmiles(mol)}')
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=p)
    except Exception:  # noqa: BLE001
        scaf = ''
    return dict(parent_id='P' + ik, family_key=ik.split('-')[0], parent_smiles=Chem.MolToSmiles(p), scaffold=scaf,
                parent_method=how)


def split_group(identity):
    """Split grouping key: scaffold if cyclic, else the family (so acyclic families never straddle splits)."""
    return identity['scaffold'] or 'acyclic:' + identity['family_key']


def assign_split(group, train=70, val=15, salt='jmcg-split-v1'):
    """Deterministic hash bucket of the group key (no RNG state, stable across runs and machines)."""
    b = int(hashlib.sha256(f'{salt}|{group}'.encode()).hexdigest(), 16) % 100
    return 'train' if b < train else 'val' if b < train + val else 'test'
