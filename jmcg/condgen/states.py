"""Explicit chemical-state handling for conditional conformer generation.

A *state* is a fully specified microstate SMILES (formal charges, H placement, bond orders, stereo). It is never a label
on top of another graph: the RDKit Mol built here IS the graph that is featurised, embedded and written to SDF.

`build_state(state, parent)` accepts the state only if ALL of these hold, else raises StateError with the reason:
  1. state and parent parse and sanitise (valence, aromaticity); single fragment; no radicals; no isotopes;
  2. same heavy-atom skeleton (elements + connectivity, ignoring bond orders, charges, H);
  3. related by proton transfer only: total dH == total dq (state and parent have the same electron count once the
     added/removed H+ are accounted for). Redox changes (e.g. ethanol -> acetaldehyde, dH = -2, dq = 0) are rejected;
  4. per site: either the site gains/loses H+ together with its charge (dH == dq) or an H moves with no charge change
     (dq == 0). A charge that moves without a proton (charge-separated resonance form, e.g. C=O -> [C+][O-]) is rejected;
  5. charge changes only on heteroatoms (N, O, S, P, Se); carbon may gain/lose H only with dq == 0 (keto/enol);
  6. stereo: every potential stereo element of the state must be specified (policy "require"), or the caller expands
     the state into all its stereoisomers with `expand_stereo` and treats each as a separate explicit state (policy
     "enumerate"). No stereoisomer is ever picked silently. Stereocentres specified in the parent and still present
     in the state must keep their configuration (compared on mapped neighbour order, not CIP labels, so a CIP change
     caused by the protonation itself does not trigger a false rejection).
Atom mapping: the skeleton automorphism with the fewest changed sites (so an amidine protonation is one site, not an H
shift plus a charge change); with several equally minimal mappings the lexicographically first is used and the number
of ties is recorded.
"""
from collections import Counter

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

CHARGE_SITES = {'N', 'O', 'S', 'P', 'Se'}


class StateError(ValueError):
    pass


def _parse(smi, what):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise StateError(f'{what} SMILES does not parse/sanitise: {smi}')
    if len(Chem.GetMolFrags(m)) > 1:
        raise StateError(f'{what} has several fragments (salts/counter-ions not supported): {smi}')
    if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
        raise StateError(f'{what} has radical electrons: {smi}')
    if any(a.GetIsotope() for a in m.GetAtoms()):
        raise StateError(f'{what} has isotope labels: {smi}')
    return m


def skeleton(mol):
    """Heavy-atom graph with all bonds single, no charges, no H, no stereo: identical for every protomer/tautomer."""
    m = Chem.RWMol(Chem.RemoveHs(mol))
    for a in m.GetAtoms():
        a.SetFormalCharge(0); a.SetNumExplicitHs(0); a.SetNoImplicit(True); a.SetIsAromatic(False)
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False); b.SetStereo(Chem.BondStereo.STEREONONE)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def skeleton_key(mol):
    return Chem.MolToSmiles(skeleton(mol), canonical=True)


def stereo_elements(mol):
    out = []
    for s in Chem.FindPotentialStereo(mol):
        kind = 'atom' if s.type == Chem.StereoType.Atom_Tetrahedral else \
            'double_bond' if s.type == Chem.StereoType.Bond_Double else str(s.type).split('.')[-1]
        out.append(dict(kind=kind, centered_on=s.centeredOn,
                        specified=s.specified != Chem.StereoSpecified.Unspecified))
    return out


def expand_stereo(state_smiles):
    """All stereoisomers of the unassigned elements, as canonical isomeric SMILES (sorted). Explicit expansion for
    policy 'enumerate'; each one is then a separate explicit state."""
    m = _parse(state_smiles, 'state')
    if all(x['specified'] for x in stereo_elements(m)):
        return [Chem.MolToSmiles(m)]
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True, tryEmbedding=False)
    return sorted({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m, options=opts)})


def _changed_sites(p, s, match):
    sites = []
    for i_s, i_p in enumerate(match):
        a_s, a_p = s.GetAtomWithIdx(i_s), p.GetAtomWithIdx(i_p)
        dh, dq = a_s.GetTotalNumHs() - a_p.GetTotalNumHs(), a_s.GetFormalCharge() - a_p.GetFormalCharge()
        if dh or dq:
            sites.append(dict(parent_idx=i_p, state_idx=i_s, element=a_s.GetSymbol(), dH=dh, dq=dq))
    return sites


def site_diff(parent, state):
    """Map state heavy atoms onto parent heavy atoms via the shared skeleton and list what changed."""
    if skeleton_key(parent) != skeleton_key(state):
        raise StateError('not a protomer/tautomer of the parent: heavy-atom skeleton differs')
    sp, ss = skeleton(parent), skeleton(state)
    p, s = Chem.RemoveHs(parent), Chem.RemoveHs(state)
    matches = sp.GetSubstructMatches(ss, uniquify=False, maxMatches=5000)  # match[i_state] = i_parent
    if not matches:
        raise StateError('skeleton match failed')
    scored = sorted((len(_changed_sites(p, s, m)), m) for m in matches)
    match = list(scored[0][1])
    ties = sum(1 for n, _ in scored if n == scored[0][0])
    sites = _changed_sites(p, s, match)
    bonds = []
    for b in s.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bp = p.GetBondBetweenAtoms(match[i], match[j])
        if bp.GetBondType() != b.GetBondType():
            bonds.append(dict(parent_atoms=(match[i], match[j]), parent=str(bp.GetBondType()), state=str(b.GetBondType())))
    return match, sites, bonds, ties


def check_transformation(sites):
    dH = sum(x['dH'] for x in sites)
    dq = sum(x['dq'] for x in sites)
    if dH != dq:
        raise StateError(f'not related to the parent by proton transfer: total dH={dH:+d} but dq={dq:+d} '
                         '(electron count differs: redox or other change)')
    for x in sites:
        if x['dq'] and x['dH'] != x['dq']:
            raise StateError(f"charge moved without a proton at {x['element']}{x['parent_idx']} (dH={x['dH']:+d}, "
                             f"dq={x['dq']:+d}): resonance/charge-separated form, not a protomer/tautomer")
        if x['dq'] and x['element'] not in CHARGE_SITES:
            raise StateError(f"charge change on {x['element']}{x['parent_idx']}: only {sorted(CHARGE_SITES)} accepted")
    return dH, dq


def classify(sites, bonds, dq):
    if not sites:
        return 'parent' if not bonds else 'parent_resonance_form'
    if dq == 0:
        return 'tautomer_zwitterionic' if any(x['dq'] for x in sites) else 'tautomer'
    if all(x['dH'] == x['dq'] for x in sites):
        return 'protonation'
    return 'tautomer+protonation'


def _perm_parity(a, b):
    """Parity (0 even, 1 odd) of the permutation taking list a to list b (same elements)."""
    a, parity = list(a), 0
    for i in range(len(a)):
        j = a.index(b[i], i)
        if j != i:
            a[i], a[j] = a[j], a[i]
            parity ^= 1
    return parity


def check_parent_stereo(parent, state, match):
    """Tetrahedral centres specified in both parent and state must have the same configuration. Compared directly on
    the mapped neighbour order (RDKit chiral tags are relative to the atom's bond order), so CIP-label changes caused
    by the protonation itself cannot trigger a false rejection. Centres whose H count changes are skipped (recorded)."""
    p, s = Chem.RemoveHs(parent), Chem.RemoveHs(state)
    tags = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    checked = []
    for i_s, i_p in enumerate(match):
        a_s, a_p = s.GetAtomWithIdx(i_s), p.GetAtomWithIdx(i_p)
        if a_s.GetChiralTag() not in tags or a_p.GetChiralTag() not in tags:
            continue
        if a_s.GetTotalNumHs() != a_p.GetTotalNumHs():
            continue
        ns = [match[b.GetOtherAtomIdx(i_s)] for b in a_s.GetBonds()]
        np_ = [b.GetOtherAtomIdx(i_p) for b in a_p.GetBonds()]
        same_tag = a_s.GetChiralTag() == a_p.GetChiralTag()
        if same_tag == bool(_perm_parity(ns, np_)):
            raise StateError(f'stereocentre {a_p.GetSymbol()}{i_p} of the parent has the opposite configuration in the state')
        checked.append(i_p)
    return checked


def build_state(state_smiles, parent_smiles=None):
    """Return (Mol with explicit H, info dict). Raises StateError with a precise reason on any incompatibility."""
    state = _parse(state_smiles, 'state')
    parent = _parse(parent_smiles, 'parent') if parent_smiles else Chem.Mol(state)
    match, sites, bonds, ties = site_diff(parent, state)
    dH, dq = check_transformation(sites)
    unassigned = [x for x in stereo_elements(state) if not x['specified']]
    if unassigned:
        raise StateError(f'unassigned stereo elements {unassigned}: specify them, or expand the state with '
                         'expand_stereo (policy "enumerate") so each stereoisomer is its own explicit state')
    kept = check_parent_stereo(parent, state, match)
    target = Chem.MolToSmiles(state, isomericSmiles=True)
    molH = Chem.AddHs(state)
    info = dict(
        state_smiles_input=state_smiles,
        parent_smiles=parent_smiles or state_smiles,
        state_smiles_canonical=target,
        skeleton_key=skeleton_key(state),
        state_class=classify(sites, bonds, dq),
        net_charge=Chem.GetFormalCharge(state),
        parent_net_charge=Chem.GetFormalCharge(parent),
        dH_total=dH,
        changed_sites=sites,
        changed_bonds=bonds,
        heavy_map_state_to_parent=match,
        n_equally_minimal_mappings=ties,
        n_atoms_with_H=molH.GetNumAtoms(),
        n_heavy=state.GetNumAtoms(),
        elements=dict(Counter(a.GetSymbol() for a in molH.GetAtoms())),
        charged_atoms=[(a.GetIdx(), a.GetSymbol(), a.GetFormalCharge()) for a in state.GetAtoms() if a.GetFormalCharge()],
        stereo_elements=stereo_elements(state),
        parent_stereocentres_checked=kept,
    )
    return molH, info


def torsions_parent_keyed(edge_index, edge_mask, heavy_map):
    """TD rotatable bonds as parent heavy-atom index pairs (TD bonds to H are never rotatable: 1-atom fragment)."""
    out = []
    for u, v in edge_index.T[edge_mask].tolist():
        out.append(tuple(sorted((heavy_map[u], heavy_map[v]))))
    return sorted(out)
