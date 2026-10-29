"""State enumeration for the sequential baseline: Uni-pKa's template enumerator (ext/unipka_orig/enumerator, Apache-2.0)
driven as in its enum_ensemble(), with one documented change and explicit validation.

- Templates: simple_smarts_pattern.tsv (the one Uni-pKa's enum_ensemble uses). Proton add/remove at template sites only;
  Uni-pKa's tautomer hook is a no-op ("TODO: Tautomerism enumeration"), so same-charge protomers appear only when
  reachable by deprotonation + reprotonation at template sites. C-H tautomers (keto/enol) are never generated.
- Charge window [lower, upper] relative to nothing (absolute net charge). Uni-pKa's enum_ensemble loops
  range(q0 - 2, lower, -1) / range(q0 + 2, upper), which never reaches `lower` / `upper` themselves (e.g. no -2 state
  from a neutral parent with the default -2..2). Here the window is inclusive; faithful_bounds=True reproduces Uni-pKa.
- Every enumerated SMILES is re-validated with jmcg.data.canonical.microstate (single fragment, no radical) and must
  carry the net charge of its level; violations are returned as failures, not dropped silently.
"""
import hashlib
import importlib.util
import os

from rdkit import Chem

from jmcg.data import canonical as cn

ENUM_DIR = os.environ.get('UNIPKA_ENUM_DIR', '/scratch/nishanth.r/jmcg/ext/unipka_orig/enumerator')
TEMPLATE = 'simple_smarts_pattern.tsv'
_MOD = None


def _enum_module():
    global _MOD
    if _MOD is None:
        spec = importlib.util.spec_from_file_location('unipka_enumerator', os.path.join(ENUM_DIR, 'main.py'))
        _MOD = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_MOD)
    return _MOD


def provenance(template=TEMPLATE):
    p = os.path.join(ENUM_DIR, template)
    sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()  # noqa: E731
    return dict(enumerator=os.path.join(ENUM_DIR, 'main.py'), enumerator_sha256=sha(os.path.join(ENUM_DIR, 'main.py')),
                template=p, template_sha256=sha(p))


def enumerate_ensemble(smiles, lower=-2, upper=2, maxiter=10, template=TEMPLATE, faithful_bounds=False):
    """Returns (states, failures). states: list of dict(smiles, state_id, charge, level) incl. the input state."""
    E = _enum_module()
    a2b, b2a = E.read_template(os.path.join(ENUM_DIR, template))
    q0 = Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))
    ens = {}
    s0 = [smiles]
    end_lo, end_hi = (lower, upper) if faithful_bounds else (lower - 1, upper + 1)
    if q0 > lower:
        s0, b1 = E.enumerate_template(s0, a2b, b2a, maxiter=maxiter, mode='A')
        if b1:
            ens[q0 - 1] = b1
        for q in range(q0 - 2, end_lo, -1):
            if q + 1 in ens:
                _, b = E.enumerate_template(ens[q + 1], a2b, b2a, maxiter=maxiter, mode='A')
                if b:
                    ens[q] = b
    if q0 < upper:
        a1, s0 = E.enumerate_template(s0, a2b, b2a, maxiter=maxiter, mode='B')
        if a1:
            ens[q0 + 1] = a1
        for q in range(q0 + 2, end_hi):
            if q - 1 in ens:
                a, _ = E.enumerate_template(ens[q - 1], a2b, b2a, maxiter=maxiter, mode='B')
                if a:
                    ens[q] = a
    ens[q0] = s0
    states, failures, seen = [], [], set()
    for level, smis in sorted(ens.items()):
        for s in smis:
            try:
                can, mol = cn.microstate(s)
            except cn.IdentityError as e:
                failures.append(dict(smiles=s, level=level, reason=str(e)))
                continue
            q = Chem.GetFormalCharge(mol)
            if q != level:
                failures.append(dict(smiles=s, level=level, reason=f'net charge {q} != level {level}'))
                continue
            sid = cn.state_id(can)
            if sid in seen:
                continue
            seen.add(sid)
            states.append(dict(smiles=can, state_id=sid, charge=q, level=level, is_input=can == cn.microstate(smiles)[0]))
    if not any(s['is_input'] for s in states):
        failures.append(dict(smiles=smiles, level=q0, reason='input state not in its own ensemble (filtered by Uni-pKa)'))
    return states, failures
