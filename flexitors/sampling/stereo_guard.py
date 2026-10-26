"""Exp 3 inference-time constraint operators and per-step stereo verification (EXP-3 §3.1-§3.3, §5, §9.2;
Q-P1-006 / Q-P2-007 answers; round5/DECISION.md fix 3).

Arms (EXP-3 §5, §9.2, §9.4):
* ``O``: TD unchanged (counters only).
* ``K1``: freezes the SMILES-specified structural stereo bonds only (§9.2 "K1: freezes specified structural bonds
  only. Unspecified bonds stay free"). Seed check: if |δ_b| > π/2 the bond is rotated by π before sampling (§3.2).
* ``K1_wrong``: the same bonds frozen to the opposite isomer (φ_b + π): negative control for "fewer degrees of
  freedom" (§5).
* ``K1_canonical``: freezes every TD-rotatable double and triple bond (VonMisesNet rule, §5; E-CHEM-043); specified
  bonds get the seed check, the others keep the seed state.
* ``K2_inf``: TD's score, then δ folded into I = [−π/2, π/2] after every step (§3.2, E-CHEM-048) for specified bonds
  from the start (the fold of the uniform seed randomisation is the K2 prior, §9.2); unspecified structural and
  non-stereo double bonds get the late planarity fold: free until the first reverse step with σ ≤ 0.1π, then folded
  around the nearest of {0, π} at that step (§3.1 revision, §9.2; the 0.1π switch is a pre-declared design choice,
  Q-P1-006 answer).
Freezing = the bond is excluded from the initial torsion randomisation AND its updates are masked (DECISION fix 3).
The operators only produce torsion-update vectors; positions are always changed by TD's own ``modify_conformer`` /
``perturb_batch``, so every arm uses TD's rotation code. With the update mask the score model sees the unmodified
graph; the equivalence with removing the bond from ``edge_mask`` (EXP-3 §3.2 wording) for the other bonds' scores is
tested on gnode118 (plan_P2 V-E3-3).

Verification (§3.3): after every step, per conformer, E/Z violations |δ_b| > π/2 against a defined target (SMILES
φ_b, its K1-wrong opposite, or the late-fold φ*) and planarity violations |δ_b| > 20° in the last 3 steps, through the
array-level ``flexitors.eval.validity.ez_state_array``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem

from flexitors.eval.validity import PLANARITY_DEG, ez_state_array
from flexitors.geometry.internal import dihedral, wrap_angle
from flexitors.geometry.torsion_types import is_structural_stereo, td_reference_quadruple, td_rotatable_bonds

# modify_conformer rotates the v side about (pos[u] - pos[v]) (utils/torsion.py:68, "positive rotation if pointing
# inwards"), which changes the atan2 dihedral (c, u, v, d) by MINUS the update (measured; tested against TD on gnode118)
UPDATE_SIGN = -1.0
LATE_FOLD_SIGMA = 0.1 * np.pi   # EXP-3 §3.1 revision / §9.2
LAST_STEPS = 3                  # EXP-3 §3.3 planarity window
ARMS = ('O', 'K1', 'K1_wrong', 'K1_canonical', 'K2_inf')


@dataclass(frozen=True)
class BondSpec:
    """One TD-rotatable double or triple bond: its slot in the molecule's torsion-update vector (masked edges are in
    bond order), TD's quadruple, class and SMILES target (from the committed table; None if unspecified)."""

    slot: int
    bond_idx: int
    quad: tuple[int, int, int, int]
    kind: str          # 'specified' | 'unspecified' | 'nonstereo' | 'triple'
    phi: float | None


def fold_half_circle(delta: np.ndarray) -> np.ndarray:
    """EXP-3 §3.2 K2: wrap to (−π, π], then x > π/2 ↦ π − x and x < −π/2 ↦ −π − x, so the result lies in I."""
    x = wrap_angle(delta)
    return np.where(x > np.pi / 2, np.pi - x, np.where(x < -np.pi / 2, -np.pi - x, x))


def constrainable_bonds(mh: Chem.Mol, targets: dict[int, float]) -> tuple[list[BondSpec], int]:
    """TD-rotatable double/triple bonds of ``mh`` (TD's featurised molecule) and the number of rotatable bonds.

    ``targets``: {bond_idx: φ_b} of the molecule's specified rows of ``round5/data/ez_targets.csv`` (EXP-3 §9.3)."""
    rot = td_rotatable_bonds(mh)
    out = []
    for slot, k in enumerate(rot):
        b = mh.GetBondWithIdx(k)
        if b.GetBondType() == Chem.BondType.TRIPLE:
            kind = 'triple'
        elif b.GetBondType() != Chem.BondType.DOUBLE:
            continue
        elif is_structural_stereo(b):
            kind = 'specified' if k in targets else 'unspecified'
        else:
            kind = 'nonstereo'
        out.append(BondSpec(slot, k, td_reference_quadruple(mh, k), kind, targets.get(k)))
    return out, len(rot)


@dataclass
class StereoHooks:
    """The arm's operator for one molecule. Arrays are per conformer: positions (b, n, 3), updates (b, m)."""

    arm: str
    bonds: list[BondSpec]
    n_slots: int
    steps: int
    late_phi: dict[int, np.ndarray] = field(default_factory=dict)   # conformer -> φ* per late-fold bond
    counts: dict[str, np.ndarray] = field(init=False)

    def __post_init__(self) -> None:
        if self.arm not in ARMS:
            raise ValueError(f'unknown arm {self.arm}')
        a, B = self.arm, self.bonds
        spec = [b for b in B if b.kind == 'specified']
        self.frozen = {'K1': spec, 'K1_wrong': spec, 'K1_canonical': B}.get(a, [])
        self.freeze_target = {b.slot: (b.phi if a != 'K1_wrong' else wrap_angle(b.phi + np.pi)) if b.phi is not None
                              else np.nan for b in self.frozen}
        self.fold = spec if a == 'K2_inf' else []
        self.late = [b for b in B if b.kind in ('unspecified', 'nonstereo')] if a == 'K2_inf' else []
        self.monitored = [b for b in B if b.kind != 'triple']
        self.counts = {'ez_violations': np.zeros(self.steps, int), 'planarity_violations': np.zeros(self.steps, int),
                       'fold_corrections': np.zeros(self.steps, int)}

    # --------------------------------------------------------------------------------- seeds (before the loop)
    def seed_updates(self, pos: np.ndarray, uniform: np.ndarray) -> np.ndarray:
        """Torsion updates for the seed stage: TD's uniform randomisation (``perturb_seeds``, sampling.py:105-108)
        with frozen bonds excluded, plus the K1 seed check (rotate by π if |δ_b| > π/2 to the freeze target)."""
        upd = np.array(uniform, dtype=np.float64)
        if not self.frozen:
            return upd
        slots = np.array([b.slot for b in self.frozen])
        upd[:, slots] = 0.0
        tgt = np.array([self.freeze_target[s] for s in slots])
        tau = dihedral(pos, np.array([b.quad for b in self.frozen]))
        with np.errstate(invalid='ignore'):
            flip = ~np.isnan(tgt) & (np.cos(tau - tgt) < 0)
        upd[:, slots] = np.where(flip, np.pi, 0.0)
        return upd

    # --------------------------------------------------------------------------------- reverse steps
    def mask_update(self, upd: np.ndarray) -> np.ndarray:
        """Frozen bonds get no update (DECISION fix 3)."""
        if self.frozen:
            upd = np.array(upd, dtype=np.float64)
            upd[:, [b.slot for b in self.frozen]] = 0.0
        return upd

    def fold_updates(self, conf_idx: list[int], pos: np.ndarray, sigma: float | None, step: int | None) -> np.ndarray:
        """Correction updates that fold every fold bond into I around its target after a step (or after the seed
        randomisation when ``sigma`` is None). Late-fold bonds join at the first step with σ ≤ 0.1π, where φ* is the
        nearest of {0, π} to the conformer's current state, kept for the remaining steps. Zero for other arms."""
        upd = np.zeros((len(conf_idx), self.n_slots))
        late_on = bool(self.late) and sigma is not None and sigma <= LATE_FOLD_SIGMA
        bonds = self.fold + (self.late if late_on else [])
        if not bonds:
            return upd
        tau = dihedral(pos, np.array([b.quad for b in bonds]))
        tgt = np.empty_like(tau)
        tgt[:, :len(self.fold)] = [b.phi for b in self.fold]
        if late_on:
            nl = len(self.late)
            for j, i in enumerate(conf_idx):
                if i not in self.late_phi:
                    self.late_phi[i] = np.where(np.cos(tau[j, -nl:]) > 0, 0.0, np.pi)
                tgt[j, -nl:] = self.late_phi[i]
        delta = wrap_angle(tau - tgt)
        change = fold_half_circle(delta) - delta
        if step is not None:
            self.counts['fold_corrections'][step] += int(np.sum(np.abs(change) > 1e-12))
        upd[:, [b.slot for b in bonds]] = change / UPDATE_SIGN
        return upd

    def verify(self, conf_idx: list[int], pos: np.ndarray, step: int) -> None:
        """EXP-3 §3.3 counters after a step: E/Z violations against defined targets; planarity in the last 3 steps."""
        if not self.monitored:
            return
        tau = dihedral(pos, np.array([b.quad for b in self.monitored]))
        base = [self.freeze_target.get(b.slot, np.nan if b.phi is None else b.phi) for b in self.monitored]
        phi = np.tile(np.array(base, dtype=np.float64), (len(conf_idx), 1))
        late_cols = [j for j, b in enumerate(self.monitored) if b in self.late]
        for r, i in enumerate(conf_idx):
            if i in self.late_phi:
                phi[r, late_cols] = self.late_phi[i]
        match, dev = ez_state_array(tau, phi)
        self.counts['ez_violations'][step] += int((~match).sum())
        if step >= self.steps - LAST_STEPS:
            self.counts['planarity_violations'][step] += int((dev > PLANARITY_DEG).sum())
