# family_D decisions (2026-10-26, UTC)

- Scope narrowed on purpose to the mechanism question: HOW are the discrete and continuous parts coupled, and is
  any discrete variable a protonation/charge state. Generic 3D-generation benchmark papers were not noted.
- Representative-per-mechanism rule applied instead of exhaustive coverage, to stay under the ~25-note cap:
  * Gaussian-relaxed one-hots + coordinates in one diffusion  -> EDM (GeoLDM/DiffSBDD inherit this; noted together).
  * True categorical (D3PM) noise + Gaussian coordinates, per-component schedules -> MiDi (primary, full text).
  * Same with x0-parameterisation and large-scale pretraining -> EQGAT-diff.
  * Single transformer over 2D+3D with self-conditioning -> JODO (MUDiff grouped with MolDiff note).
  * CTMC discrete flow + ODE coordinates, shared denoiser -> SemlaFlow (full text) and FlowMol3.
  * Discrete modality diffused, continuous modality a frozen predictor readout -> ProteinGenerator.
  * Torus-valued continuous + simplex categorical -> PepFlow.
  * Theory -> Campbell et al. (DFM/Multiflow), Gat et al., D3PM.
- GeoLDM and DiffSBDD were NOT separately noted: both are EDM-lineage (latent / pocket-conditioned variants) and no
  source read here showed a different discrete-continuous coupling. Their charge handling is recorded as Unknown
  rather than asserted.
- ESM3 dropped: searches did not return usable mechanism detail in this pass and it is sequence-token-only at the
  discrete level; ProteinGenerator and Multiflow already cover protein co-design coupling.
- Variable atom count: not resolved for any model read here. Recorded as a key unknown (it matters because adding or
  removing a proton changes the atom count under explicit hydrogens).
- No paper in this family mentions pH, protonation state or tautomers. Verified by targeted full-text reads of
  MiDi, SemlaFlow and FlexiFlow, not just by abstract absence.
