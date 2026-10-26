# family_A summary - joint protonation/tautomer + conformer generation
Scouted 2026-10-26 (UTC). Backends: OpenAlex (search_works, citation_neighbours), general WebSearch, WebFetch.
GOOGLE SCHOLAR WAS NOT SEARCHED (no tool available in this session).

## Headline finding
No paper found generates a protonation (or tautomer) microstate z JOINTLY with 3D geometry x, and none conditions a
conformer generator on pH. Every system that touches both sides is a SEQUENTIAL PIPELINE: enumerate or predict the
state, then generate conformers for each enumerated state. This is an absence in the indexes searched, not proof of
non-existence - OpenAlex title-word search was useless for these concepts (see decisions.md) and ChemRxiv/arXiv were
reached only through general web search.

## Closest work, with what it actually does
1. AcepKa (PlayMolecule, JCIM 2026 / arXiv:2604.00841) - FULL TEXT READ. Closest overall because it is the only
   system whose INTERFACE matches p(z,x|G,h): input is SMILES (or 3D conformers) plus a pH value; output is the
   pH-weighted protonation ensemble plus conformers. Mechanism is the opposite of joint: SMARTS-rule combinatorial
   microstate enumeration (net charge in [-2,+2]) -> Uni-Mol transformer predicts standard Gibbs free energy per
   microstate -> Boltzmann weights at the requested pH -> AceConfgen (RDKit-equivalent distance geometry + MMFF94)
   generates a conformer ensemble FOR EACH microstate. z is never sampled from a learned distribution, and geometry
   never feeds back into state selection. This is the concrete sequential baseline to beat.
2. Uni-pKa (2023) - abstract only. The underlying ensemble formulation: enumerate microstates, learn per-microstate
   free energy, Boltzmann-weight by pH and net charge. Defines the z-space and the thermodynamic target that any
   learned p(z|G,h) must reproduce. Scoring model, not a generator; its internal conformer handling is Unknown.
3. FlexiFlow (arXiv:2511.17249) - FULL TEXT READ. Closest on the generative-mechanism side: flow matching that
   samples a molecule together with MULTIPLE conformations, with a per-atom CHARGE-TYPE categorical sampled by a
   discrete-flow (Campbell) update alongside coordinate flow matching, trained with explicit hydrogens. But the
   categorical is formal charge for valence bookkeeping, not a pH-dependent protonation state, and the paper never
   mentions protonation, tautomers or pH.
4. LoQI / ChEMBL3D (ChemRxiv 2025) - abstract only. Conformer diffusion trained on 250M+ AIMNet2-optimised
   geometries in implicit solvent that reportedly span multiple protonation states, with STEREOCHEMISTRY injected as
   a discrete graph-level condition. The conditioning trick is the right shape for conditioning on z, but protonation
   is supplied in the input graph, never generated, and there is no pH input.
5. pKAce (arXiv:2407.11103) - abstract only. Equivariant NN for protonation-site prediction trained on several
   conformers per molecule, reported to reduce prediction variability across conformations: direct published evidence
   that z-prediction depends on x, i.e. motivation for coupling.
6. Epik v7 (ChemRxiv 2023) - abstract only. The production rule+ML enumerator at a target pH (7.4 +/- 2.0); the
   other half of the sequential baseline.
7. OMol25 (arXiv:2505.08762) - abstract only. Large QM dataset whose benchmark suite includes BOTH conformer and
   protonation optimisation tasks, with protonation reported as the hardest energy-based evaluation. Candidate
   external scorer for generated (z,x) pairs.
8. CORDIS project 101210115 - project record, not a paper. Proposes flow matching for missing protonation states
   together with multiconformer PROTEIN structures. Adjacent competitor on the protein problem; no outputs found.

## Also checked, nothing found
Diffusion docking / pocket-conditioned generation (DiffDock, SurfDock, DynamicBind, SigmaDock) and co-folding models
(AF3, Boltz-2, Chai-1): no result showed any of them sampling ligand protonation or tautomer states; ligand
protonation is a fixed preprocessing choice. Protein CpHMD / pKa ML (DeepKa, PypKa, Titr-DMD) is state prediction,
not generation. No "tautomer-aware" 3D generator surfaced under four distinct phrasings.

## Key unknowns
- Whether an unindexed or very recent ChemRxiv/arXiv preprint already does joint z/x; ChemRxiv's own search API was
  not queried directly in this pass (only via general web search).
- Uni-pKa's internal conformer treatment (does the free-energy model see a conformer ensemble per microstate?).
- Whether ChEMBL3D exposes a protonation-state label that LoQI could be conditioned on.
- Whether any commercial tool (Epik/Schrodinger, OpenEye) couples state and conformer generation internally.
