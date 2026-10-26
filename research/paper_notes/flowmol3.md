# FlowMol / FlowMol3
## Citation
Dunn I., Koes D.R., "Mixed Continuous and Categorical Flow Matching for 3D De Novo Molecule Generation" (FlowMol),
arXiv:2404.19739, 2024; and "FlowMol3: Flow Matching for 3D De Novo Small-Molecule Generation", arXiv:2508.12629, 2025.
Code: github.com/dunni3/FlowMol. PDF: research/papers/arxiv_2508.12629.pdf. Checked 2026-10-26.
## Fulltext status
abstract_only (abstracts + secondary summaries).
## Method (as reported)
FlowMol3 generates atom types, FORMAL CHARGES and bond orders by simulating continuous-time Markov chains (discrete
flow matching) while coordinates follow an ODE; modality-specific losses are summed with weights. FlowMol (v1)
compared simplex-based categorical flows (SimplexFlow) against plain continuous flows on one-hot vectors and
reported the simpler non-categorical treatment performed better - a directly relevant negative result for how to
parameterise a discrete z.
## Relevance to p(z,x|G,h)
Same coupling pattern as SemlaFlow but with charges inside the CTMC. The FlowMol v1 ablation (simplex vs naive
continuous relaxation of categoricals) is evidence to cite when choosing the parameterisation of z.
## Unknowns
Whether charges are sampled or auxiliary in FlowMol3 exactly; resolvable from arXiv:2508.12629 full text.
