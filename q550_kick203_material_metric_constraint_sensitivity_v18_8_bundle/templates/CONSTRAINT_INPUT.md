# Optional constraint Jacobian contract (not included in measured Q550 data)

An optional external analyst may supply **genuine, independently documented**
Jacobian evaluations of one or more Q550 material constraints on the same
unperturbed reference trajectory and at the same 641 xi nodes as v18.7.

`--constraint-npz /path/constraint_jacobians.npz` (NumPy archive with **both**):

- `xi_nodes`: float64 shape `(641,)`, original `0..8 xi` coordinates.
- `constraint_jacobians`: float64 shape `(641, m, 11)`, every row of `DC(y0(xi))`
  in the **same original state coordinates and units** as the v18.7 checkpoints.

`--constraint-provenance /path/constraint_provenance.json` with all mandatory
non-empty string keys: `constraint_name`, `constraint_formula`,
`coordinate_mapping`, `jacobian_construction`, `source_or_reference`.

Example *schema only* (NO Q550 invariants provided):

```json
{
  "constraint_name": "SPECIFY_VERIFIED_CONSTRAINT",
  "constraint_formula": "SPECIFY_FORMULA_C(y)=constant",
  "coordinate_mapping": "SPECIFY_EXACT_y0_through_y10_MEANINGS_AND_UNITS",
  "jacobian_construction": "SPECIFY_ANALYTIC_OR_VALIDATED_DERIVATIVE",
  "source_or_reference": "SPECIFY_VERIFIABLE_ORIGINAL_SOURCE"
}
```

Run: `python3 scripts/analyze_metric_v18_8.py --constraint-npz ... --constraint-provenance ...`

For every node the script computes orthogonal tangent projector `Pi = ker DC`,
its rank, normal leakage of the existing off-constraint source vectors, and
sign crossings in metrics `W_eta = Pi + eta * (I-Pi)` with eta=0,0.01,0.1,1.
At eta=0 this is only a **semidefinite projection**. The script labels all
results **EXTERNAL_UNVERIFIED_CONSTRAINT_JACOBIANS**; it does not validate
constraint conservation or perform a true constrained perturbation.

**Mathematical warning:** a constrained first variation `u` satisfies
`DC(y0)u=0`, but the second variation `w` typically satisfies
`DC(y0)w + D²C(y0)[u,u]=0`, so its **normal** part can be physically required.
Projecting the v18.6 quadratic second variation into ker DC would omit this
curvature contribution. To establish an on-constraint physical result, construct
retracted ± perturbations, check conservation under integration, and compare
with the covariant second variation of the verified constraint manifold.
