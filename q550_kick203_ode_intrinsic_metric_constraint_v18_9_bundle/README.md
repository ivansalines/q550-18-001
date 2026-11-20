# Q550 kick203 — v18.9 ODE-Intrinsic Metric & Constraint Audit

**Status:** mandatory frozen-Q550 11D tangent/source audit EXECUTED on October 18, 2026; optional exact-original-ODE replay **NOT YET EXECUTED**; genuine physical/material metric and constraint manifold **NOT IDENTIFIED**. Report no unwarranted physical claims.

## Motivation and primary result

Original Q550 v18.7: E(xi) and K2(xi) show three (grid-resolved) Euclidean orthogonality crossings at xi ~= 4.135198957, 4.997336555, 7.174435814. The v18.8 continuous SPD metric audit established that these roots can disappear under a change of metric. The v18.9 present audit reproduces the three original crossings from actual Q550 full11D checkpoints in both integration tolerances and compares them against the **homogeneous ODE linearization's inverse-cocycle transported metric** `W(s)=P(s)^-T P(s)^-1`. No source-frame crossings are observed on xi=1..8, and source-frame cosine remains positive. This metric is **not** asserted to be a physical metric: its initial Euclidean choice remains arbitrary.

## Quantities and units

Q550 original state dimension 11; `s=s0+slope*xi`, with slope read from the v18.6 source (not guessed). Let

```
P' = J(s) P   in support time s
E=P A_early
K2=P A_total
A_total'=P^-1 S2,  in support time s
W=P^-T P^-1

E^T W K2 = A_early^T A_total
(E^T K2)'_xi = slope * E^T(J^T+J)K2 + slope * E^T P g
```

- State-coordinate `I` metric: a coordinate-dependent diagnostic.
- `W` metric: time-dependent inverse transport of **the same arbitrary initial** state-unit choice, also diagnostic.
- A hypothetical material-energy metric `G(y)`: only a documented physical structure can justify this label.
- `DC(y)`: true constraint tangent only if an actual function `C(y)` and its invariance have been checked; a Jacobian sample by itself cannot certify it.

## What is checked now

1. Original v18.7 and nested v18.6 data member hashes, 641 nodes, 11D source+transport identities.
2. Reproduction of original Euclidean 11D root locations, with standard/tight numerical comparison.
3. ODE linearized cocycle `dP/dxi=(ds/dxi) J P` finite-grid regression; strain vs fresh source derivative decomposition, source-frame pullback identity.
4. SPD metric unit-coordinate covariance test: transforming **the metric** along with the state preserves the dot, while changing units and silently resetting the metric to `I` does not preserve root locations.
5. Linear necessary screens for an affine first integral (`c^T J=0`) and a constant quadratic tangent metric (`J^T G+GJ=0`). **Neither establishes an invariant, physical SPD metric, or excludes nonlinear constraints.** The Jacobian contains an exact zero row at component 10 in original 0-based indexing, producing a rank-deficient semidefinite quadratic candidate, but it does not determine whether the corresponding state component itself is conserved or advances autonomously.
6. The two 11D transported vectors at each Euclidean orthogonality root are linearly independent (measured sine near 1), so their mutual sign cannot be required to be the same for every SPD metric.

## Run: self-contained archived-data audit (no original server dependencies)

```bash
unzip q550_kick203_ode_intrinsic_metric_constraint_v18_9_bundle.zip
cd q550_kick203_ode_intrinsic_metric_constraint_v18_9_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_9_server.log
```

Defaults to `python3.13`, falls back to `python3`, respects `PYTHON`. Requires numpy/scipy (`requirements.txt`). BLAS threads fixed to 1. Produces:

- `results/summary.json`
- `results/observable_geometry.csv`
- `results/REPORT.md`
- self-contained `../q550_kick203_ode_intrinsic_metric_constraint_v18_9_sendback.zip`

## Optional original ODE replay (no fresh 400-orbit integration)

Once original Q550 source modules are accessible under the pinned server base directory:

```bash
Q550_REPLAY=1 bash scripts/run_on_server.sh 2>&1 | tee v18_9_ode_replay.log
```

This invokes `scripts/model_structure_check_v18_9.py`: checks exact historical hashes for `audit/support/tube` modules, integrates **directly from the archived v18.6 source state at kick 203** to xi=8, compares the 11D Jacobian to the historical checkpoint, and writes `results/model_replay_summary.json`. No full baseline orbit or 400-period computation is needed. If the source is missing, this optional mode fails rather than inventing an RHS; run the default mode for the fully independent frozen-data audit.

For physically grounded metric/constraint candidate, see `templates/STRUCTURE_API.md` and pass both `--structure-module` and `--structure-provenance` to the optional original-ODE replay script. Scripts and manuscript/report comments are entirely in English.

## Interpretation boundaries

This audit is not a rigorous interval proof of root count, a proof of true conserved quantities, an intrinsic Riemannian geometry, a rigid-body flip certificate, or an on-manifold second-order perturbation. Specifically `DC K2+0.5 D2C[u,u]` can contain a physically necessary normal component; naively projecting the 11D second-order response away from constraints is unjustified. Initial state source kicks were off-constraint.

**Next real gate:** reproduce original f(y) and check actual proposed `C(y)`, `DC`, `D2C`, and `G(y)` with physically traceable provenance, after confirming original coordinate definitions. Only then test on-manifold perturbations using retraction and independent trajectory integration.
