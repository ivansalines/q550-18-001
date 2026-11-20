# Q550 v18.8 — Material Metric & Constraint Sensitivity Audit

Measured from the user-provided v18.7 full11D checkpoints; no new Q550 ODE solve.

## Falsifiable result

- Euclidean (11D state-coordinate) metric retains exactly three sign-crossing roots.
- Continuous positive-definite pullback mixture removes the first pair, then the final root, as alpha increases.
- A constant positive-definite metric obtained by whitening at a different reference xi may yield a different number of roots.
- All these experiments measure observer/metric dependence, **not** a physical material metric.

## Data

| Metric | standard roots (grid-resolved) |
|---|---|
| mixed_pullback_alpha_0 | 3 |
| mixed_pullback_alpha_0.05 | 3 |
| mixed_pullback_alpha_0.1 | 3 |
| mixed_pullback_alpha_0.15 | 3 |
| mixed_pullback_alpha_0.2 | 3 |
| mixed_pullback_alpha_0.25 | 1 |
| mixed_pullback_alpha_0.3 | 1 |
| mixed_pullback_alpha_0.4 | 0 |
| mixed_pullback_alpha_0.5 | 0 |
| mixed_pullback_alpha_0.75 | 0 |
| mixed_pullback_alpha_1 | 0 |
| fixed_reference_whitening_xi_1 | 2 |
| fixed_reference_whitening_xi_2 | 4 |
| fixed_reference_whitening_xi_4 | 3 |
| fixed_reference_whitening_xi_8 | 4 |

## Diagonal scale sensitivity

Deterministic seed 184, 200 draws per maximum weight ratio. No ensemble constitutes a probability distribution over plausible physical metrics.

- Maximum ratio 1.1: 3 crossings=200
- Maximum ratio 1.5: 3 crossings=200
- Maximum ratio 2.0: 3 crossings=200
- Maximum ratio 4.0: 1 crossings=19, 3 crossings=181
- Maximum ratio 10.0: 1 crossings=32, 3 crossings=161, 5 crossings=5, 7 crossings=2

## What is and is not constrained

Constraint input: NOT SUPPLIED (off-constraint diagnostic).
Optional user-supplied rank-stable constraint Jacobians produce tangent/normal leakage and projected-metric diagnostics. A projector onto ker(DC) is NOT a second-order on-manifold construction: the true second variation can have normal curvature prescribed by D²C.

The original v18.7 source data came from unconstrained kicks. Even a genuine DC cannot retroactively make those trajectories constrained.

## Reproduction

Run `bash scripts/run_on_server.sh` from the bundle root (no server-specific Q550 files needed).
For optional constraint arrays see `templates/CONSTRAINT_INPUT.md`. The original v18.7 archive retains scripts, inputs, data, manifests and complete v18.6/v18.5 history.

## Numeric qualifications

Cubic Hermite reconstruction resolves roots bracketed by the 0.0125-xi grid. Closely spaced, unbracketed pairs may be missed. Mixed-metric annihilation thresholds are interpolated diagnostics, not rigorous bifurcation certificates.
