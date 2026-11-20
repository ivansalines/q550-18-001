# Q550 v18.10 — Material plane transport and constraint preflight

All data sourced from the v18.9 sendback and its hash-verified nested v18.7/v18.6 checkpoints.

The geometric objects are B=E wedge K2 and b=a wedge A; source-plane departure is a wedge A wedge g.

## Findings (numerical, not interval certified)

- **standard**: rank-2 counts {'2': 560}; rank-3 counts {'3': 560} (SVD cutoff 1e-8 in original source coordinates).
  B'=transport + birth finite-difference median residual 1.621e-04.
  Source-frame median departure indicator 0.901835 (NOT invariant as magnitude).
  Reparametrized GL test: 7/7 matching numerical ranks; tensor rate discrepancy 1.045e-15.
- **tight**: rank-2 counts {'2': 560}; rank-3 counts {'3': 560} (SVD cutoff 1e-8 in original source coordinates).
  B'=transport + birth finite-difference median residual 1.621e-04.
  Source-frame median departure indicator 0.901835 (NOT invariant as magnitude).
  Reparametrized GL test: 7/7 matching numerical ranks; tensor rate discrepancy 1.718e-15.

## Interpretation

The rank of [a,A] and [a,A,g], and the fact of wedge vanishing, are invariant under exact invertible coordinate changes. Singular values, dot products, angles, normalized volumes, and the numerical thresholds used to estimate rank are NOT invariant. The zero at xi=1 is constructed (A(1)=a).
A nonzero a wedge A wedge g means the incoming source escapes the existing source-memory 2-plane (in the linearized 11D state). It does not prove material body torsion or constitutive hysteresis.
The exterior derivative has two algebraic parts: homogeneous cocycle transport and birth. Their norms and cancellation depend on metric; the equation transforms covariantly under fixed invertible coordinate changes.
No source-derived mechanical energy, nonlinear material constraint, or admissible second-order shape manifold has yet been certified.
Run the optional original-ODE diagnostic to inventory field structure; an ODE RHS alone does not imply physical metric uniqueness.
