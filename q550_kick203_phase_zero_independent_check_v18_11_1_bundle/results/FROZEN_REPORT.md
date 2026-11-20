# Q550 v18.11.1 — Independent phase-zero countercheck (frozen part)

## Status

- Frozen J→P propagation: **PASSED** for standard and tight.
- Tests with intentionally corrupted structural zeros: **3/3 PASSED**.
- Independent nonlinear initial-condition perturbation: **NOT RUN LOCALLY**; original ODE is on the Q550 server.

## Quantitative results

| Tolerance | RK4 1 substep | RK4 2 substeps | RK4 4 substeps | Error-reduction ratios |
|---|---:|---:|---:|---:|
| standard | 3.318716065e-08 | 2.070640153e-09 | 1.295092915e-10 | 16.02749, 15.98835 |
| tight | 3.318717149e-08 | 2.070652035e-09 | 1.295106698e-10 | 16.02740, 15.98827 |

Checks: original Jacobian last row exactly zero across 641 samples; original cocycle lower-left block exactly zero and bottom-right corner exactly one across the same samples; upper-right phase-to-shape coupling is nonzero. RK4 propagation is a new numerical implementation using a cubic-interpolated J, with archived P used *only* for final comparison.

These checks establish structural consistency, **not** a physical metric, nonlinear material constraint, rigorous interval theorem, or autonomous 10D flow.
