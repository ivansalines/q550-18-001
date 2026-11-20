# Q550 v18.18 — Source Recurrence Cancellation & Transported Echo Residual Audit

## Scientific conclusion

The exact weighted quadrature of instantaneous projected source pairs reproduces the v18.17 integrated kernel to numerical roundoff. **This is an algebraic identity**: it proves the observable kernel is fully determined by the supplied source on the same frozen trajectory. It neither certifies nor excludes physical hereditary dynamics inside that source.

Nonzero discrepancies of amplitude-blind, endpoint-only and age-scrambled constructions show limitations of those observer baselines, not constitutive memory. The source-only identity is a high-strength falsification control against overinterpreting an apparent residual.

## Diagnostic table

| Gauge | Closure max abs | Independent-g max abs | Amplitude-gated unit residual RMSE | Endpoint holdout RMSE | Endpoint fitted holdout RMSE |
|---|---:|---:|---:|---:|---:|
| standard/euclidean | 2.220e-16 | 1.867e-04 | 0.0129 | 0.0727 | 0.0732 |
| standard/strong_diag | 2.220e-16 | 2.085e-04 | 0.0151 | 0.0740 | 0.0722 |
| tight/euclidean | 2.220e-16 | 1.867e-04 | 0.0129 | 0.0727 | 0.0732 |
| tight/strong_diag | 2.220e-16 | 2.085e-04 | 0.0151 | 0.0740 | 0.0722 |

## Definitions

- Pair direction `c(a)=<q(z-a),q(z-a-L)>_G/(||q(z-a)||_G ||q(z-a-L)||_G)`.
- Exact `K=∫r(a)s(a)c(a) da/sqrt(∫r(a)^2 da ∫s(a)^2 da)`.
- Magnitude gate `H=∫rs da/sqrt(∫r² da ∫s² da)`; directional-only predictor `mean(c)`; amplitude-gated predictor `H*mean(c)`.
- Covariance residual `K-H*mean(c)=H*(weighted_mean(c)-mean(c))` is completely determined by age-varying instantaneous source norms and orientations.
- Endpoint-only prediction is `cos(q(z),q(z-L))`; calibrated endpoint has 2 parameters trained on xi=4..6, tested on xi>6.
- Age scrambled pairs roll the older source across the age grid; not a physical counterfactual.

## Checks

- Authenticated CRC/SHA256 for the original v18.17 server sendback and the nested genealogy through v18.7.
- Recomputed v18.17 kernel from its independent checkpoint and compared against original archived surface values.
- Compared derivative of cumulative A to separately saved original source integrand.
- Standard/tight; alternate diagonal metric; quadrature n=40,80,160; delayed-source age permutation.

## Next falsifiable step

Reconstruct independently perturbed ODE trajectories with matched reduced observables and phase, and explicitly separate linear predicted differences from nonlinear response. Do not claim memory merely from same-trajectory source correlation.
