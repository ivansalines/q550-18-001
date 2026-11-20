# Q550 — v18.17 Moving Ridge, Instantaneous-Source Recurrence & Window-Support Audit

Original archived Q550 kick203 state, no ODE rerun. All interpretations conditional on declared 10D metric/projection.

## Fixed-support ridge (time >=4, lag 1..2, paired ages 0..1)

| Metric | Median best lag | Lag range | Median K at minimum | Broad minimum span (K <= min+0.02) | Edge hits |
|---|---:|---|---:|---:|---:|
| euclidean | 1.3861 | 1.2635–1.4154 | -0.8809 | 0.175 | 0.0% |
| moderate_diag | 1.3929 | 1.2702–1.4251 | -0.8823 | 0.175 | 0.0% |
| strong_diag | 1.4133 | 1.2890–1.4504 | -0.8865 | 0.175 | 0.0% |

## Sensitivity of the minimum to how recurrences are observed

| Construction | Median lag | Range | Median K |
|---|---:|---|---:|
| matched_1xi | 1.3861 | 1.264–1.415 | -0.8809 |
| recent_half_xi | 1.3651 | 1.283–1.456 | -0.8920 |
| older_half_xi | 1.3669 | 1.283–1.456 | -0.8863 |
| instantaneous_pair | 1.3691 | 1.279–1.493 | -0.9062 |
| variable_support_3xi | 1.3842 | 1.301–1.416 | -0.8685 |

## Independent checks

- Maximum standard vs tight fixed-support K drift: `1.390e-10`.
- Maximum Euclidean vs strong-diagonal ridge location change: `0.04599 xi`.
- Maximum derivative-of-A vs original source-g fixed-support K disagreement: `2.085e-04`.
- Maximum 20-to-80 node K quadrature difference: `1.152e-03`.
- Maximum 40-to-80 node K quadrature difference: `7.660e-04`.

## Interpretation

A negative paired-age kernel may be a recurring pattern of a smooth forcing vector, not a delayed response of material degrees of freedom. The constant paired-age support removes one lag-dependent overlap bias from v18.16; the variable-support and instantaneous comparisons quantify other observer dependencies. No causal propagation law is inferred.

## Audit trail

Original original v18.16 server sendback is retained verbatim under inputs/. Every earlier source script, raw output and checkpoint is nested inside. CRC and SHA256 are checked at runtime.

## Next falsifiable question

Repeat the ridge analysis on independently restarted physical trajectories with matched instantaneous phase and state but controlled, different valid histories; only then attempt a constitutive memory claim.
