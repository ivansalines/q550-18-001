# Q550 v18.6 — Material quadratic source transport

Actual Q550 data. Off-constraint diagnostic. The analytic 11D Jacobian was differenced to obtain the local second directional derivative.

## Global metrics

- Median immediate relative **11D vector** prediction error: 8.1696e-05
- Max immediate relative 11D vector prediction error: 0.0018638
- Median Simpson coarse/fine error: 2.1554e-07
- Median Hessian step/half-step error: 1.329e-08
- Median extrapolated zero-amplitude comparison error: 0.00035818

## Errors by delay and tolerance

| Tolerance | Delay xi | Median immediate vector error |
|:--|--:|--:|
| standard | 0.25 | 0.000225885 |
| standard | 0.5 | 0.000140411 |
| standard | 1 | 6.70013e-05 |
| standard | 2 | 2.90617e-05 |
| standard | 4 | 1.69473e-05 |
| standard | 8 | 2.07329e-05 |
| tight | 0.25 | 0.000164235 |
| tight | 0.5 | 0.000136208 |
| tight | 1 | 0.000110527 |
| tight | 2 | 8.99339e-05 |
| tight | 4 | 7.34581e-05 |
| tight | 8 | 3.96867e-05 |

## Interpretation

- The source is S2(s)=0.5*D²f(y0(s))[P(s)v,P(s)v]. Its transported age-bin contributions form a signed 11D vector sum, not a sum of positive memory percentages.
- Increasing K2 magnitude with delay does not imply monotonically accumulating source; source contributions can rotate and cancel after propagation.
- All material kicks remain off the certified constraint manifold; no independent hysteresis mechanism is inferred.
- Use the CSV/NPZ for directional sign, source-birth bins, and error checks.
