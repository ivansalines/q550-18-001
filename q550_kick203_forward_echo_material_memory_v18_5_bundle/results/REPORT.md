# Q550 kick203 v18.5 forward echo — measured audit

Checkpoints: 36; sampled readouts: 180.
Lag xi-equivalent: [0.25, 0.5, 1.0, 2.0, 4.0, 8.0].
Relative kick amplitudes: [1e-08, 3e-08, 1e-07].
Range of mean echo fractions at reversal: [2.0356152114060137e-05, 0.0011372314355344389].
Median even scaling exponent (descriptive only): 2.00013.
Median odd scaling exponent (often close to numerical floor): 2.56422.

## Immediate nonlinear memory left after tangent-matched anti-kick

| Delay xi | Median echo / uncancelled write |
|---:|---:|
| 0.25 | 6.10723e-05 |
| 0.5 | 0.000107908 |
| 1 | 0.000180427 |
| 2 | 0.000315908 |
| 4 | 0.000341157 |
| 8 | 0.00028276 |

## Scientific limit

This is an off-constraint nonlinear **echo / residual state** diagnostic; the reverse operation is a forward-time anti-kick, not backward-time ODE integration nor a thermodynamic/control-parameter hysteresis loop. Control for apparent memory: compare full 11D matched states and solver tolerances.

