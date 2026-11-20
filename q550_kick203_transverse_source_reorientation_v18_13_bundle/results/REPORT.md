# Q550 v18.13 — Transverse Source Reorientation & Radial Turning Audit

## Reproduction
Original v18.12 sendback independently checked using CRC and SHA256. Data: standard/tight; 560 nodes each; previous genealogical chain retained inside original nested archive.

## Four radial crossings and directional coherence

| event | xi (standard, PCHIP) | local cos ±0.05 | broad cos ±0.25 | broad cos ±0.50 |
|---|---:|---:|---:|---:|
| 1 | 2.625653077 | +0.988751 | +0.721113 | +0.026649 |
| 2 | 3.949609768 | +0.964669 | +0.154173 | -0.989975 |
| 3 | 5.560114949 | +0.991532 | +0.783075 | +0.118767 |
| 4 | 6.958293257 | +0.975762 | +0.400263 | -0.839421 |

- Maximum standard/tight root drift: 2.84e-10 xi.
- Maximum linear/cubic/PCHIP method drift: 7.48e-05 xi.
- Maximum even/odd decimation root drift: 6.74e-06 xi.
- Maximum standard/tight coherence drift: 4.89e-10.
- Continuous tangent directions were compared using spherical Levi-Civita parallel transport in a **declared Euclidean initial 10D gauge**.

## Metric sensitivity
| fixed observer | number of radial turnings | locations (xi) |
|---|---:|---|
| identity | 4 | 2.62565, 3.94961, 5.56011, 6.95829 |
| uniform_scale_1p7 | 4 | 2.62565, 3.94961, 5.56011, 6.95829 |
| mild_diagonal | 4 | 2.63812, 3.97453, 5.57111, 7.03827 |
| strong_diagonal | 4 | 2.67232, 4.13643, 5.55706, 7.41154 |

## Interpretation and hard limits
At the four radial turnings, rho_prime passes through zero while residual line turning remains strictly positive at sampled/interpolated points. The transverse fraction approaches one at each event algebraically, since radial source projection is zero there; this is **not** separate physical evidence.
Local tangent-direction coherence is strong at ±0.05 xi, while ±0.25 and ±0.50 xi resolve progressive reorientation, including sign reversal for some events at the broadest window.
The effect is a property of the co-moving source-frame residual under the chosen Euclidean metric. True coordinate-free claims are limited to rank/nonvanishing of the appropriate wedge and exact tensorial equations. No distinguished material metric, physical torsion, constraint, or hysteresis loop has been established.

## Next decisive experiment
To claim history dependence beyond source accumulation, construct distinct admissible histories sharing the same reduced material state and phase (with verified tolerance) and measure unequal residual response, or derive a material metric from the actual ODE constitutive structure.
