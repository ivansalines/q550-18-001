# Q550 v18.16 — Source-Lag Interference Kernel & Echo Persistence Audit

## Scope

This audit evaluates all xi checkpoints (sampled 0.05 xi) rather than only four selected extrema. It compares the full post-anchor history u>=1 against a fixed recent 3 xi window, then varies slab width and diagonal metric. Frozen original Q550 kick203 v18.7 source inputs remain unchanged.

## Event summary (Euclidean metric, source from xi>=1, slab width 0.0625)

| Event | xi | K(lag=1.5) | Most negative lag in [1,2] | K at that lag |
|---|---:|---:|---:|---:|
| 1 | 2.625653 | n/a | 1.0000 | -0.0889 |
| 2 | 3.949610 | -0.6957 | 1.5000 | -0.6957 |
| 3 | 5.560115 | -0.7871 | 1.2500 | -0.8009 |
| 4 | 6.958293 | -0.8020 | 1.5000 | -0.8020 |

## Persistence (xi>=4, 0.05 xi sample spacing)

| Window | Lag xi | Negative fraction | Max negative run xi | Median K |
|---|---:|---:|---:|---:|
| post_anchor_full | 1.25 | 1.000 | 4.05 | -0.7949 |
| post_anchor_full | 1.50 | 1.000 | 4.05 | -0.7778 |
| post_anchor_full | 1.75 | 1.000 | 4.05 | -0.5353 |
| fixed_3xi | 1.25 | 1.000 | 4.05 | -0.8214 |
| fixed_3xi | 1.50 | 1.000 | 4.05 | -0.8399 |
| fixed_3xi | 1.75 | 1.000 | 4.05 | -0.5960 |

## Checks

Largest standard/tight difference of kernel = 1.750e-10.
Maximum standard/tight event drift = 2.841e-10 xi.
Largest independent source-integrand A-versus-g kernel disagreement = 1.261e-06.
Closure checks are in summary.json. Resolution and metric sensitivity are explicit in event_lag_kernel.csv; changes in strongest lag are not a failure but evidence against metric-independent localization.

## Interpretation

The signed K is a normalized sum of dot products of age-offset source slabs in a declared 10D metric. A negative value describes opposing vectors across a time separation, not anti-causality, work, physical reversal or hysteresis.
The previous ~1.5 xi value was selected from the single most opposing slab pair at two preselected extrema; this kernel averages matched pairs over their available overlap and can refute that value as a broadly persistent lag.
For genuine constitutive memory, compare physically realizable trajectories with matched state and phase but differing histories; this audit does not perform that test.

## Genealogy

Original authentic v18.15 server sendback SHA256: `ad19c1ab78314d6c87130624d58e8f2c63b747c5099c777d9633bede6feb33cd`.
Nested original v18.14 -> v18.7/6 data and their source scripts are inside inputs/v18_15_original_sendback.zip. The v18.15 anchor/lag addendum script is retained separately under reference/.
All code and report text are in English. No baseline 400 orbit replay.
