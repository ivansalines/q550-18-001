# Q550 v18.15 — Anchor-Sensitive Interference and Age-Lag Review Addendum

This is a **post-hoc diagnostic** of the original v18.15 server sendback. It is not a modified original result, a new ODE integration, or a proof of constitutive material hysteresis. All source references and scripts/comments are in English.

## Reproduction

The archive contains a single directory root. In that root run:

```bash
python3 scripts/analyze_anchor_lag.py \
  --input inputs/v18_15_original_sendback.zip \
  --output results/anchor_and_lag.csv
```

Python 3.10+, numpy and scipy (versions compatible with v18.15). The script validates the v18.15 manifest and every nested archived manifest through v18.7. It independently reads original source A(xi) from the standard/tight v18.7 checkpoints. No 400-orbit baseline replay.

## Diagnostic motivation

The transverse projection is `Q_e = I - ee^T/(e.e)` where `e=A(1)`, and `A(0)=0`. Therefore `Q_e[A(1)-A(0)] = 0` algebraically. Breaking source history [0,1] into age slabs nevertheless creates sub-vectors that can mutually cancel even though their complete contribution is zero. To distinguish that decomposition effect, compare the original full-history slab partition with a second partition that **ends at original source xi=1**. Both have the identical resultant `Q_e A(z)`.

## Quantitative findings

At age width 0.0625 xi (standard tolerance):

| Event | Full longitudinal cancellation | Post-anchor cancellation | Full vector cancellation | Post-anchor vector cancellation | Strongest opposing slab-center age separation |
|---|---:|---:|---:|---:|---:|
| 1 | 0.053123 | 0.000000 | 0.310939 | 0.243881 | 1.1250 |
| 2 | 0.068940 | 0.000000 | 0.752993 | 0.741076 | 1.5000 |
| 3 | 0.435110 | 0.426461 | 0.656795 | 0.649763 | 1.2500 |
| 4 | 0.346949 | 0.340401 | 0.857603 | 0.855705 | 1.4375 |

Pairs at events 2 and 4 retain an approximate 1.44-1.50 xi separation across the fine resolutions inspected. This is a candidate lag in the signed *Euclidean source-slab dot-product kernel*, not a measured mechanical propagation delay, invariant time, or new constitutive law. The stronger post-anchor vector cancellations at events 2 and 4 **survive** removal of the zero-net pre-anchor source segment.

## Limits

Every norm, angle, dot product, cancellation percentage, and selected strongest pair depends on the chosen Euclidean co-moving 10D representation and on binning/interpolation. The rank argument of previous audits has a different invariance class. An exact zero-net pre-anchor projected contribution does not imply the early physical history is dynamically irrelevant: the projection itself uses `A(1)` and the earlier state determines the later trajectory. The source history is diagnosed on a fixed trajectory and is not physically ablated.

## Genealogy

`inputs/v18_15_original_sendback.zip` contains the unmodified server sendback, including its nested v18.14 -> v18.7 source archives and scripts. `scripts/analyze_anchor_lag.py` generates `results/anchor_and_lag.csv` with standard/tight and four resolutions, while the manifest confirms file integrity. The current addendum is an auxiliary review artifact, not a replacement of v18.15.
