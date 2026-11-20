# Q550 v18.16 — Source-Lag Interference Kernel & Echo Persistence Audit

**Frozen-source observer diagnostic, kick 203.** This bundle advances v18.15 (including the subsequent anchor/lag correction) and uses authenticated source history from the original v18.7 checkpoints embedded in the **unchanged** original v18.15 server sendback. No 400-orbit baseline replay is required. It is not a new ODE integration or an experimental demonstration of physical constitutive hysteresis.

## Reproducible server commands (exact ZIP root checked)

Copy `q550_kick203_source_lag_interference_v18_16_bundle.zip` to the hysteresis directory, then:

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_source_lag_interference_v18_16_bundle.zip
cd q550_kick203_source_lag_interference_v18_16_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_16_server.log
```

It creates `/root/elab1/hysteresis/q550_kick203_source_lag_interference_v18_16_sendback.zip` with **one top-level directory, exactly matching the bundle root**. A single root is enforced by `pack_sendback.py`. For a repeat run, extract to a clean directory rather than merging old results. Python 3.10+ with numpy/scipy/matplotlib.

## Question

In v18.15 the single most antagonistic individual source-slab pair at events 2 and 4 suggested a separation near 1.5 xi. Does an **aggregated** negative source-lag correlation persist at 1.5 throughout the path, or is it an event-conditioned maximum from a few bins? Is the value robust to (a) resolution, (b) source interpolation, (c) the chosen positive-definite coordinate metric, and (d) a fixed versus growing historical observation window?

## Mathematical definition

Given the co-moving accumulated source A(xi) and the original anchor e=A(1), fix a positive-definite diagonal metric G. Project with `Q_e^G v = v - e (e^T G v)/(e^T G e)`. For an observation z and age edges a_i=i Δ:

```
q_i(z) = Q_e^G [A(z-a_i) - A(z-a_(i+1))]
K(z,kΔ) = sum_i <q_i,q_(i+k)>_G /
            sqrt( sum_i ||q_i||_G² * sum_i ||q_(i+k)||_G² )
```

Here each summation uses only valid matched pairs (i, i+k) and source birth `u>=1`. There is a minimum paired support of 0.5 xi. A final partial age bin smaller than Δ is **omitted**, producing a controlled O(Δ) truncation error; the four resolutions test stability. We test post-anchor full history versus a trailing age window of length at most 3.0 xi. Both share the same 0.5 xi minimum paired support. The support shrinks when the lag increases and must be read alongside `K`. This normalized signed quantity is a correlation of source contributions, not a material response kernel deduced from experiments.

## Validation and key findings

- The 7 synthetic/algebraic tests cover positive/negative sources, boundedness, support rejection, additive slab closure, the metric-orthogonal projection, and covariance under *consistent* diagonal coordinate transformations.
- Both numerical tolerance inputs are used; their maximum K difference in the event audit is about `1.75e-10`.
- A separate source-integrand quadrature builds slabs from the archived g rather than differences of interpolated A. Differences in K are saved in `independent_source_check.csv`.
- At source-event 2 (xi ~3.950), K at lag 1.50 xi is about -0.696 (Euclidean, Δ=0.0625); at source-event 4 (xi ~6.958), about -0.802.
- In all 81 xi-profile samples from xi=4 through 8, the lag-1.5 kernel is negative with both full post-anchor and fixed 3.0 xi windows. **This is a property of the chosen gauge and fixed trajectory, not a verified physical propagation delay.**
- The broader negative region includes lag 1.25 xi. At event 3, the most negative lag in [1,2] is 1.25 in the Euclidean metric and 1.5 in the strongest tested diagonal metric. Thus a unique gauge-independent lag cannot be claimed.
- Event 1 lacks the stipulated 0.5 xi overlap for lag 1.5 after the anchor; its K(1.5) must be reported as unavailable, **not zero**.

## Outputs

- `results/summary.json` — structured findings, tolerances, checks, provenance, limitations.
- `results/REPORT.md` — compact results report.
- `results/event_lag_kernel.csv` — event × lag × resolution × metric × history (includes overlap).
- `results/event_lag_minima.csv` — minimum lag in 1–2 xi for each configuration.
- `results/kernel_xi_lag_profile.csv` — signed correlation along xi in 0.05 steps, for two history windows.
- `results/lag_persistence_stats.csv` — negative-run diagnostics, event/nonevent comparisons.
- `results/independent_source_check.csv` — independent g-integral versus A-difference source construction.
- `results/source_lag_kernel_heatmap.png` and `results/echo_persistence_comparison.png` — generated graphics.
- `tests/` and `scripts/` — all actual executable generators and tests.
- `inputs/v18_15_original_sendback.zip` — unchanged prior server sendback. Full original earlier source/input genealogy is nested inside.
- `reference/analyze_anchor_lag_v18_15.py` and `reference/ANCHOR_LAG_REVIEW.md` — preserved pre-anchor correction logic and its provenance; original v18.15 was not rewritten.
- `SHA256SUMS_OUTPUT.json` — digests of all substantive files, checked before sendback packing.

## Interpretation guardrail

A time-lag correlation of sources along a **single** computed trajectory, regardless of sign persistence, cannot distinguish intrinsic mechanical memory from a history-dependent observer construction. The lag `xi` is the experiment's co-moving parameter, not a calibrated laboratory time. In particular, a negative kernel is not negative energy, an anti-causal response, or proof of a torsional flip. Physically identifiable memory would require multiple feasible histories with matched final state and phase, and justified material constraints and metric.
