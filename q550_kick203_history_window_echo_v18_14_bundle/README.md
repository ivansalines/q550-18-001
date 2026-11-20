# Q550 v18.14 — History-Window Ablation & Transverse Source Echo Audit

This is a **reproducible analysis of the real frozen Q550 v18.13 checkpoint genealogy**, not an independent ODE re-simulation or a physical hysteresis certificate. The original v18.13 server sendback is embedded verbatim as `inputs/v18_13_original_sendback.zip`; its complete nested v18.12→v18.11→v18.10→v18.9→v18.7→v18.6 source genealogy is SHA256/CRC checked before use. The scripts used here are included under `scripts/`, with mathematical tests, a self-contained input, outputs, plots and SHA256 manifest.

## Single command on the user's server

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_history_window_echo_v18_14_bundle.zip
cd q550_kick203_history_window_echo_v18_14_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_14_server.log
```

The generated sendback will be found at `/root/elab1/hysteresis/q550_kick203_history_window_echo_v18_14_sendback.zip`. No earlier server files or 400-orbit regeneration needed. Python 3.10+ with numpy/scipy/matplotlib.

## What is actually ablated

Let `A(xi)` denote the historical source integral in the **initial co-moving 10D frame** and `g(xi)` the corresponding source integrand, with `s = ds_ODE/dxi`, `A(0)=0`, `A'=s*g`. The same real underlying ODE trajectory is retained. At every evaluation xi for a retained recent-history horizon `h` (xi units), compute

- `A_recent(xi;h) = A(xi) - A(xi-h)`
- `A_old(xi;h) = A(xi-h)`
- **Moving-boundary derivative** `A_recent' = s*(g(xi)-g(xi-h))`, `A_old' = s*g(xi-h)`.

All history vectors are already pulled back into the same initial frame, so subtracting them is well-defined. `e=A(1)` is the fixed early-memory anchor of the preceding audits. With `w=(I-eeᵀ/||e||²)A`, measure the unsigned angle between residual lines, their projector speeds, and the portion of old history *outside the plane formed by e and recent history*. At the four original full-history radial turnings, compute rank([e,A_recent,A_old]) and direction changes on local/broad windows. The rank and wedge non-vanishing are invariant under invertible coordinate changes; **angle, speed, magnitude and event position are not** without a physically derived metric.

## Caution

*This removes old terms only from the observer's integral*, **not from the actual material ODE**. It measures which portions of source history are needed to reconstruct the **cumulative observable**, NOT a proof of material memory or hysteresis. Source history is automatically additive; a second physical trajectory at matched state and phase is still needed for a stronger physical claim.

## Files

- `scripts/checkpoint_io_v18_14.py`: recursive SHA256/CRC checkpoint provenance, using original reader scripts.
- `scripts/history_window_math.py`: fixed/recent/old history decomposition and exact boundary terms.
- `scripts/analyze_v18_14.py`: standard/tight event comparisons, profiles, rank-3 tests and numerical stability.
- `scripts/plot_v18_14.py`: graphs derived exclusively from source checkpoints.
- `scripts/pack_sendback.py`: hashed sendback archive with original genealogy.
- `tests/test_echo.py`: synthetic positive/negative control tests.
- `results/summary.json`, `results/REPORT.md`, CSVs and PNGs: measured outputs.
- `requirements.txt`, `ENVIRONMENT.json`, `SHA256SUMS_OUTPUT.json`: reproducibility metadata.
