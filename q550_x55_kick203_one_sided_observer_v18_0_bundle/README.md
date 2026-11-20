# Q550 kick 203 — one-sided observer lens v18.0

Goal: disentangle the v17.9 crossing band caused by a **centered finite
 difference** that straddles the +/- restart, from a surviving small directed
 residual that is robust to the *observer* choice.

The test is **B-frame only**: leave section A, its 11D frame, and cocycle M
unaltered (from original v17.9 server arrays). Reintegrate original +/- Q550
nonlinear histories solely to obtain tangent probes t(s-h),t(s),t(s+h) for B.
Calculate 4 alternatives: centered, forward, backward, and non-straddling
(backward before B crossing, forward after). On all variants reuse the same
8-term off-shell Shapley allocation, no new physical model.

## Run on server

From the extracted bundle directory:

```bash
bash scripts/run_on_server.sh 2>&1 | tee v18_0_server.log
```

Requires the **same server-side base bundle and historical x500 source ZIP** as
v17.9, at their validated paths. Override `--base-root` and
`--source-bundle-zip` in the Python invocation if those folders were moved.

Before accepting any results the script verifies original source SHA and
compares regenerated central B frames and Shapley vectors with the **v17.9
archived arrays** at exactly the same xi, not interpolated points.

Output `q550_kick203_one_sided_observer_v18_0_sendback.zip`, containing scripts,
original source references, audit CSVs, figures, regression and hashes.

## Mathematical comparison

Let `t(s)` be the normalized full11D velocity along one nonlinear history.

- `d_C=(t(s+h)-t(s-h))/(2h)` (centered; historical reference)
- `d_F=(t(s+h)-t(s))/h` (forward)
- `d_R=(t(s)-t(s-h))/h` (backward)

All are Gram-projected with the same current `t(s)` to build the oriented
11x2 local frame. `s=s_B(xi)`, with `s_B(xi_B)=s_kick`.

By geometry alone, the centered estimator's straddling region is
`abs(xi-xi_B)<h/b1`; forward straddles for `xi_B-h/b1<xi<xi_B`;
backward for `xi_B<xi<xi_B+h/b1`. The non-straddling variant excludes
both, but at the exact crossing there is no unique two-sided definition.

**Only B is changed**, not A or M; this prevents confounding a frame change
with a changed tangent-cocycle transport. Test includes hhalf/hquarter,
standard/tight and relative kick amplitudes 1e-8/3e-8.

The generated figures show the distribution of `||B||`, angle(A,B), and
signed final delta against a frozen 11D A axis. The off-shell Shapley
allocations are mathematical diagnostics and do NOT certify a material flip,
a physical impulse, or a causal three-component mechanism.

## Reproducibility and files

- `scripts/run_kick203_one_sided_observer_v18_0.py`: actual generator.
- `scripts/pack_sendback_v18_0.py`: packer.
- `inputs/q550_kick203_event_micro_v17_9_sendback.zip`: original previous
  full11D arrays and exact v17.9 xi grid (not re-created).
- `reference/history_bridge_v16.py`, `reference/q550_state_impulse_response_v13.py`,
  `scripts/frame_geometry.py`: original code used for the algebra and Q550
  trajectories.
- `reference/run_kick203_event_micro_v17_9.py`: prior generator source.
- `tests/test_one_sided_frame.py`: pure geometry smoke tests.
- `SOURCE_MANIFEST.json`, `SHA256SUMS.json`, `requirements.txt`.

## Important scientific boundary

This is NOT a new full11D propagation result yet. Source/syntax/unit tests can
be checked offline, but actual one-sided Q550 samples must be obtained on the
user's server from the original Q550 orbit and derivatives.
