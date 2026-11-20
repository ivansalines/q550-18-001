# Q550 v18.17 — Moving Antagonism Ridge & Instantaneous-Source Recurrence Audit

**Kick 203; authenticated original source; fixed-trajectory observer audit.** This is a direct continuation of v18.16. It is not a reintegration of the ODE, a constitutive hysteresis certificate, or a physical time-delay measurement.

## Commands — exact ZIP hierarchy

The archive has **one** directory root `q550_kick203_moving_antagonism_ridge_v18_17_bundle/`. Place the archive in `/root/elab1/hysteresis` and run:

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_moving_antagonism_ridge_v18_17_bundle.zip
cd q550_kick203_moving_antagonism_ridge_v18_17_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_17_server.log
```

It produces `/root/elab1/hysteresis/q550_kick203_moving_antagonism_ridge_v18_17_sendback.zip`. Send this ZIP back for scientific review.

**Log fix:** `v18_17_server.log` is a live session transient and is intentionally excluded from both the hash manifest and packaged sendback. SHA256 covers only immutable code, input, and completed reports. This avoids the `tee`/manifest race found in v18.16. The `SHA256SUMS_OUTPUT.json` file is regenerated after computation, then verified before writing the ZIP. For repeated runs, use a clean extracted directory.

## Scientific question

The v18.16 pooled source-slab correlation displayed a broad negative lag band with a moving minimum around 1.25–1.5 xi. Is the moving ridge retained when **paired ages have equal support for every lag**, or is it induced by shrinking pair overlap as lag grows? How does its location compare with direct instantaneous-source recurrence, short-age subsets, arbitrary diagonal metric changes and different source reconstructions?

## Source and gauge

The input is the **unchanged v18.16 server sendback**, with its complete embedded v18.15 → v18.7 authentic source genealogy. The input archive's CRC and 30-file SHA256 manifest are independently checked at runtime. The original co-moving cumulative source is `A(xi)` in the phase-fixed 10D sector, measured in a *declared* coordinate metric. For each tested positive-definite diagonal metric `G`, project orthogonal to `e=A(1)`:

`q(t)=Q_e^G[dA(t)/dξ]`.

For each observation `z in [4,8]` and lag `L in [1,2]`, the primary matched-support kernel is:

```text
C(z,L) = integral[a=0..1] <q(z-a),q(z-a-L)>_G da
         / sqrt( integral[0..1] ||q(z-a)||_G² da
               * integral[0..1] ||q(z-a-L)||_G² da )
```

The source is post-anchor for every compared time (`z-a-L >= 1`). It has **fixed 1 xi paired-age support** for every lag, unlike the variable-support kernel of v18.16. The integrals are evaluated with fixed midpoint quadrature counts `N = 20,40,80`, and the source interpolator is PCHIP. The minimum of `C(z,L)` over the scanned grid (step 0.025 xi) is optionally refined by a local quadratic fit. All endpoints, finite support, and broad minima are recorded. Quadratic interpolation does not imply a physical timing precision of 0.001 xi.

## Independent checks

- Original archived `source_integrand_11D`, multiplied by `slope`, is independently interpolated; its fixed-support correlations are compared with the PCHIP derivative of cumulative `A`.
- `standard` and `tight` profiles must agree.
- Midpoint 20/40/80 source quadratures are compared; do not infer rigorous convergence beyond the observed precision.
- Alternate tests: source-age halves `[0,.5]` and `[.5,1]`; instantaneous source cosine `q(z)` versus `q(z-L)`; variable-support 3 xi pooled correlation (overlap shrinks with lag).
- Three metric diagonals (Euclidean, moderate, strong) and three reference directions (`e=A(.75), A(1), A(1.25)`) are compared. The source start restriction remains `u>=1` for all anchor-reference experiments.
- Nine synthetic mathematical tests cover signs, support, metric covariance, projection and ridge localization.

## Expected preflight finding

In the original Euclidean declared gauge, the moving minimum has median ~1.386 xi, observed range ~1.264–1.415 xi across 81 points, and median kernel about -0.881. This is a **broad**, moving opposition trough, not a sharp phase-locked delay. Its location moves in other metrics and under other comparisons. Endpoint-source recurrence provides a similar but distinct pattern; these checks cannot decide whether the opposition is propagated memory or regenerated forcing on a single trajectory.

A previous radial event near xi=2.626 has **insufficient post-anchor history** for the fixed-paired-age `L=1..2` scan, so it is explicitly labeled unavailable in the event comparison. The event near xi=3.950 is represented by the boundary sample xi=4.0 and labeled a proxy with its offset; it is not presented as an exact event measurement.

## Outputs and script genealogy

- `results/summary.json` and `results/REPORT.md` — full methods, findings, limitations, provenance.
- `results/moving_ridge_surface.csv` — xi × lag × method × metric × tolerance correlations.
- `results/ridge_tracking.csv` — sampled and quadratic-refined lag minima, broadness and scan boundary flags.
- `results/four_event_ridge_comparison.csv` — event-anchored crosschecks (with missing-support handling).
- `results/quadrature_and_independent_source_checks.csv` — convergence and independent integrand verification.
- `results/source_integrand_quality.csv` and `results/anchor_sensitivity.csv` — source and gauge checks.
- `results/moving_ridge_heatmap.png` — moving ridge over the signed fixed-support kernel.
- `results/ridge_methods_comparison.png` — comparison of source recurrence definitions.
- `results/fixed_support_echo_persistence.png` — selected lag profiles.
- `scripts/` — *actual* generating scripts, historical input loaders, packager and hash script.
- `tests/` — nine synthetic checks.
- `reference/` — historical source interpretation and mathematical loader code.
- `inputs/v18_16_original_sendback.zip` — unmodified original server sendback, preserving all deeper genealogy.
- `SHA256SUMS_OUTPUT.json` — immutable-file manifest.

## Epistemic limits

The strongest negative lag, sign, ridge trajectory, and cosine are observer- and metric-dependent. No unique intrinsic lag or physical torsional flip is established; nor does any source autocorrelation from a single trajectory establish constitutive memory. One must compare independently realized trajectories at matched instantaneous state and phase, with mechanically justified constraints and metric, to distinguish a transported internal memory from a recurrent forcing pattern.
