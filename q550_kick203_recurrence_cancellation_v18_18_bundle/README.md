# Q550 v18.18 — Source Recurrence Cancellation & Transported Echo Residual Audit

**Original Q550 kick 203, frozen-source observer diagnostic.** Uses the authenticated *original server* v18.17 sendback, with the nested v18.16–v18.7 genealogy. It does not rerun 400 orbits or independently integrate the ODE. No constitutive memory, causal delay, or physical flip is inferred from a same-trajectory source correlation.

## Server commands — verified single ZIP root

Place `q550_kick203_recurrence_cancellation_v18_18_bundle.zip` under `/root/elab1/hysteresis` and execute:

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_recurrence_cancellation_v18_18_bundle.zip
cd q550_kick203_recurrence_cancellation_v18_18_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_18_server.log
```

The sendback will appear at:

`/root/elab1/hysteresis/q550_kick203_recurrence_cancellation_v18_18_sendback.zip`

The script always resolves its root from its own location. Only immutable input, code and completed result files are included in SHA256; active `*.log` and the transient `.runtime/` directory are **excluded**, preventing the historical `tee`/hash race.

## Scientific test and falsifiable distinction

The exact v18.17 source-lag kernel, already defined as a finite-age source correlation with fixed paired support, is:

\[
K(z,L)=\frac{\int_0^1\langle q(z-a),q(z-a-L)\rangle_G\,da}
{\left[\int_0^1\|q(z-a)\|_G^2\,da\int_0^1\|q(z-a-L)\|_G^2\,da\right]^{1/2}}.
\]

Here `q=Q_e^G dA/dxi` is the projected 10D co-moving source; `e=A(1)` is the declared reference and `G` the declared metric. The v18.17 kernel is reproduced from authenticated checkpoints and independently checked against archived CSVs.

Define pairwise instantaneous cosine `c(a)` and norms `r(a),s(a)`. Then algebraically

`K = sum_a r(a) s(a) c(a) / sqrt(sum_a r(a)^2 sum_a s(a)^2)`.

The **complete source-weighted closure must equal K** at floating-point precision. A nonzero closure error would indicate a numerical or implementation inconsistency, *not material memory*. The full closure cannot independently falsify existence of constitutive memory, because the input source may already contain the past influence.

Compare deliberately incomplete observer models:

1. `instantaneous`: endpoint cosine `cos(q(z),q(z-L))`, with no age integration.
2. `unit_direction`: unweighted age average of pairwise cosines.
3. `amplitude_gated_unit`: unweighted cosine mean multiplied by the norm-overlap gate.
4. `source_weighted_closure`: full amplitude and age-dependent cosine quadrature; must close exactly.
5. `scrambled`: rotate the old source's age samples by 1/4 window; diagnostic null (not a realizable trajectory).
6. `endpoint_calibrated`: a two-parameter endpoint affine prediction trained only over `xi=4..6`, tested on `xi>6` (a dependent temporal holdout, not independent experiment).

The `amplitude_gated_unit` residual is a source-defined covariance between pair amplitudes and pair orientation. **Its nonzero value is not proof of hereditary dynamics.** Also compare original archived source integrand and dA/dxi; standard/tight tolerances, Euclidean/strong diagonal metrics, and quadrature resolutions 40/80/160.

## Generated contents

- `results/summary.json` — all numeric diagnostics, provenance and limitations.
- `results/REPORT.md` — scientific interpretation and controls.
- `results/source_kernel_decomposition.csv` — full surfaces for closure, incomplete source models, scrambled null.
- `results/time_residual_and_ridge.csv` — lag ridge and residual norm versus xi; temporal holdout marked.
- `results/event_probe.csv` — available radial-event probes (insufficient xi is explicitly marked).
- `results/source_recurrence_residual_maps.png` — comparison of the source-only kernels and observer residual.
- `results/source_ridge_baseline_comparison.png` — crests of the full, endpoint and unit-pair kernels.
- `scripts/` — all current Python implementations, archived loader hierarchy and sendback packager.
- `tests/` — synthetic tests of closure, weighting, metric covariance, age registration, failure detection.
- `reference/` — unmodified mathematical source and prior v18.17 analysis scripts.
- `inputs/v18_17_original_sendback.zip` — unmodified authenticated server sendback with deeper genealogy.
- `SHA256SUMS_OUTPUT.json` — immutable manifest generated after outputs are closed.

## Interpretive checkpoint

The single-trajectory source kernel cannot by itself distinguish an autonomous material history variable from the nonlinear state-dependent forcing that produced the source. The next discriminating experiment would compare independently integrated perturbed ODE trajectories at controlled phase and reduced observables, with a predicted linear response and clearly defined limitations on state matching. For a finite-dimensional deterministic ODE, equal full states at the same phase necessarily have identical futures; such a pair cannot establish different hidden histories without an extended state description.
