# Q550 kick 203 — v18.1 second-order one-sided B observer

## Target

Determine whether the small 11D residual direction surviving the near-opposition
of Shapley A and B after the kick 203 B-section restart converges under
refinement of the **B observer only**. This is an operational observer test,
not proof of a physical flip.

Reference B crossing: `xi_B = 210.861136603178`.

## Controlled experiment

- Preserve **A section frames from the archived v17.9 hquarter stencil** for
  all comparisons. Do not re-estimate A when B's derivative step changes.
- Preserve all full11D cocycle matrices M from v17.9 without recomputing them.
- Recreate exactly the same nonlinear `+/-` histories and their common restart
  using the historically pinned x500 source archive and dense-sampling bundle.
- Change only the **B frame tangent difference**, formed from the same unit
  tangent samples on the same history.
- Four B-step sizes: `hhalf`, `hquarter`, `heighth`, `hsixteenth`.
- One-sided first order:
  `F1 = [t(s+h)-t(s)]/h`, `R1 = [t(s)-t(s-h)]/h`.
- One-sided second order:
  `F2 = [-3t(s)+4t(s+h)-t(s+2h)]/(2h)`,
  `R2 = [3t(s)-4t(s-h)+t(s-2h)]/(2h)`.
- The **non-straddling** observer uses backward differences on the left of B
  crossing and forward differences on the right. The exact event sample is
  censored for these observers.
- Recalculate the exact symmetric 3-factor Shapley decomposition (`A, B, M`)
  from the observed section frames, and verify `A+B+M=Delta` in full11D.

**Terminology boundary:** The section-A frame is frozen, but the Shapley *A
contribution* may change with B, because the Shapley allocation is formed from
mixed counterfactual descriptors. Do not call Shapley A an invariant vector.

## Convergence diagnostics

The script compares entire 11D `Delta` vectors at successive B-step sizes,
including their norms, angles, and signed component along a reference axis
frozen from v17.9. It saves the actual 11D Shapley arrays. The geometric
error-ratio tests are descriptive, particularly near a restart; numerical
convergence does not certify a material torsion, force, or physical impulse.

The `A/B` and `S/M` angles are **unsigned 11D vector angles**. Signs of scalar
components require their named reference axes.

## Dependencies / provenance

Embedded minimum input archives:

- `inputs/q550_kick203_event_micro_v17_9_sendback.zip` — original 320-point
  micro-lens arrays with frozen BA, BB and M.
- `inputs/q550_kick203_one_sided_observer_v18_0_sendback.zip` — original
  centered versus first-order one-sided result, used for regression.
- `inputs/singular_direction_vectors.csv` — kick direction.
- `reference/*.py`, `scripts/frame_geometry.py` — actual original mathematical
  code of the v18.0 / v17.9 chain, copied unchanged for functional genealogy.

External archived dependencies (unchanged, same as the successful v18.0 run):

- `/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle`
- `/root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip`

The source ZIP SHA256 is checked against the historical pinned source by the
underlying `q550_state_impulse_response_v13` module. All archive x-grid and
centered frame reproductions are checked. For the original hquarter B step,
first-order one-sided residuals are checked against the v18.0 sendback.

## Run on server

From the extracted bundle directory:

```bash
bash scripts/run_on_server.sh 2>&1 | tee v18_1_server.log
```

The runner executes unit tests, integrates the same nonlinear +/- histories,
performs the observer diagnostics, and packs:

`q550_kick203_second_order_observer_v18_1_sendback.zip`

**If the calculations succeed but packing fails**, do not rerun the orbit:

```bash
python3.13 scripts/pack_sendback_v18_1.py \
  --results "$(pwd)/results/kick203_v18_1" \
  --out "$(pwd)/q550_kick203_second_order_observer_v18_1_sendback.zip"
```

## Main output files

- `observer_orders.csv` — B observer order, step, signed scalar components,
  unsigned A/B and S/M angles, norm and cancellation diagnostics.
- `residual_convergence.csv` — full11D Delta agreement across B step sizes.
- `order1_vs_order2.csv` — full11D comparison at the *same* B step size.
- `near_event_convergence_summary.csv` — descriptive convergence metrics.
- `second_order_full11D_arrays.npz` — actual A, B, M, Delta 11D arrays and
  fixed reference axes.
- `regression_checks.csv` — original v17.9 central frame and v18.0 first-order
  regression errors.
- `summary.json`, `REPORT.md`, PNG charts, and source manifests.

## Interpretation boundary

An observer-dependent near-90 degree band is an operational pseudo-flip; the
existence of such a band is **not** a certified physical torsional reversal.
The center-versus-non-straddling result in v18.0 established the need to test
how the **small 11D residual**, not only a large A/B angle, converges.
