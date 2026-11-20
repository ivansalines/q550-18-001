# Q550 kick203 v18.6 — Second-Order Material Source Transport Audit

**Protocol version:** `2026-10-16-q550-v18.6`  
**Status:** executable bundle; **new Q550 second-variation solver must run on the original server**. Existing v18.5 Q550 echo measurements are archived; *source reconstruction numbers are not precomputed*.  
**State:** actual Q550 support ODE, 11D, `Q=550`, `kappa0=(55/6)*(10/550)` and `gamma_rho=.3`.  
**Constraints:** `OFF_CONSTRAINT_DIAGNOSTIC`. No invariant manifold retraction/certification provided; no physical hysteresis conclusion.

## Scientific question

Exactly *when* does the nonlinear (even-in-kick) material response originate, and how is that response transported to a later time? v18.5 identified a quadratically scaling 11D echo after an anti-kick that cancels the first-order tangent response. v18.6 reconstructs its source from the original vector field rather than inferring it from echo samples.

Let `f(y)` be the original autonomous 11D Q550 ODE. At the frozen material source `s0=sA(xi=203)` with state `y0(s0)` and unit kick `v`, define:

```
y0' = f(y0)
J(s) = Df(y0(s))
P' = J P                  P(s0) = I
u(s) = P(s) v
q' = J q + S2(s)           q(s0) = 0
S2(s) = (1/2) D²f(y0(s))[u(s),u(s)]
q(s) = P(s) integral_s0^s P(t)^(-1) S2(t) dt
```

The symmetric history satisfies `1/2 (y_plus+y_minus)-y0 = amp_abs² q + O(amp_abs⁴)` for a smooth flow. The tangent-compensated **even echo immediately after the anti-kick** is exactly this symmetric history at its pulse time, because the compensations of opposite signs cancel. At subsequent read times it is linearly transported to leading order; the material source was accumulated **before** the anti-kick.

Thus `q(s)` is an independently predicted **11D vector**, not a scalar fitted to the previous memory curve.

## Source construction and control gates

1. Validate the SHA256 of the original v18.5 sendback and check its 36 checkpoint entries. Original Q550 ODE, Jacobian, orbit, event map, and original x500 ZIP are validated using the frozen v13 SHA256 gate.
2. Reconstruct original baseline from 400-period orbit, 70-period tail and original dense solve, separately for `standard` (`rtol=1e-12,atol=1e-14`) and `tight` (`rtol=3e-13,atol=3e-15`). Recover `P` from the original analytic 11×11 Jacobian. The original source state must match v18.5 within an explicit relative threshold.
3. Obtain `S2` from **centered finite differences of the original analytic Jacobian**, along the *unit full-state tangent* `w=u/||u||`: `S2=||u||² [J(y+h*w)-J(y-h*w)] w /(4h)`; compare `h=2e-4` and `h/2` (configurable). This is not an interval-certified Hessian.
4. Integrate `P^-1 S2` by composite Simpson with `0.025 xi` fine cells; independently compare Simpson cells of `0.050 xi`, sharing function evaluations. Reconstruct total and six source-age bin contributions. Fail when the cocycle inverse is ill-conditioned (`cond(P)>1e8` by default).
5. Reconstruct `q` at delays `0.25,0.5,1,2,4,8 xi`, and propagate *only previously born source* for echo read ages `0,0.25,0.5,1,2 xi`. Compare each predicted **full 11D vector** to the original v18.5 echoes, at eps `1e-8`, `3e-8`, `1e-7`, for both tolerances; report relative norm, absolute error and angle. Also extrapolate observed `r_even/amp_abs²` to zero amplitude with a diagnostic linear fit versus `amp_abs²`.
6. Diagnose how much the source born in `[0,.25], [.25,.5], [.5,1], [1,2], [2,4], [4,8] xi` contributes to a later `q`. **Source-bin norm ratios do not sum to 100%**, as signs and vector cancellations matter. Report signed projection onto the eventual q, plus the full 11D arrays.

**Important distinctions:** a quadratic material echo can arise in any nonlinear smooth ODE and is not itself hysteresis. The reverse kick is not a reversal of time or a parameter cycle. The source-time decomposition is deterministic finite-amplitude perturbation theory in 11D, still without invariant-manifold constraints.

## Server invocation

Copy the bundle ZIP into `/root/elab1/hysteresis`:

```bash
unzip q550_kick203_second_order_material_source_transport_v18_6_bundle.zip
cd q550_kick203_second_order_material_source_transport_v18_6_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_6_server.log
```

Python defaults to `python3.13` and can be overridden with `PYTHON=...`. Server-specific original Q550 path overrides follow the original v13 defaults, for example:

```bash
python3.13 scripts/run_source_transport_v18_6.py \
  --base-root /root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle \
  --source-bundle-zip /root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip
```

For already completed source checkpoints, **postprocess without another Q550 replay**:

```bash
python3.13 scripts/analyze_source_transport_v18_6.py
python3.13 scripts/pack_sendback_v18_6.py
```

Checkpoints can resume per tolerance (`checkpoint/standard_source_transport.npz`, `tight_source_transport.npz`); add `--force` to regenerate both. `--hessian-step` and `--cell-xi` expose finite-difference and quadrature sensitivity. Changing these controls requires `--force`, otherwise existing checkpoints are reused.

## Output and readout

- `checkpoint/*_source_transport.npz`: `xi_nodes`, `source_11D`, `source_integrand_11D`, fine and half-Hessian-step arrays, complete source-age-bin vectors, source coefficient `q` at each target delay, inverse cocycle condition numbers, original material kick.
- `results/echo_vs_second_variation.csv`: comparison with **all 180 paired v18.5 readings** (6 delays × 5 read ages × 3 amplitudes × 2 tolerances), including 11D vector error and angle.
- `results/source_age_bins.csv`: how contributions of six birth intervals transport to each delay.
- `results/source_profile.csv`: local source magnitude and cocycle conditioning versus xi.
- `results/source_tolerance_comparison.csv`, `results/extrapolated_zero_amplitude.csv`, `results/summary.json`, `results/REPORT.md`, three plots.
- `q550_kick203_second_order_material_source_transport_v18_6_sendback.zip`: all scripts actually used, source notebooks/scripts, original v18.5 archive with historical ancestry, numerical inputs/checkpoints/results, README, manifests and checksums.

## Acceptance criteria

- No source result asserted until real Q550 run succeeds.
- All 36 original echo input checkpoints present and integrity-verified.
- Observed versus independently predicted 11D material response should converge toward agreement at decreasing kick amplitudes, and under tighter solver tolerances / finer quadrature / smaller Hessian step. Deviations are reported, not hidden.
- Componentwise reconstructed age-bin contributions must sum to q to numerical precision.
- A persistent, signed age-bin contribution can be labeled a *transported nonlinear response* only; constitutive hysteresis or time-reversal asymmetry requires a different physical protocol.

Source-code comments, file labels and technical documents are in English. The seven local tests are **synthetic checks of the method**, not Q550 source findings.
