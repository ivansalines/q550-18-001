# Q550 kick203 — v18.4 Five-Amplitude On-Flow Interaction Scaling Audit

## Exact objective

Test whether the **v18.3 geometric-descriptor mixed-factor residual**
`I_ABC(eps)` is approximately cubic and whether `I_ABC / eps^3` stabilizes
across **five** impulse amplitudes. Keep the real nonlinear support ODE,
initial kick direction, **47 xi samples**, original `sA`, `sB`, full11D
cocycle `M`, A `hquarter` section frames, frozen A projection axis and solver
configurations. Historical levels `1e-8` and `3e-8` are replay-free pinned
anchors, consumed unchanged from the embedded v18.3 archive; compute only
three missing amplitudes on the server.

**Amplitudes:** `5e-9`, `1e-8`, `1.5e-8`, `2e-8`, `3e-8`.

**Two solver tolerances:** standard (`rtol=1e-12, atol=1e-14`) and tight
(`rtol=3e-13, atol=3e-15`). For each missing amplitude and tolerance, integrate
**two** genuinely perturbed ODE histories (+epsilon, -epsilon). Six new
pairs = **12** nonlinear ODE branches, plus one common baseline per tolerance.
The archived 400-period orbit is loaded and its 70-period tail rebuilt only once; completed `(tol,eps)`
checkpoints are reused by default, including after a server interruption.

## Numerical definition

Let `F=(BA,BB,M)` be the three factors in the two-section descriptor
`d(F)` (an 11-component geometric **descriptor**, *not* the 11D state vector).
For each amplitude use actual `F-`/`F+` from the perturbed nonlinear ODE and
write `Fbar=(F-+F+)/2`:

```
S_i(eps) = d(Fbar with Fi=Fi+) - d(Fbar with Fi=Fi-)
Delta(eps) = d(F+) - d(F-)
I_ABC(eps) = Delta(eps) - S_A(eps) - S_B(eps) - S_M(eps)
```

This is an **exact algebraic identity** in the geometric-descriptor space.
The factor hybrid configurations are computational substitutions, not claimed
physical histories. With smoothness and `F+-F-` odd at leading order, the
mixed interaction starts at cubic order, but a numerical cubic exponent is
**not** a proof of real mechanical hysteresis or irreversible memory.

For each `xi` and tolerance, fit both:

- `Delta(eps) = C1*(eps/1e-8) + C3*(eps/1e-8)^3 + C5*(eps/1e-8)^5` in fixed 11D;
- `I_ABC(eps) = K3*(eps/1e-8)^3 + K5*(eps/1e-8)^5` and cubic-only control;
- independent `h16` non-straddling second-order B observer vs centered `hquarter`;
- `I_ABC(eps)/(eps/1e-8)^3` across five amplitudes (frozen units, orientation
  and post-B physical-support sampling).

Report leave-one-amplitude-out errors, local two-level effective exponents,
relative cubic-coefficient spread, fifth/cubic term ratio, frozen signed
projection and **per-amplitude standard/tight solver discrepancies**. If noise contaminates the
smallest amplitude, report an **unresolved numerical floor** rather than
manufacturing an asymptotic law.

## Physical and mathematical distinctions

1. **ODE flow is not an invariant-manifold certificate.** Original kick vectors
   may be off-constraint; there is no new `invariants(y)` certificate in v18.4.
2. The support material state 11D `y(s)` and geometric descriptor 11D `d(F)`
   are separate observables. They are **not** subtracted or compared componentwise.
3. A and B samples correspond to different physical support parameters, not
   simultaneous measurements. Causal masks use `sB-s0`, not just `xi`.
4. The same **single frozen projection axis** from v18.3 `standard eps=1e-8`
   is used for *every* epsilon/tolerance in v18.4.
5. A forward/reverse path or nested-loop measurement is still necessary to
   establish history dependence; cubic response alone does not imply hysteresis.
6. Near an event, the **centered B stencil** spans the restart; only 19 clean
   post-B and 17 clean pre-B samples are used for exponent summaries.
7. Values are kept in original observable units; odd-polynomial coefficients
   use dimensionless `t=eps/1e-8` to improve numerical conditioning.

## Input genealogy and server dependencies

Contains the complete original `v18.3` ZIP (which embeds v18.2, v18.1, v17.9),
original singular-direction CSV, active Python scripts and the historical
source implementations. The archived ZIP is **SHA256 pinned** in
`SOURCE_MANIFEST.json` and every historical source dependency is SHA256
validated by the original v13 script at runtime.

Original server paths (same as v18.1/v18.2):

```
/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle
/root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip
```

The active scan uses the original `support_rhs_float`, analytic
`float_support_jacobian`, DOP853 dense nonlinear history integrator, and the
historical **locally rebased full11D variational cocycle** reconstruction.

## Run on the Q550 server

From `/root/elab1/hysteresis` after copying the ZIP:

```bash
unzip q550_kick203_five_amplitude_scaling_v18_4_bundle.zip
cd q550_kick203_five_amplitude_scaling_v18_4_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_4_server.log
```

The wrapper runs tests, preflight, 12 new nonlinear histories with recoverable
checkpoints, five-amplitude analysis, and packs:

`q550_kick203_five_amplitude_scaling_v18_4_sendback.zip`

If the ODE stage completed but plotting or packaging failed, **do not replay**:

```bash
python3.13 scripts/analyze_scan_v18_4.py
python3.13 scripts/pack_sendback_v18_4.py
```

`python3.13 scripts/preflight_v18_4.py` requires no Q550 server inputs.
Pass `--base-root` and `--source-bundle-zip` to the wrapper to change upstream
Q550 paths. `--force` regenerates all six new checkpoints.

## Expected numerical output (server not yet run)

- `checkpoint/{standard,tight}_eps{5e-9,1.5e-8,2e-8}.npz` (six full11D pairs)
- `results/five_amplitude_samples.csv`: 470 rows = 47 xi x 5 eps x 2 tolerances
- `results/samplewise_odd_fits.csv`: 94 xi/tolerance fits, with leave-one-out
- `results/tolerance_discrepancy.csv`: standard/tight noise-floor diagnostic
- `results/summary.json` and `results/REPORT.md`
- `results/{interaction_scaled,interaction_raw,signed_cancellation,tolerance_sensitivity}.png`
- `run_provenance.json`, tests, console logs, output SHA256 manifest, and
  self-contained sendback ZIP.

No uncomputed amplitudes are represented as measured results. Software and
source documentation is in English for executable provenance.
