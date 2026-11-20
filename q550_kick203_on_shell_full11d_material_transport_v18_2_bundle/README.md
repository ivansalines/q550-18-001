# Q550 kick203 — v18.2 On-Shell Full11D Material Transport Audit

**Scientific distinction:** The archived v18.1 result tests an 11D **two-section dlocal observer**. The new v18.2 experiment tests an 11D **support-model state response** at the *same physical support time* in three nonlinear ODE histories. They have dimension 11 but are NOT interchangeable observables. Raw 11D vectors from the two instruments are never subtracted or dotted together. Only the xi timing and separate scalar magnitudes can be compared as context.

## What the audit actually tests

Keep the kick 203 input direction, normalized amplitude convention, nonlinear Q550 ODE, orbit initialization, and original A/B support maps unchanged. The B event remains `xi_B=210.8611366031777`; the physical impulse occurs at `s0=sA(203)`. Read 47 xi samples from v18.1 and explicitly include the exact event (already present in that archive). Keep tolerances `standard/tight` and amplitudes `1e-8/3e-8`.

At identical physical support parameter `s` for the + and - histories measure

`u(s) = [y_+(s)-y_-(s)] / (2*epsilon*max(1,|y0|))`.

On `s<s0`, the histories are identical by causal construction and `u(s)=0`. At `s=s0`, record the instantaneous source jump separately, never interpret it as a smooth material rotation. Evaluate `u(s_A(xi))` and `u(s_B(xi))`; these are the *same history at two DIFFERENT physical times*, not a simultaneous A/B state difference.

Check independently a true 11D **baseline variational transport** `Phi(s,s0)` from `dy/ds=f(y)` and `J=∂f/∂y`. After the restart, compare `u(s)` to `Phi(s,s0)*direction`. When B is post-restart, test homogeneous B→A prediction `Phi_A * Phi_B^{-1} * u_B` with a matrix-conditioning gate. When B is pre-restart, **no homogeneous B→A prediction is legitimate**: the kick is an intervening source and A must instead be compared to its predicted source response `Phi_A * direction`.

Use a mobile 11D Euclidean reference built directly from the **baseline** flow tangent `t=f/|f|` and its analytic curvature `t'=(I-ttᵀ)Jf/|f|`, independent of B stencil `h`. Construct canonical minimum ambient rotations between successive **rank-1 tangent** and **rank-2 tangent-curvature** subspaces, accumulate `R(s)` from the kick, and independently refine physical integration-grid spacing `0.01 → 0.005`.

In the body-locked observer convention `U(s)=R(s)ᵀ u(s)`, the increment has an **exact algebraic split**

```
u(s_next)-u(s_prev)
  = R(s_prev) [U(s_next)-U(s_prev)]          # co-moving (intrinsic-to-gauge)
  + [R(s_next)-R(s_prev)] U(s_next)          # baseline frame advection
```

The two components are reported in full11D and their vector closure is validated. They need **not** have additive norms. Do not call this a rigid physical torsion decomposition: the Euclidean rank-1/rank-2 framing is an operational gauge, and no SO(11) mechanical body symmetry has been established. The split deliberately excludes the instantaneous kick source.

## On-shell honesty (important)

The `+/-` perturbations are exact numerical **ODE histories**, but the inherited kick direction need not lie on a conserved-invariant manifold. Unless an independent `invariants(y)` evaluator is supplied, summary reports `FLOW_ON_SHELL_ONLY_WITHOUT_CONSTRAINT_CERTIFICATE` / `NOT_CERTIFIED_OFF_CONSTRAINT` and no material physical flip conclusion is permitted. This is the same scientific boundary as v13–v18.1.

For an independently derived invariant callback, use `--invariant-module /path/to/invariants.py`, containing `invariants(y)` that returns a 1D NumPy-compatible vector of model invariants. In that optional branch, the kick is projected tangent to invariant levels and retracted, and invariant drift is measured at all output support points. **This sampled numerical check is not a rigorous interval certificate** and modifies the kick relative to the frozen v18.1 off-constraint direction; report the two branches separately. No invariant is guessed or silently invented.

## Run on Q550 server

The bundle includes the **actual archived** v18.1 sendback (which itself embeds v17.9/v18.0 archives), original singular kick direction, functional historical scripts, audit code, synthetic controls, SHA256 manifests, and output packer. It still needs the same **external archives used successfully in v18.1**:

- `/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle`
- `/root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip`

From the unpacked v18.2 bundle folder:

```bash
bash scripts/run_on_server.sh 2>&1 | tee v18_2_server.log
```

Alternative to override defaults (the wrapper forwards extra CLI options):

```bash
bash scripts/run_on_server.sh --transport-step 0.005 \
  --eps 1e-8 3e-8 --tolerance-sets standard tight
```

After completed calculations the wrapper saves `q550_kick203_on_shell_full11d_material_transport_v18_2_sendback.zip`. If **packing alone** fails, do **not** reintegrate:

```bash
python3.13 scripts/pack_sendback_v18_2.py
```

For embedded archives without ODE integration:

```bash
python3.13 scripts/preflight_v18_2.py
python3.13 -m unittest discover -s tests -v
```

## Deliverables

- `results/kick203_v18_2/material_responses.csv`: same-s 11D material sensitivity at A and B, linear response, phase-tangent and transverse components, body-locked signed projections, rank-1/rank-2 control, original v18.1 scalar norm context.
- `results/kick203_v18_2/moving_frame_increment_decomposition.csv`: exact intrinsic/convective 11D increment norms, rank dependence and algebraic closure; no artificial smoothing or event-straddling.
- `results/kick203_v18_2/full11D_material_arrays.npz`: original full 11D vectors, 11D variational predictor, even-in-epsilon remainder, exact xi/A/B/s0 maps and source vector.
- `results/kick203_v18_2/summary.json`, `REPORT.md` and two PNGs; output source hash manifest packaged with scripts, reference code and archived input zips.

## How to interpret

An enduring co-moving deviation across changes of rank, physical-grid refinement, tolerances, and eps would favor **a material support-state directional response** over a pure section-observer pseudo-flip. It would *not*, by itself, prove a physical torsional flip, conserved constraint, or a dynamical bifurcation. Conversely, a large apparent turn absorbed by reference transport, or an abrupt change only in the dlocal observation, favors an observer contribution.

All source comments and document content are in English; numerical results are generated on the original server, not fabricated locally.
