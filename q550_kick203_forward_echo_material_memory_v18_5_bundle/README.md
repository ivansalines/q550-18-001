# Q550 kick203 v18.5 — Forward-Echo Full11D Material Memory Audit

**Version:** `2026-10-16-v18.5`  
**Status:** Prepared/locally tested; **not yet executed on the actual Q550 server**.  
**Scope:** Deterministic Q550 support ODE in 11D, historical kick direction rank 1, fixed parameters Q=550 and `kappa0=(55/6)*(10/550)`. Original Q550 code and data are SHA256 validated at runtime.  
**Constraint status:** `OFF_CONSTRAINT_DIAGNOSTIC`. There is **no** verified invariant-manifold retraction or invariant certificate in the historical kick protocol.

## The question

How much of the initially written history survives **as the material flow advances** after the ideal *linear part* of the write has been undone? Distinguish the surviving nonlinear state from simple linear transport, observer rotation, and solver error. Importantly this does not assume that the small cubic geometric-descriptor interaction from v18.4 is already material hysteresis.

### Forward-time two-pulse protocol

At physical support `s0=sA(xi=203)` (also `sB(xi=210.8611366031777)`), initialize paired actual nonlinear 11D Q550 histories:

`y_sigma(s0+) = y0(s0) + sigma * amp * v`,  `sigma = +/-1`.

The original rank-1 physical unit vector `v` is frozen. `amp=eps_rel*max(1,||y0(s0)||)`.

Evolve both using the **original nonlinear Q550 support ODE**. The exact unperturbed fundamental matrix `P(s)` is computed with the historical analytic 11D Jacobian along the common baseline, using `P(s0)=I` and the same DOP853 tolerance.

After a delay `L` (expressed as an equivalent xi increment via the original A-section slope), at `s1=s0+L*(sA(433)-sA(390))`, apply a **second, tangent-matched, forward-time anti-kick**:

```
y_echo_sigma(s1+) = y_write_sigma(s1-) - sigma * amp * P(s1) v
r_sigma(s1+) = y_echo_sigma(s1+) - y0(s1)
```

**If the flow were linear, this would cancel the entire perturbation exactly**. Hence nonzero `r_sigma(s1+)` is the finite-amplitude nonlinear material remainder at the reversal time. This is not a physical reversal of external control `Q` or `kappa0`, and it is not backward integration.

Then **advance again** for ages `0,0.25,0.5,1,2` xi-equivalents. At each read time `s`, measure

```
r_sigma(s) = y_echo_sigma(s) - y0(s)
r_even(s) = (r_plus(s) + r_minus(s))/2
r_odd(s)  = (r_plus(s) - r_minus(s))/2
```

The paired even response is expected to start at order amp² for a smooth nonlinear flow; the paired odd component generally starts at order amp³. Both exponents are *tested*, not assumed. The actual full 11D vector is preserved; `dlocal` descriptor vectors are **never** compared componentwise to material state vectors.

### Delay and amplitude sweep

- **Six delays:** `0.25, 0.5, 1, 2, 4, 8` xi-equivalents from initial material kick; unlike observer xi samples, these are *physical support durations* obtained with the historical A map.
- **Three relative amplitudes:** `1e-8`, `3e-8`, `1e-7`, with original v18.4 levels acting as context. A larger third amplitude is included so the quadratic residue can be measured above the numerical floor.
- **Two solver tolerances:** standard `rtol=1e-12, atol=1e-14`; tight `rtol=3e-13, atol=3e-15`.
- **Five forward read ages** per two-pulse experiment: `0,0.25,0.5,1,2` xi-equivalents.
- **72 forward echoed branches** = 6 delays x 3 amplitudes x 2 tolerances x two initial signs, plus 12 shared write branches and common baselines / variational dynamics. Repeated branches and baseline are reused where possible. 36 independently recoverable `(tolerance, amplitude, delay)` checkpoints.

### Memory quantities (not percentages of material constitutive hysteresis)

**Immediate nonlinear state memory relative to the non-erased history:**

```
M_even(L) = ||r_even(s1+)|| / ||(y_write_plus(s1)-y_write_minus(s1))/2||
M_total(L)= mean_sigma ||r_sigma(s1+)|| / ||(y_write_plus(s1)-y_write_minus(s1))/2||
```

**Persistence after attempted cancellation:**

`Pi_even(L,age) = ||r_even(s1+age)|| / ||r_even(s1+)||`.

`Pi>1` may mean ordinary dynamical amplification; it is NOT by itself a memory or hysteresis certificate. We therefore compare the echo to the 11D variational prediction `P(s) P(s1)^(-1) r(s1+)`, conditioning-gated. Also report the control using a naive anti-kick along the **original**, untransported direction `v`, which would generally fail to undo the linear part.

## Run on the original Q550 server

Place the ZIP in `/root/elab1/hysteresis`, then:

```bash
unzip q550_kick203_forward_echo_material_memory_v18_5_bundle.zip
cd q550_kick203_forward_echo_material_memory_v18_5_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_5_server.log
```

Uses `python3.13` by default, overridable with `PYTHON=...`. To rerun *postprocessing only* after the 36 checkpoints exist:

```bash
python3.13 scripts/analyze_echo_v18_5.py
python3.13 scripts/pack_sendback_v18_5.py
```

The wrapper emits `q550_kick203_forward_echo_material_memory_v18_5_sendback.zip`. Full control of original solver locations:

```bash
python3.13 scripts/run_echo_v18_5.py --base-root /root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle \
  --source-bundle-zip /root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip
```

Existing checkpoints are skipped on rerun (unless `--force`). Preserve `checkpoint/` on the server. The expensive original orbit/tail is computed once per invocation, not once per delay, and baseline/tangent/write histories are shared among delay experiments.

## Results emitted after Q550 execution

- `checkpoint/{standard,tight}_eps{1e-8,3e-8,1e-7}_delay{0.25,0.5,1,2,4,8}.npz`: **real material** 11D trajectories, before and after anti-kick; variational predictions; source vectors; physical support times; controls and condition numbers.
- `results/forward_echo_samples.csv`: 180 rows (2 tolerances x 3 eps x 6 delay x 5 age), including even/odd memory, persistence, naive-control residual, solver conditioning.
- `results/echo_amplitude_scaling.csv`: nominal second/third-order fits over the **three epsilon levels** at each delay/read time.
- `results/tolerance_discrepancy.csv`: same physical experiment standard vs tight.
- `results/summary.json`, `results/REPORT.md`, four diagnostic figures.
- `checkpoint/run_provenance.json`, logs, active Python scripts, tests, original v18.4 sendback ZIP with its earlier ancestry, and SHA256 manifests.

## Interpretation gates

1. **Causal gate:** all post-kick states are integrated forward; at zero delay, write and opposite correction restore the identical initial full state within floating arithmetic.
2. **Linear transport gate:** original analytic 11D Jacobian must agree with local RHS finite differences; propagated 11x11 matrices have explicit conditioning gates.
3. **Scaling gate:** `r_even` should approximately scale quadratically and `r_odd` cubically only if data exceed solver tolerance uncertainty; no claim if a component is noise-dominated.
4. **Observer independence:** this v18.5 does *not* reuse geometric B-section frames or dlocal vectors as material coordinates; memory metrics use original material state 11D and Euclidean norm.
5. **History vs hidden memory:** for a deterministic autonomous 11D ODE, identical full states with the same future model produce identical futures. Nonzero echo is a *state-encoded nonlinear residue*, not proof of independent hysteretic memory. Genuine physical hysteresis requires a specified reversible parameter/drive cycle and matched controls/invariants.
6. **Off-constraint limitation:** absent invariant projection/certification, the kick protocol remains a numerical dynamical diagnostic, not a constrained physical experiment.

All documentation and executable comments are in English. The tests include analytical linear/weakly nonlinear echoes; local postprocessing was also dry-run on **external synthetic data**, which are **not distributed as Q550 results**.
