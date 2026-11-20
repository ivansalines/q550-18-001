# Q550 kick203 — v18.3 Material-to-Observer Full11D Bridge

## Objective

Separate the causal **11D material-history sensitivity** at physical support samples
from the **11D geometric two-section descriptor** built by A and B tangent/curvature
frames and their intervening full11D cocycle `M`. Test whether the negative
signed geometric residual after B meets the kick is a robust descriptor output,
or an artifact of the centered B derivative stencil.

## Numerical sources and provenance

This complete standalone bundle embeds `inputs/q550_kick203_on_shell_full11d_material_transport_v18_2_sendback.zip`.
The latter archive embeds the original complete v18.1 sendback, which in turn
embeds the original complete v17.9 sendback. The v18.3 runner extracts the
archived 11D material, observer, and historical factor arrays from these exact
nested inputs with NO new ODE solver run. The pinned outer input SHA256 is
specified in `SOURCE_MANIFEST.json`; archived historical source scripts
are included in `reference/` for auditability and further replay.

The **material coordinate 11D** and **geometric descriptor 11D** remain distinct
quantities even if both are represented by 11-component vectors. No direct
subtraction or dot product of these different observables is permitted.
The A and B section states are measured at different physical support values
`s_A(xi)` and `s_B(xi)`; they are not simultaneous state samples.

## Exact bridge construction

For each `xi`, tolerance `standard/tight` and epsilon `1e-8/3e-8`, recover from
v17.9 the actual opposite history factors
`(BA-,BB-,M-)` and `(BA+,BB+,M+)`. Let `Fbar=(F++F-)/2` be their componentwise
midpoint. Reevaluate the historical nonlinear descriptor `d(BA,BB,M)` and
form the exact **directional factor secants**:

```
S_A = d(BA+,BBbar,Mbar)-d(BA-,BBbar,Mbar)
S_B = d(BAbar,BB+,Mbar)-d(BAbar,BB-,Mbar)
S_M = d(BAbar,BBbar,M+)-d(BAbar,BBbar,M-)
I_ABC = [d(F+)-d(F-)] - (S_A+S_B+S_M)
```

Thus the historical centered descriptor's observed `+/-` response has the exact
identity `Delta_centered = S_A+S_B+S_M+I_ABC`. Here the midpoint need not be a
physical frame or trajectory; it serves ONLY as a coordinate for conditional
secants. We also independently reconstruct the exact symmetric 3-factor
Shapley allocation, validating it against archived v17.9 factors.

The separate **observer correction** is defined directly in the same geometric
11D descriptor space:

`Delta_observer_correction = Delta_v18.1_one_sided_second_order - Delta_v17.9_centered`

`Delta_v18.1_one_sided_second_order = S_A + S_B + S_M + I_ABC + Delta_observer_correction`

Unlike a complete Fréchet material-to-descriptor Jacobian, these are
**one-direction historical secants** from the existing plus/minus experiment.
They do not identify independent responses to arbitrary perturbations of all
11 material coordinates. The v18.1 A frame and M cocycle remain frozen when
only the B derivative stencil changes, as in the archived experiment.

The material histories `u_A=(y_A+-y_A-)/(2 amp)` and similarly `u_B` are carried
along as causal controls, with their original physical-support coordinates.
We require `u_B=0` before B reaches the restart, and check the v18.2 linear
transport predictor separately. No artificial conversion of material to
geometric descriptors is assumed.

## Causal masks and validity

The centered v17.9 B observer hquarter spans `s_B+-hquarter`, where
`hquarter=0.0004448784722222222` support units. Its contaminated strip is
`|s_B-s0| <= hquarter` (about `0.00347131494 xi` on each side).
The exact event sample is censored for the v18.1 non-straddling observer.
Signed comparisons on clean left/right samples are audited separately.

This is an **on-flow** history audit and not a physical invariant certificate,
rigorous interval enclosure, or proof of an actual torsional flip.

## Already-computed results and repeatability

The `results/kick203_v18_3` folder ships with:

- `summary.json`, `REPORT.md`, `console.log`
- `bridge_samples.csv` — 188 rows, 47 xi samples times four experiments
- `robust_sign_audit.csv` — sign/gain/interaction statistics on clean samples
- `bridge_full11D_arrays.npz` — actual aligned material, geometric descriptor,
  Shapley, secant, interaction, and observer-correction arrays
- Four plots showing signed correction, A/B/M compensation, causal state
  response and conditional secant non-additivity

Reproduce everything **without the 400-orbit solve**:

```bash
cd q550_kick203_material_to_observer_bridge_v18_3_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_3_server.log
```

Override the Python interpreter with `PYTHON=python3 bash scripts/run_on_server.sh`.
The same wrapper regenerates and packages
`q550_kick203_material_to_observer_bridge_v18_3_sendback.zip`.

If packaging fails after the audit succeeds, run
`python3.13 scripts/pack_sendback_v18_3.py` alone.

## Interpretation boundary

A robust *signed geometric residual* after B crosses the source is observed,
with nearly cancelling A and B Shapley contributions and a small negative M
term. The A and B contributions are conditional Shapley allocations, not
individual experimentally separable physical forces. The nonlinearity of the
three-factor map must be retained explicitly. Material state continuity of
A across the B event and causality at B disfavor interpreting the abrupt
observer switch as a material state flip at B. A full material->descriptor
linear operator, independently perturbing the 11D state basis directions,
would be needed for stronger causal differential interpretation.
