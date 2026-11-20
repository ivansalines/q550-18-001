# v18.3 — Material-to-Observer full11D Bridge

## Source alignment
- 47 identical xi samples; physical B crossing xi=210.861136603177698.
- B support kick s0=693.612977136753784; hquarter stencil half-width=0.00347131494 xi.
- Clean left/right: 17/19; centered straddling samples: 11.

## Exact finite-history identity

`Delta_dlocal_cent = S_A + S_B + S_M + Interaction` (conditional midpoint factor secants).
`Delta_dlocal_order2 = Delta_dlocal_cent + ObserverStencilCorrection`.

These are exact vector identities, not causal identification of an intrinsic torsion.

## Signed residual away from centered-stencil contamination

- standard, eps=1e-08: pre positive=True; post negative=True; post median=-1.20413405e-04.
- standard, eps=3e-08: pre positive=True; post negative=True; post median=-3.55287408e-04.
- tight, eps=1e-08: pre positive=True; post negative=True; post median=-1.20411769e-04.
- tight, eps=3e-08: pre positive=True; post negative=True; post median=-3.55038424e-04.

## Controls

- source_frame_descriptor: 0.0000e+00
- archived_shapley: 5.7259e-17
- factor_split: 0.0000e+00
- observer_shapley: 5.6207e-17
- observer_identity: 1.1301e-17
- material_causality_before_B: 0.0000e+00

## Cautions

The directional `+/-` material impulse is an on-flow perturbation; its invariant manifold is not independently certified.
The section-A and section-B states are sampled at different physical support coordinates.
The descriptor and state sensitivities are not interchanged, subtracted, or implicitly dot-producted.
The factor midpoint construction does not identify the complete Jacobian from material states to dlocal.
Centered stencils near B crossing mix both sides of the restart; the one-sided sample exactly at the crossing is censored.
An observable residual sign flip is not automatically a material torsional flip.
