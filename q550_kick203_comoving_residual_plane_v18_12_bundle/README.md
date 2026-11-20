# Q550 kick203 — v18.12 Co-moving Residual Material-Plane Transport Audit

Goal: quantify the motion of the material memory plane after subtracting the shared homogeneous transport in the phase-fixed 10D variational sector.

## Source and reproducibility

- `inputs/v18_11_original_sendback.zip`: original FULL v18.11 sendback, including ODE replay, with an unbroken nested chain back to v18.6. All nested manifests verified by checkpoint loader.
- `reference/`: original v18.9 checkpoint reader and math references, carried forward to make v18.12 standalone.
- `scripts/`: actual computation, chart rendering, provenance verification, sendback creation.
- `tests/`: independent mathematical sanity checks, including negative controls with no residual plane turning.
- `results/`: committed frozen-data report and tables from current local run; server should reproduce them.

## Run

    cd /root/elab1/hysteresis
    unzip q550_kick203_comoving_residual_plane_v18_12_bundle.zip
    cd q550_kick203_comoving_residual_plane_v18_12_bundle
    bash scripts/run_on_server.sh 2>&1 | tee v18_12_server.log

Send back `../q550_kick203_comoving_residual_plane_v18_12_sendback.zip`.

## Crucial interpretation restrictions

`Phi^{-1}` removes the homogeneous **linearized transport**, not every geometrical frame effect. The Euclidean source-frame angles/path lengths are gauge-dependent. Algebraic non-collinearity and wedge covariance are invariant under invertible coordinate changes. A nonzero turn indicates the accumulated source changes the co-moving plane; it does **not** certify a physical torsion mechanism, non-holonomic constraint, or hysteresis from path-dependent repeated states. `theta` is a prescribed evolving phase, not a constant physical coordinate.

No 400-orbit replay. This version uses the on-server replay already certified in v18.11 and frozen full-source data, so it can be repeated quickly.
