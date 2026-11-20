# Q550 v18.10 — Material 2-Plane, Source Escape and Original-ODE Inventory

**Protocol:** 2026-10-18-Q550-kick203-v18.10-plane-constraint-gate

This bundle is an executable continuation of the actual v18.9 sendback. It embeds the entire v18.9 sendback, including its original v18.7/v18.6/v18.5 checkpoint genealogy and solver scripts, together with all newly executed Python scripts, tests, original-v18.9 historical readers, manifest and report. The original 11D state-space source remains **OFF_CONSTRAINT_DIAGNOSTIC** until independent constraint/metric provenance is provided.

## Mathematical question

The Euclidean v18.7 sign-crossings at xi ~4.1352, 4.9973, 7.1744 are metric dependent. Instead study two transported 11D vectors E=P a (the early source born by xi=1) and K=P A (the total quadratic source); a=A(1) and A'=slope*g.

Define B=E wedge K. Under fixed invertible changes of 11D coordinates the bivector transforms by the induced exterior-square map; B=0 (exact collinearity) is coordinate invariant. Its xi-derivative separates algebraically into the cocycle-transport term and the source-birth term:

    dB/dxi = slope [(J E) wedge K + E wedge (J K)] + slope [E wedge (P g)]

A second rank test follows a wedge A wedge g: if nonzero the fresh source enters a third linearly independent state-space direction. All **numerical** ranks use an explicitly declared 1e-8 SVD gauge; their thresholded classification is not strictly GL-invariant. Norms, angles, sine values and relative magnitudes are also NOT intrinsic without a metric.

The two-plane B=0 at xi=1 by construction. Never count this initial point as a genuine dynamical degeneracy. A source-frame plane or an off-constraint 11D phase space is not automatically a **physical body plane**.

## Server

From /root/elab1/hysteresis after copying this ZIP:

    unzip q550_kick203_material_plane_constraint_v18_10_bundle.zip
    cd q550_kick203_material_plane_constraint_v18_10_bundle
    bash scripts/run_on_server.sh 2>&1 | tee v18_10_server.log

This runs synthetic unit tests, reconstructs real 641-node standard+tight exterior measurements and, by default, verifies the original hash-pinned Q550 ODE using the v18.6 archived kick state. No 400-orbit baseline solve. It writes one directory higher:

    q550_kick203_material_plane_constraint_v18_10_sendback.zip

To explicitly omit original ODE verification when pinned modules are unavailable:

    Q550_SKIP_REPLAY=1 bash scripts/run_on_server.sh

To set a different original model location:

    Q550_BASE_ROOT=/path/to/original/q550/model bash scripts/run_on_server.sh

## Artifacts

- `results/summary.json`: per tolerance exterior-gauge diagnostics and coordinate-change tests.
- `results/material_plane_profile.csv`: measurable xi series and rank estimates.
- `results/REPORT.md`: numerical interpretation and limitations.
- `results/ode_structure_summary.json`: reloaded *original ODE* module SHA256 and Jacobian regression.
- `results/ode_coordinate_inventory.csv`: Jacobian-row and RHS coordinate screens, not mechanical identifications.
- `results/ode_jacobian_replay.csv`: original-v18.7 replay regression.
- `SHA256SUMS_OUTPUT.json` and nested immutable input history.

## Scientific limits

This is not a rigorous interval enclosure, physical constitutive hysteresis, a certified ODE energy law, or a material-frame torsion proof. The ODE replay verifies model identity and its local Jacobian but does not produce a material G or nonlinear C on its own. The fact of exact exterior rank is coordinate-invariant; measurement thresholds are not. Report all null results and ill-conditioned segments.

## Visual profiles

Three regenerated figures (source-frame normalized plane separation, source-frame three-direction volume, and Euclidean laboratory plane separation) are produced by the included `scripts/plot_plane_v18_10.py`. All magnitudes are explicitly gauge-dependent, even though exact wedge vanishing is covariant.
