# Q550 kick 203 — v18.7 Source-Birth Sign Crossing & Co-Moving Transport Audit

**Protocol:** `2026-10-16-q550-v18.7-material-birth-crossing`  
**Status:** executable Q550 server bundle; numerical Q550 crossing **not precomputed**.  
**Mode:** `OFF_CONSTRAINT_DIAGNOSTIC`. A signed 11D dot product is **not** proof of a material rigid-body torsional flip, physical memory kernel, or constitutive hysteresis.

## Research question

v18.6 independently reconstructed the material quadratic residual of the v18.5 echo from the source
`S2 = 1/2 D²f(y0)[u,u]`, transported through the original 11D variational cocycle `P`.
It established that the source born in `[0,1] xi` contributes a signed **+3.1236%** at `4 xi` and **-13.3297%** at `8 xi` to `K2`, with strong convergence of the v18.6 source reconstruction.

v18.7 locates **every grid-resolved sign crossing** of that previously born contribution at `xi >= 1`, within the complete `xi=1..8` interval, and decomposes its sign-change rate into:

1. deformation / anisotropic strain of the full 11D tangent cocycle;
2. accumulation of newly generated material source.

It additionally distinguishes the old vector's laboratory direction, its direction after orthogonal polar co-rotation, and its initial-source-frame direction (the exact inverse-cocycle gauge). These frames are *diagnostic choices*, not assumed rigid-body coordinates.

## Mathematical construction

Let `s=s0+xi*slope`, `P' = J(s) P`, `P(s0)=I` and `g(s)=P(s)^(-1) S2(s)` be the original v18.6 source integrand. Define

```text
A_early = integral_s0^(s0+1xi) g(s) ds
A_total(xi) = integral_s0^s g(s') ds'
E(xi) = P(s) A_early
K2(xi) = P(s) A_total(xi)
L(xi) = K2(xi) - E(xi)
N(xi) = E(xi)·K2(xi) = A_early^T G(xi) A_total(xi), G=P^T P.
pct_early(xi) = 100 N(xi) / ||K2(xi)||².
```

The crossing is **exactly** `N(xi*)=0` in this chosen 11D Euclidean metric. For normalized angle it is an old-vs-total 90° crossing, not necessarily a 180° old-vs-new flip.

In terms of the original xi coordinate (not the ODE support time),

```text
N'(xi) = A_early^T G'(xi) A_total(xi) + A_early^T G(xi) A_total'(xi)
G'(xi) = slope * P^T [J^T+J] P
A_total'(xi) = slope * g(xi)
```

**First term:** full cocycle strain/deformation, including metric rotation effects from non-isometric stretching. **Second term:** new material source birth. They sum algebraically to the complete scalar sign-change rate. The split is well-defined for the specified original 11D coordinates and inner product; it is **not invariant under arbitrary non-orthogonal changes of state units**.

The orthogonal polar factor `R` in `P=R U` supplies a reproducible co-rotating diagnostic: `R^T E` compared to `E` in lab coordinates. Applying `R^T` to *both* E and K2 does not change their mutual angle or crossing. `P^-1 E=A_early` is exactly constant by construction and is a non-isometric source-coordinate transformation.

## What is reused, and what the server computes

- `inputs/v18_6_original_sendback.zip` is frozen by SHA256 and retains **all original v18.6 code, Hessian/source measurements, 641-node 11D arrays, v18.5 checkpoint history, v16/v13 original references and manifests**. The sendback's own 32 SHA256 entries are verified.
- `reference/q550_state_impulse_response_v13.py` and `reference/history_bridge_v16.py` are byte-identical archival copies from v18.6; they validate original Q550 source dependencies and provide the same 400-period, 70-period-tail baseline construction and analytic Jacobian.
- The only Q550 server computation is **replaying the unperturbed baseline and the 11×11 variational cocycle** separately for `standard` and `tight`; no new +/- kicks, Hessian calculations, or echo simulations.
- The 641 archived source points are spaced by `0.0125 xi`. The fine integral uses Simpson cells `0.025 xi`, with consistent half-cell values and cubic Hermite interpolation between nodes. Roots are bracketed on the sampled grid and refined by Brent's method using dense `P(s)`.
- Controls: source Jacobian step halved (archived from v18.6), coarse Simpson cells of `0.05 xi`, standard/tight agreement, full-vector regressions at all six v18.6 delays, scalar rate finite differences around each root, complete 11D closure, and inverse-cocycle constancy.

**Important:** interpolation can refine a root *within* a sign-changing interval but cannot exclude an unresolved pair of crossings between sampled nodes. Source Hessian estimates remain finite-difference diagnostics, not interval enclosures. A separate constraint-preserving perturbation protocol is still needed to claim a physical memory law.

## Server command

Copy the ZIP to `/root/elab1/hysteresis` and run:

```bash
unzip q550_kick203_source_birth_sign_crossing_v18_7_bundle.zip
cd q550_kick203_source_birth_sign_crossing_v18_7_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_7_server.log
```

The runner defaults to `python3.13`; set `PYTHON=python3` if necessary. External original Q550 dependencies must already exist in their historical server locations (exact hash checks fail if missing/modified).

Explicit original-path overrides can be forwarded:

```bash
bash scripts/run_on_server.sh \
  --base-root /root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle \
  --source-bundle-zip /root/elab1/hysteresis/q550_x55_multiloop_transport_runner_v1_1_bundle/q550_x55_multiloop_transport_x50000_rebased_checked_v1/q550_x55_multiloop_transport_law_v2_bundle.zip
```

The final sendback is:

`q550_kick203_source_birth_sign_crossing_v18_7_sendback.zip`

To resume a partial run, rerun the same command. The `standard` and `tight` checkpoints are reused; use `--force` only to regenerate them. **Existing checkpoints are trusted after version checks:** if original source dependencies, source hashes, or solver settings change, explicitly use `--force`.

After both tolerances are done, regenerate only analyses/figures, without another Q550 replay:

```bash
python3.13 scripts/analyze_sign_crossing_v18_7.py
python3.13 scripts/pack_sendback_v18_7.py
```

## Outputs

- `checkpoint/{standard,tight}_sign_crossing.npz`: complete sampled 11x11 P and J, material 11D early/late/total vectors, polar co-rotated early 11D vectors, source-frame cumulative integral, roots, v18.6 regression errors.
- `results/{standard,tight}_dense_projection.csv`: 561 post-1xi observations each, signed percentages, old/total angles, lab-vs-polar angles, three derivative columns and conditioning.
- `results/{standard,tight}_crossing_roots.csv`: locations, brackets' signs, source/transport contributions to derivatives, and finite-difference check.
- `results/{standard,tight}_crossing_summary.json`: all roots, source-Hessian and quadrature shifts, six historical 11D anchors.
- `results/summary.json`, `results/REPORT.md`, and three PNG plots: dense signed crossing, source vs transport rate, polar vs lab co-moving angles.
- Full archival ZIP with script genealogy, nested v18.6/v18.5 history, SHA256 and raw checkpoints.

## Verification status before server

Ten local tests passed, including full synthetic 11D source-to-crossing regression. Offline v18.6 history SHA256 and all 32 embedded member hashes verified. **The numerical v18.7 Q550 zero and its decomposition remain to be measured on the server.**

Comments, technical report labels, and source-file documentation are in English.
