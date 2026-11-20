# Q550 v18.20 — Reduced Jet Hierarchy & Hidden-State Observability Audit

**Goal.** Determine which hidden tangent-state directions become locally distinguishable after observing the material state projection and its first, second, and third derivatives. Continue v18.19 using the **original server sendback** and pinned 11D Q550 ODE. No 400-orbit recomputation.

## Run on the server

Copy `q550_kick203_reduced_jet_hidden_observability_v18_20_bundle.zip` into `/root/elab1/hysteresis`, then:

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_reduced_jet_hidden_observability_v18_20_bundle.zip
cd q550_kick203_reduced_jet_hidden_observability_v18_20_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_20_server.log
```

The result archive appears one directory above:
`/root/elab1/hysteresis/q550_kick203_reduced_jet_hidden_observability_v18_20_sendback.zip`.

If the original ODE modules are installed elsewhere, pass `--base-root /path/to/original/Q550/bundle` after `scripts/run_on_server.sh`. The module files are independently checked against pinned SHA256 from the original v13 script. **Do not relabel the frozen preflight as nonlinear server confirmation.** If the original modules are unavailable, `bash scripts/run_preflight.sh` performs only the frozen audit, and does NOT package a scientific server sendback.

## Mathematical setup

Write the phase-driven 11D ODE as material `z' = f_z(z,theta)` and phase `theta'=2pi ds/dxi`, with phase-fixed perturbations obeying `delta z'=A(xi) delta z` where `A=ds/dxi*J_{z,z}`. Here the actual phase convention is inherited from the original Q550 code. The base trajectory still evolves in phase: **phase-fixed variations are not a 10D autonomous system.**

For any fixed diagnostic projection `H`, set `u(xi)=H z(xi)`. The tangent sensitivities of its first four time jets at `xi_m=4` are:

```
M0=I
M1=A
M2=A' + A^2
M3=A'' + 2 A' A + A A' + A^3
O_k = vertical_stack(H M0, tau H M1, ..., tau^k H Mk)
```

The derivatives `A'` and `A''` are *total along the original base trajectory*, computed from 5-point centered stencils on the archived Jacobian samples. Scans at xi=2,3,4,5,6,7 with finite-difference stride=1,2,4 and numerical rank thresholds `1e-7,1e-8,1e-9`; jet-time scales `tau=.05,.1,.2`. SVD-based rank is a **local linear numerical result** and cannot by itself certify global nonlinear state recovery, robust noise-conditioned inversion, a privileged physical metric, or constitutive hysteresis. Full rank does not imply good conditioning. Material3 is a constructed H from three archived source directions at xi=4; raw4/raw6 are coordinate projections.

## Experimental design

* Frozen preflight authenticates each nested generation back through v18.18 and v18.7, inspects exact rank and singular values of the phase-fixed material Jacobian jets, computes 7 hidden directions with jets 0, 1 or 2 matched, and compares third-order Taylor predictions with the archived cocycle near xi=4. A *local frozen result* is included in this bundle and labeled distinctly.
* **Server-only validation:** integrate the original nonlinear 11D ODE from the v18.6 kick203 state to xi=4, reproduce the original Jacobian independently from the RHS; then integrate y+(4)=y(4)+eps(v,0), y-(4)=y(4)-eps(v,0) *independently* to xi=3.9 and 4.5, eps=`3e-6,1e-5`. Verify first departure orders of the reduced observations and tangent comparisons. This tests 7 seed directions times 2 amplitudes = 14 nonlinear trajectory pairs.
* The matching of jets is **linearized**, not exact at finite eps. The centered nonlinear difference introduces higher-order terms in eps. Equality of full 11D state determines identical deterministic future; this is a test of the **chosen reduction**, not an independent causal memory variable.
* Original-model physical constraints have not been demonstrated. The experiment remains **OFF-CONSTRAINT** until a physically established admissible manifold is supplied.

## Results / provenance

`results/frozen_jet_summary.json` and related CSV/PNG files are frozen source-only tests. `results/nonlinear_jet_summary.json`, `results/nonlinear_jet_validation.csv`, and optional `results/nonlinear_jet_growth.png` are generated only on successful original-ODE server run. `inputs/v18_19_original_server_sendback.zip` holds the actual v18.19 nonlinear results and embedded historical genealogy, preserving functional executability. `SHA256SUMS_OUTPUT.json` hashes immutable files; active `*.log` files and `.runtime` are excluded to avoid the previous `tee` synchronization error. The ZIP must have exactly one root directory and follow the commands above.
