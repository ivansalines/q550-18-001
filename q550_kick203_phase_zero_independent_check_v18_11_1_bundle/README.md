# Q550 v18.11.1 — Independent Phase-Zero Countercheck

The audit is an **addendum** to v18.11, not a replacement or a new material-constraint certificate.

## Scientific gates

1. Read historical standard/tight arrays from the embedded, immutable **original v18.11 user sendback**, descending the nested v18.10 -> v18.9 -> v18.7 -> v18.6 ZIP chain. Verify CRCs and the SHA256 manifests for v18.11 and inner chain. The readers and generator scripts remain present in the original embedded ZIP.
2. Independently integrate the full 11x11 cocycle from **archived J alone**, without reading archived P during propagation: RK4 at one, two and four substeps per xi interval, using cubic interpolation of J. Compare with archived P after propagation. Check expected fourth-order convergence and block structure, including phase-off diagonal coupling.
3. On the server, load the hash-pinned **original nonlinear RHS** and run centered, nonlinear finite-amplitude initial-state tests, perturbing phase and three representative material coordinates. Compare the resulting centered differences with archived P columns. Check that theta evolves at 2pi for off-trajectory initial perturbations. These are off-constraint *numerical diagnostics*, not mechanical perturbations proven physically admissible.
4. Inject artificial errors into the structural-zero gate in unit tests to verify that it detects corruption.

## Run (server)

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_phase_zero_independent_check_v18_11_1_bundle.zip
cd q550_kick203_phase_zero_independent_check_v18_11_1_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_11_1_server.log
```

Default pinned ODE source root: `/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle`.
Override with `Q550_BASE_ROOT=/other/path bash scripts/run_on_server.sh`.
If and only if intentionally skipping original-model checks, use `Q550_SKIP_NONLINEAR=1 bash scripts/run_on_server.sh`; this is **not** nonlinear verified.

Send back: `../q550_kick203_phase_zero_independent_check_v18_11_1_sendback.zip`.

## Limitations

- The last Jacobian row equal to zero and the cocycle block form are exact identities only **conditional on** the stated ODE clock equation `dtheta/ds=2pi`. These tests do not prove an independently known material constraint or metric.
- Finite-difference nonlinear tests are numerical, involve off-constraint probes, and cover chosen initial perturbations, not every phase-space state.
- Classical RK4 verification uses cubic interpolation between measured Jacobian nodes; small residuals verify consistency, not rigorous interval enclosure.
- Coordinate index `10` is zero-based and is the eleventh component.

Python dependencies are listed in requirements.txt; source genealogy and hashes remain inside `inputs/v18_11_original_sendback.zip`.
