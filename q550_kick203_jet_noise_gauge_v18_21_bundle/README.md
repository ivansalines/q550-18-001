# Q550 v18.21 — Jet Sensitivity, Noise Floor & Observer Gauge Audit

## Main question

Does the hidden tenth material direction detected by the third jet remain *operationally* observable with finite time windows, hypothetical measurement noise, and coordinate/observer gauge changes? This continues the **authentic original-server v18.20 sendback**, not just the prior frozen bundle.

## Exact server commands

Copy `q550_kick203_jet_noise_gauge_v18_21_bundle.zip` to `/root/elab1/hysteresis`:

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_jet_noise_gauge_v18_21_bundle.zip
cd q550_kick203_jet_noise_gauge_v18_21_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_21_server.log
```

Expected output: `/root/elab1/hysteresis/q550_kick203_jet_noise_gauge_v18_21_sendback.zip`.

For original model modules in a different folder:

```bash
bash scripts/run_on_server.sh --base-root /path/to/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle \
  2>&1 | tee v18_21_server.log
```

Frozen-only (no nonlinear original ODE): `bash scripts/run_preflight.sh`. It does **not** produce a server-certified sendback.

## Inputs and mathematical distinction

The original v18.20 server sendback is embedded in `inputs/`. It includes real original-ODE nonlinear trajectories and the genealogy back to the kick203 source. `scripts/checkpoint_io_v18_21.py` checks CRC, exact inner root, immutable SHA256 manifest, nonlinear-completion status and the embedded predecessor; no orbit regeneration.

At the reference point `xi=4`, each fixed diagnostic projection `H` has jet matrices:

```
M0=I
M1=A
M2=A' + A^2
M3=A'' + 2 A' A + A A' + A^3
O_k = stack_{j=0..k}(tau^j H Mj)
```

Here `A=ds/dxi * J[0:10,0:10]` and derivatives of `A` are total along the original 11D base orbit. The phase-fixed 10D *variation* is invariant; its underlying trajectory is still phase-driven, not a 10D autonomous system. Jets are numerical and belong to selected observer gauges, not to a proven material constraint manifold.

The audit contains **three different, noninterchangeable experiments**:

1. **Mathematical jets with synthetic noise:** independent Gaussian noise injected directly into each *scaled jet* component (`sigma=1e-6,1e-4,1e-3`). This is an optimistic thought experiment, NOT experimentally measured sensor noise. Noise can be much worse after numerically differentiating position data.
2. **Nine directly sampled observables** across symmetric time windows (halfwidths `.05,.1,.2,.5,1` xi), reconstructed from the already-certified full10D tangent propagator. This avoids taking numerical derivatives of noisy data. An isotropic fixed per-sample noise level defines a comparison signal-to-noise ratio; changing window duration genuinely changes information.
3. **Gauge and covariance audit:** invertible 10D coordinate and output transformations. Raw SVD condition numbers are not invariants. After transforming observables and noise covariance coherently, physical predictions remain unchanged. Likewise, multiplying jet blocks by `tau^k` while also transforming their uncertainties leaves whitened Fisher information unchanged. Treating all newly scaled jet blocks as identically noisy introduces a *new noise assumption*, and cannot itself manufacture measured information.

`run_on_server.sh` adds a **fresh original-ODE spot check**, six nonlinear trajectory pairs built from three v18.20 hidden seeds at two amplitudes. It authenticates original model Python source hashes, checks their nonlinear response against the archived tangent cocycle at `xi=4` to `4.5`, and verifies phase matching. This is a real ODE re-evaluation, not a claim of measured instrument noise or physical hereditary memory. It does not integrate 400 orbit periods.

## Files and controls

- `scripts/preflight_v18_21.py`: actual archived Jacobian/cocycle, jet SVD, 256 Monte Carlo noise draws per case, direct time window, gauge/covariance controls.
- `scripts/nonlinear_server_v18_21.py`: independent original Q550 model spot check and model SHA256 pinning.
- `scripts/noise_gauge_math.py`: testable gauge/SVD/noise formulas.
- `results/`: CSV, JSON reports and plotted results. Any `original_ode_check_*` files are **server-only**.
- `scripts/build_manifest.py` and `scripts/pack_sendback.py`: hash immutable files only; logs in progress excluded; reject any sendback without original-ODE check.
- `tests/test_noise_gauge.py`: numerical and synthetic invariance tests.

## Explicit limitations

This audit cannot certify global nonlinear injectivity, physical constraints, a physically preferred metric, real sensor-noise statistics, or hereditary constitutive memory. Numerical observability full rank should not be confused with stable reconstruction or with a physical proof of absence/presence of hysteresis. Coordinate rescalings are passive gauges only when all metrics/noise models are transformed consistently.

The ZIP has **exactly one top-level directory**. Paths, README and shell command are checked after clean extraction.
