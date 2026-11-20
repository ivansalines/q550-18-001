# Q550 — v18.23 Model-Blind Hidden-State Reconstruction & Predictive Closure Audit

## Central falsifiable question

Do past measurements of a reduced observation predict an **unseen future** better
than its current value, when the predictor does **not** receive the original ODE,
Jacobian, cocycle, or true hidden material state?

The estimator is deliberately **sensor-only**, not an equation-informed filter.
`blind_math.py` contains the actual estimator and cannot import model operators.
The authentic original Q550 ODE is used **exclusively to generate** independent
training/validation/test trajectories on the server and to evaluate their future.
The locally computed preflight is a separate frozen-tangent experiment.

## Exactly reproducible server commands

Place the archive in `/root/elab1/hysteresis`, then run:

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_model_blind_predictive_closure_v18_23_bundle.zip
cd q550_kick203_model_blind_predictive_closure_v18_23_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_23_server.log
```

Result: `/root/elab1/hysteresis/q550_kick203_model_blind_predictive_closure_v18_23_sendback.zip`.

If the original pinned model was relocated, use:

```bash
bash scripts/run_on_server.sh --base-root /path/to/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle \
  2>&1 | tee v18_23_server.log
```

Frozen-only local validation: `bash scripts/run_preflight.sh` (no original
nonlinear Q550 server certification).

The ZIP contains exactly ONE folder root. The live `tee` log and transient
`.runtime` directory are excluded from immutable SHA256 manifests.

## Experimental protocol

- Source: **authenticated original v18.22 server sendback** with SHA256 and
  historical chain. Original Q550 kick 203 state and **pinned original ODE** on
  server; `q550_state_impulse_response_v13.py` checks the original dependencies.
- Meeting time `xi=4.0`, fixed original phase. Observation histories at 17
  distinct times in `[3.5,4.0]`, and future targets at `xi=4.2,4.5,5.0`.
- Reduced sensors: `material3`, `raw4`, `raw6`, all declared *diagnostic* gauges.
  The choice of `material3` was ODE-aware historically; sensor-only status
  refers to **the fit and deployment**, not the design of that sensor.
- 96 training directions, 32 independently generated validation directions,
  36 independent random test directions plus one specially selected hidden
  direction **never seen in training** per sensor.
- Nonlinear training and validation amplitude `epsilon=1e-5`, independent
  held-out amplitude `epsilon=3e-6`, trajectories compared to a reference orbit
  and normalized by epsilon. Perturbation amplitudes are known to the data
  preparation process, not the predictor.
- All predictors use only noisy reduced measurements in one of four forms:
  (i) instantaneous endpoint; (ii) recent 0.10 xi; (iii) complete 0.50 xi
  history; (iv) separately shuffled time order in each trajectory (destructive
  timing control). Same sample-noise specification and clean future targets.
- Ridge penalty is selected by validation ONLY. Test is held out throughout.
  `sigma=[0,1e-4,1e-3]` refers to **normalized perturbation observations**,
  not known physical sensor noise.
- Future prediction is the primary criterion. SVD of observed history is a
  latent-dimensionality diagnostic, not a claim of unique reconstruction of
  physical hidden coordinates.
- No true hidden state labels are supplied to the predictor. A hypothetical
  state-labeled decoder would be supervised, hence **not** genuinely blind;
  this audit intentionally does not train one.

## Audit and scope distinctions

`frozen_*`: the authentic previous Q550 cocycle generates synthetic tangent
observations; **NOT** a fresh nonlinear validation. It tests split discipline,
sensor-only fitting, regularization and history ablation.

`original_*`: only the user's server can generate independent original-ODE
training and test trajectories and issue `original_ode_blind_summary.json`.
Packaging refuses a server-certified sendback without that result.

If the full-history model predicts an unseen future substantially better than
instantaneous sensors, it establishes **observable predictive value of history
for that reduction**, not an independent constitutive memory or unique hidden
state. Since the full 11D ODE is deterministic, distinct hidden states can
fully explain different futures without any additional hereditary force.

## Genealogy and complete executable provenance

- `inputs/v18_22_original_server_sendback.zip`: verified original server
  lineage with Python scripts, reports and nested historical bundles.
- `scripts/checkpoint_io_v18_23.py`: verifies prior SHA and original ODE status.
- `scripts/frozen_preflight_v18_23.py`: frozen source-propagator data generation.
- `scripts/nonlinear_server_v18_23.py`: independent original ODE integration,
  strict source hash checks, disjoint train/validation/test trajectories.
- `scripts/blind_math.py`: model-blind estimator; never touches ODE,
  Jacobian, cocycle, hidden states, simulator internals.
- `scripts/blind_experiment.py`: controlled splits, matched noise,
  out-of-sample performance comparisons, diagnostics.
- `tests/test_blind_v18_23.py`: leakage checks and independent synthetic tests.
- `scripts/plot_v18_23.py`, `results/*`: plots and data, marked frozen/original.
- `scripts/build_manifest.py`, `scripts/pack_sendback.py`: immutable file hashes,
  fail-closed server packaging, live-log exclusion.

## Limits

This does not certify a unique physical metric, material twist or hysteresis.
Two different hidden 10D states with the same reduced observation and phase
may exhibit different predictions; this may reflect ordinary incomplete-state
observation. Predictions are made relative to a known reference trajectory,
not from uncalibrated absolute sensor readings. Local coordinate and sensor
choices matter, as does synthetic noise, and results must be validated on
independent trajectories outside the fitting interval.

## Sensor-learned latent delay states (additional diagnostic)

`blind_math.delay_latent_rank_scan()` learns principal delay-coordinate factors
from **training observations only**, then validates a linear readout to future
sensor values. It compares latent ranks 3, 4, 6, 8, 10, 12, 15, 20;
**rank is selected by validation data, not test data**. No physical state labels
are supplied. The fitted latent coordinates are not identified with material
state variables and a noise-influenced best numerical rank can exceed 10.
