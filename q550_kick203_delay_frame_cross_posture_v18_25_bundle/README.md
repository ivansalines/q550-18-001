# Q550 v18.25 — Delay-Frame Alignment & Cross-Posture Latent Transport Audit

## Question and scope

Does the predictive advantage of a 20-component time-delay latent representation survive **sensor-history-only** alignment of different orbital postures? Can a rank-10 predictor transport adequately once a reduced frame calibration is learned? Are observed gains robust to synthetic observation noise, and do they disappear if paired calibration identities are deliberately wrong?

**No constitutive memory, physical extra 10D state, or invariant torsional flip is assumed.** This is a local off-constraint diagnostic for the material 10D sector of an 11D Q550 system (11th coordinate: phase). Section changes from xi=4 to 3.5/4.5 change *both* phase and material posture: not pure phase experiments.

## Exact server instructions

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_delay_frame_cross_posture_v18_25_bundle.zip
cd q550_kick203_delay_frame_cross_posture_v18_25_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_25_server.log
```

The expected completed sendback is:

`/root/elab1/hysteresis/q550_kick203_delay_frame_cross_posture_v18_25_sendback.zip`

If the SHA256-pinned original ODE tree was relocated, use

```bash
bash scripts/run_on_server.sh --base-root /absolute/path/to/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle \
    2>&1 | tee v18_25_server.log
```

For the checkpoint-only preflight (which **cannot** certify an original nonlinear result): `bash scripts/run_preflight.sh`. The packager refuses to package without `results/original_ode_transport_proof.json` with the correct nonlinear completion status.

## Data separation and observer's blindness

- **Reference predictor:** trained on 96 nonlinear perturbations at xi=4; selects ridge/latent rank only on 32 reference validation trajectories; evaluated on 36 independent directions plus one held-out matched hidden-direction probe. Train/validation epsilon = 1e-5, test epsilon = 3e-6.
- **Calibrated section transport:** uses 32 *paired* training-direction histories at xi=3.5 and 4.5 with the matching histories from xi=4. Pairing identity is known from calibration design. No future outputs of these shifted calibration trajectories are passed to the learning algorithm, no test input is used to select the alignment, and no future of shifted test trajectories is used in any fit.
- **Calibration train/holdout:** 24 of these 32 direction pairs fit the pure lag and affine regularization selection; 8 choose the lag or affine penalty using **past/current histories only**. The full 32 may then fit the final mapping.
- **Held-out shifted test:** the same 36 unseen random directions at all sections and one previously isolated hidden-direction probe transported through the archived cocycle solely by the separate data-generating oracle. These are not treated as locally jet-matched at shifted sections.
- **Inputs given to blind predictor/aligner:** three-, four-, or six-component reduced histories of 17 past-to-present samples at offsets [-0.5,0] xi; calibration correspondence; and observed future reduced targets at offsets [0.2,0.5,1.0] xi **only at xi=4 for training and validation**. No model state, governing RHS, tangent propagator, Jacobian, or true material coordinates are supplied.
- **Noise:** independent Gaussian noise on each amplitude-normalized observed history, sigma in {0,1e-4,1e-3}, no noise on the target future. These are not physical experimental uncertainties.

## Comparisons and invariant/null controls

For each sensor `material3`, `raw4`, `raw6`, for both ranks **10 and 20**, three noise levels, and three sections:

1. `none`: reference predictor without a transported frame (zero-shot).
2. `past_lag`: single time shift of past sampling positions, selected solely by paired calibration histories. The terminal sample is **bitwise preserved** in this intervention (the shift acts only on earlier samples).
3. `orthogonal`: orthogonal Procrustes map between calibrated reduced *latent history frames* (not a physical rotation of the Q550 body).
4. `lag_orthogonal`: previous two mechanisms together.
5. `affine_ridge`: a regularized learned affine map between the two time-delay frames; regularization selected on held-out **calibration histories**, not held-out prediction results.
6. `broken_pair_orthogonal`: negative control with independently permuted calibration correspondence. A false positive must be detectable.

The learned latent frame transformations, unlike the pure past-lag operator, may transform the representation of the present along with the past. They are **not** a claim of preserving the raw endpoint components. More importantly, they should not be interpreted as physical rigid rotations.

A common global orthogonal change of latent gauge leaves the correctly transformed predictions unchanged; this is covered by an explicit unit test. The results record every candidate, including alignment methods that **worsen** the held-out future prediction.

> Identifiability caveat: calibrated transport uses a richer experimental protocol than zero-shot transfer (identities of paired training perturbations at the two sections). Any gain must be attributed to *this additional calibration*, not to discovery of an autonomous material law.

## What is generated

- `results/frozen_transport_comparison.csv`, `frozen_past_lag_selection.csv`, `frozen_calibration_diagnostics.csv`, `frozen_latent_rank_reference_selection.csv`, `frozen_summary.json` and `frozen_transport_comparison.png`: authenticated archived frozen tangent **preflight only**.
- `results/original_*`: the corresponding independently integrated **original Q550 nonlinear server** validations; the original proof is `results/original_ode_transport_proof.json`.
- Original nonlinear server run independently integrates 312 local perturbation/reference trajectories, reusing one trajectory for all three observer projections; **it does not regenerate 400 Q550 reference orbits**.
- All Python scripts, English docstrings/comments, tests, parameterization, pinned RHS references, minimal authentic nested ZIP input, SHA256 manifest and environment information are included. Logs, caches and `.runtime/` are excluded from both manifest and packaged sendback.

## Input genealogy, authentication and limitations

`inputs/v18_24_original_server_sendback.zip` is the **unaltered** original nonlinear v18.24 sendback. `scripts/checkpoint_io_v18_25.py` checks its CRC, root, every SHA256 checksum, the successful original nonlinear proof, and recursively delegates authentication to the older v18.24→v18.23→… checkpoint loaders. The original Q550 RHS is loaded on the server from the *independently pinned* source-tree digests recorded under `reference/q550_state_impulse_response_v13.py`.

* The transport is between **local sections of one reference orbit**. Section-dependent posture and phase are confounded. We do not claim invariance along distinct global orbits.
* The noise stress test depends on the chosen synthetic noise family and norm/coordinate units.
* The first sensor `material3` is a diagnostic projection previously constructed using knowledge of the ODE. The **predictive estimator** is model-blind, not the experiment-design process.
* Rank 20 describes an observed time-delay coordinate representation; it does **not** imply 20 physical material degrees of freedom.
* Calibrated frame transport should be assessed against the zero-shot baseline and negative controls at each noise level. Improving one is not proof of a constitutive material transport rule or uniqueness of the reduced state.
