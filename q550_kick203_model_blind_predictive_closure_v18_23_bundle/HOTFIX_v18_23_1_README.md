# Q550 v18.23.1 — Latent rank validation-tie hotfix

## Symptom
`AssertionError: 20 != 10` from `test_latent_reconstructs_unseen_linear_futures`.

## Diagnosis
The synthetic dataset has a 10D generator, but the learned latent delay
rank is a validation-selected predictor parameter. On different LAPACK/BLAS
builds, near-zero validation errors at rank >=10 may change order at floating
point precision, causing rank 20 to win an arbitrary near-tie. This does NOT
establish 20 physical dimensions or invalidate the numerical experiment.

## Correction
- `blind_math.py` now chooses the smallest rank within a validation-error
  numerical-tie threshold, 1e-12*max(1,RMS(validation targets)).
- The rule uses validation data ONLY and explicitly retains rank 20 whenever
  its error is meaningfully smaller than a simpler candidate.
- The synthetic test checks both exact noiseless reconstruction at rank 10
  and the parsimonious numerical-tie choice. Three extra tests cover tie
  handling, meaningful rank-20 improvement, and no test-set selection.
- Output latent-rank selection metadata records the tolerance and rule.
- No physics or source data and no server Q550 ODE were modified.

## In-place application, starting from /root/elab1/hysteresis
```
unzip -o q550_v18_23_1_latent_rank_tie_hotfix.zip
cd q550_kick203_model_blind_predictive_closure_v18_23_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_23_server_retry.log
```

The first unzip is run in the parent directory. Do NOT unzip a second copy of
v18.23 itself or delete the original input/trace. The normal run script will
recompute the immutable-output manifest and prepare the sendback.

The Q550 server has not yet performed its nonlinear validation in this run;
unit tests and frozen preflight alone are not server certification.
