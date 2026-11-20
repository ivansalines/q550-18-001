# Q550 v18.8 — Material Metric & Constraint Sensitivity Audit

**Status: FULLY EXECUTED postprocessing on actual Q550 v18.7 full11D checkpoints.**
No new 400-orbit Q550 replay required. All stored source/cocycle data are the
original frozen Q550 model outputs, not synthetic results. Unit tests use
separate synthetic data and do NOT assert facts about Q550.

## Scientific question

Do the three v18.7 sign crossings at xi = 4.135199, 4.997337, 7.174436
persist under changes of the metric used to compare early material source
`E=P a_old` with accumulated quadratic source `K2=P a_total`?

Observable for each SPD state metric W(xi): `N_W(xi)=E(xi)^T W(xi) K2(xi)`.
Sign changes depend on W unless the metric has a separately established
mechanical/physical interpretation.

## Families

- `W=I`: original 11D Euclidean diagnostic, reproduces v18.7 roots.
- `W_alpha=(1-alpha)I + alpha P^-T P^-1`, alpha in
  `[0,.05,.10,.15,.20,.25,.30,.40,.50,.75,1]`:
  an SPD comparison between laboratory Euclidean and source-frame transport
  pullback (the latter is **time dependent, not a mechanical energy metric**).
- `W_ref=P(xi_ref)^-T P(xi_ref)^-1` for xi_ref=1,2,4,8: constant SPD
  metrics induced by freezing a cocycle frame; **diagnostics only**.
- 200 deterministic random diagonal SPD reweightings per maximum weight-ratio
  bound `[1.1,1.5,2,4,10]` and per tolerance (`seed=184`), to detect
  coordinate unit sensitivities. This ensemble is not a distribution over
  physically plausible metrics.
- Eleven leave-one-channel-downweightings with weights 0.01 or 0.1,
  explicitly **not true constraints**.
- Optional, externally supplied constraint Jacobians, see
  `templates/CONSTRAINT_INPUT.md`. Unless supplied and independently verified,
  no on-constraint status is claimed.

Every dot derivative for a FIXED metric is evaluated from original `J` and
source-integrand (not a numerical differencing of a sign chart). Pullback
mixtures have analytic derivatives that account for the metric variation. Each
root is bracketed on the archived xi grid (`0.0125` spacing) and refined with
cubic Hermite interpolation. There is no interval exclusion of invisible root
pairs inside a grid interval.

## Run

```bash
unzip q550_kick203_material_metric_constraint_sensitivity_v18_8_bundle.zip
cd q550_kick203_material_metric_constraint_sensitivity_v18_8_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_8_server.log
```

Python 3.13 is the default; set `PYTHON=python3` to override. No external Q550
sources are needed because all necessary 11D arrays and the historical
v18.7→v18.6→v18.5→v18.4 executable genealogy are archived.

Output: `q550_kick203_material_metric_constraint_sensitivity_v18_8_sendback.zip`.
For local re-analysis: `python3 scripts/analyze_metric_v18_8.py`.

## Included

- `inputs/v18_7_original_sendback.zip`: **original user v18.7 sendback**, itself
  including the complete v18.6, v18.5 history, original Q550 Python scripts,
  inputs, reports, manifestations and full 11D checkpoints.
- `scripts/*.py`: complete v18.8 computation and packing script; `tests/` tests;
  `templates/CONSTRAINT_INPUT.md`: explicit physically sourced constraint API.
- `results/`: computed CSV, JSON, plots, report; `SHA256SUMS.json`: manifest;
  pinned `requirements.txt`.

## Interpretation limit

This bundle determines **metric sensitivity** of a nonlinear transported 11D
observable. It does *not* by itself find a material physical metric or certify
conservation, on-manifold perturbations, a rigid-body torsional flip, or
constitutive hysteresis. Real constraint invariants for this specific support
ODE are not supplied by the available archived v18.7 model inputs.

## Actual Q550 results (computed from frozen v18.7 checkpoints)

- Reconstructed original Euclidean roots in standard mode:
  `4.135198957200635`, `4.997336555349751`, `7.174435813608548`.
  Maximum absolute discrepancy from v18.7 server-refined roots: `7.16e-10 xi`.
- Pullback mixture: the first pair of roots merges/disappears at the
  approximate diagnostic parameter `alpha=0.24646845`; the last root
  disappears near `alpha=0.37911434`. Agreement of standard and tight
  thresholds is approximately `3e-11` in alpha. These are NOT model thresholds.
- Random bounded diagonal weight-ratio 1.1, 1.5, or 2: three roots in all
  200 deterministic tests per weight bound; maximum bound 10:
  161 three-root, 32 one-root, five five-root, and two seven-root cases.
- Fixed whitening at xi=1,2,4,8: respectively 2,4,3,4 resolved roots.
- Downweighting channel index 6 to 0.01 yields five roots. Coordinate
  suppression is NOT a conserved-quantity or constraint projector.
- No verifiable Q550 constraint equations or physical metric were delivered
  with the original v18.7 files; consequently the constrained physical
  sub-experiment remains explicitly **not performed**.

Results and their caveats are included in `results/REPORT.md` and
`results/summary.json`. All sources/comments/reports inside the artifact are in English.
