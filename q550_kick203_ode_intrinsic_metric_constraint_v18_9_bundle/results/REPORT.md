# Q550 v18.9 — ODE structure, metric, and constraint audit

Frozen *real* Q550 v18.7/v18.6 full11D source and Jacobian/cocycle arrays; standard and tight independently inspected.

## Observed facts (numerical; not rigorous interval certificates)

- standard: Euclidean crossing xi=[4.135198957, 4.997336555, 7.174435814]; inverse-cocycle-source frame crossing xi=[]; source frame minimal cosine=0.14838696494810638
- tight: Euclidean crossing xi=[4.135198957, 4.997336555, 7.174435814]; inverse-cocycle-source frame crossing xi=[]; source frame minimal cosine=0.14838696493410847

## Coordinate-free statement

For linearly independent vectors E and K the sign of E^T W K is not intrinsic without choosing a positive-definite metric W. This is mathematical, not a numerical defect. Co-moving pullback W=P^-T P^-1 only preserves an arbitrary initial metric.

## Tests of the original ODE tangent structure

Each archived J gives a local linearization; the audit checks P evolution, the transport-vs-source derivative decomposition, a necessary nullspace screen for a constant linear first integral, and a constant-symmetric-quadratic-form screen. These screens cannot establish or refute nonlinear conserved quantities or a state-dependent material metric.

## Constraint gate

No genuine C(y) or material energy G(y) accompanied the archived checkpoint. The v18.6 source was generated from off-constraint kicks. Even if C(y) is found, projecting K2 onto a tangent plane would incorrectly erase required 1/2 D2C[u,u] normal curvature.

## Optional original-ODE replay and candidate verification

Use `python3 scripts/model_structure_check_v18_9.py --base-root /path/to/original/q550/base` to replay from the **archived kick state** without repeating 400 orbits. To inspect independently sourced constraints and metric also pass `--structure-module /path/to/candidate.py --structure-provenance /path/to/metadata.json`. Report labels remain candidate diagnostics pending external physical verification.

## Reproduction

`bash scripts/run_on_server.sh`; no external Q550 files needed for mandatory frozen-checkpoint audit; full original v18.7→v18.6→v18.5 genealogy retained under inputs.

## Important limitation

The computed crossing count is only grid-bracketed (xi step 0.0125), not an exclusion proof for missed crossing pairs. Absence of roots in one chosen metric is NOT absence of memory, torsion, or another physical mechanism.
