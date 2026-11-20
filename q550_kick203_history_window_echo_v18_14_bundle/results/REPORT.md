# Q550 v18.14 — History-Window Ablation & Transverse Source Echo Audit

## Premise

The exact source history A(xi) and integrand g(xi) are retained from the verified Q550 v18.13 genealogy; standard and tight checkpoints verified recursively. No 400-orbit replay.

For h=0.25,0.5,1,2: A_recent=A(xi)-A(xi-h); A_old=A(xi-h). At fixed xi A_total=A_recent+A_old.
Under a moving window A_recent prime = s*(g(xi)-g(xi-h)); A_old prime=s*g(xi-h).
All orientation and scalar comparisons occur in the *declared initial Euclidean 10D frame*.

## Event-centred echo — angular distance from full history

| event | xi | h=.25 | h=.50 | h=1.0 | h=2.0 |
|---|---:|---:|---:|---:|---:|
| 1 | 2.6256531 | 72.66 deg | 54.27 deg | 23.92 deg | 3.07 deg |
| 2 | 3.9496098 | 83.98 deg | 78.94 deg | 75.19 deg | 63.20 deg |
| 3 | 5.5601149 | 77.87 deg | 64.92 deg | 35.67 deg | 20.04 deg |
| 4 | 6.9582933 | 81.08 deg | 74.53 deg | 67.87 deg | 69.22 deg |

## Independent old-history direction outside the recent-history plane

The vector triple [e, A_recent, A_old] has numerical rank 3 for all four event times and all four windows, using SVD of separately normalized columns. This is a coordinate-invariant *rank* property, not a coordinate-invariant SVD magnitude.
Smallest event rank-3 singular value in declared gauge: 0.465830.
Additive history and moving-window derivative identity residual: 1.64e-14.
Max standard/tight event-angle discrepancy: 1.942e-08 deg.

## Scope and interpretation

Old-history contributions have distinct orientations at all four events. In this gauge the minima (events 2 and 4) retain a large angle to the complete-history line even with h=2. This is *history support of a cumulative observable*; it is not dynamical hysteresis by itself.
Ablation of A at fixed physical trajectory does NOT model deletion of constitutive memory, because no source perturbation was fed back into the ODE. For physical hysteresis, compare two distinct admissible histories meeting at a common phase-fixed material state and then measure their responses.

## Provenance

Original v18.13 SHA256: 3968e63515e98b60616ad710d3472932b7f0b1134e5245d464be1f8124a59198
All scripts, prior original sendback, summary, profiles, event audit, plots, tests and SHA256 manifest retained.
