# Q550 v18.15 — Age-Slab Source Interference & Echo Localization Audit

## Test definition

At each pre-existing four radial event times, the unchanged cumulative source A(xi) in the initial co-moving 10D representation is decomposed by **age** a=xi-u, most recent first:

S_i=A(xi-a_i)-A(xi-a_{i+1}); w_i=Q_e S_i; Q_e=I-ee^T/(e.e), e=A(1).
The sum of all w_i equals w=Q_e A(xi). For v=w/||w||, the signed axial contribution is L_i=v.w_i, the lateral contribution is w_i-L_i*v, and the pairwise interaction is w_i.w_j.
By construction sum L_i=||w|| and sum lateral=0; sign, angles and interference magnitude are Euclidean-gauge dependent.

## Event summary (standard tolerance)

| Event | xi | min/max | old versus recent cosine (h=2) | cancellation (h=.25) | cancellation (h=.0625) | 10-degree stable age horizon (h=.0625) |
|---|---:|---|---:|---:|---:|---:|
| 1 | 2.625653 | radial_max | -0.501 | 0.053 | 0.053 | 1.250 |
| 2 | 3.949610 | radial_min | -0.533 | 0.069 | 0.069 | 2.688 |
| 3 | 5.560115 | radial_max | 0.151 | 0.432 | 0.435 | 3.875 |
| 4 | 6.958293 | radial_min | -0.726 | 0.336 | 0.347 | 5.562 |

## Controls

Max standard/tight event drift: 2.841e-10 xi.
The CSV `resolution_convergence.csv` measures signed-negative mass, cancellation, angular stabilization and independent source-integrand quadrature for slab sizes 0.5, 0.25, 0.125, 0.0625 xi.
The CSV `strongest_negative_interference_pairs.csv` identifies pairs with negative dot products. These are gauge-specific interference diagnostics, not forces or material energy flows.

## Interpretation limits

This is **not** an ablation of physical memory in the ODE. A different source-history decomposition from the same fixed trajectory is a diagnostic of the observer. Proof of physical hysteresis requires distinct physically realizable histories meeting at matched state and phase.
Age axes refer to the original co-moving xi parameter, not directly to laboratory time.

## Source genealogy

Embedded, authenticated original v18.14 sendback SHA256: `6e3f5b31ee2b42c6d40c3933613cb8e0103fa3fd41f87e9629349126cf6be25a`.
All original v18.14 and nested v18.13 -> v18.6 generation scripts remain inside the preserved original sendback archive.
No 400-orbit recomputation.
