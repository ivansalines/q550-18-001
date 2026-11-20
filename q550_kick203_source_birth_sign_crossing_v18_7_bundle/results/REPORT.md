# Q550 v18.7 — Source-Birth Sign Crossing & Co-Moving Transport

**Constraint status: OFF_CONSTRAINT_DIAGNOSTIC.** The Euclidean 11D metric is model-coordinate dependent and not a certified physical rigid-body metric.

## Definition

Early source `E=P(s)*a_early`, total quadratic material source `K2=P(s)*a_total`; roots satisfy `E·K2=0`, with `a_early=int_0^1 P^-1 S2 ds` and `a_total=int_0^xi P^-1 S2 ds`.

`d_xi(E·K2) = a_early^T (d_xi G) a_total + a_early^T G (d_xi a_total)`, where `G=P^T P`.

The first term is tangent cocycle deformation (including stretch/shear); the second is newly born source. Neither is a phenomenological memory kernel.

## Crossings and numerical controls

### standard

- Crossings found: 3
- Worst v18.6 11D regression error: 5.391e-14
- Max cocycle condition: 333.4120
- Largest old-source inverse-cocycle drift: 1.745e-15

- Root 1: xi=4.1351989579; material source s=694.168680622; angle old/total=90.000000 deg; signs at ±0.05 xi: 1.17598% / -1.18712%.
  - d(E·K2)/dxi: cocycle -19.6714751, new source -2.66276904, total -22.3342441; finite-difference -22.3342395.
  - Root difference with Hessian half step -6.7311e-09 xi; with coarse Simpson -9.5585e-08 xi.

- Root 2: xi=4.9973365558; material source s=694.284537893; angle old/total=90.000000 deg; signs at ±0.05 xi: -1.99915% / 2.16836%.
  - d(E·K2)/dxi: cocycle +38.303306, new source -0.262167588, total +38.0411384; finite-difference +38.0411341.
  - Root difference with Hessian half step +7.9109e-09 xi; with coarse Simpson +3.5407e-08 xi.

- Root 3: xi=7.1744358136; material source s=694.577104625; angle old/total=90.000000 deg; signs at ±0.05 xi: 0.970125% / -0.977194%.
  - d(E·K2)/dxi: cocycle -71.6762848, new source -0.780365177, total -72.45665; finite-difference -72.4566439.
  - Root difference with Hessian half step -6.2469e-09 xi; with coarse Simpson -6.3028e-08 xi.

### tight

- Crossings found: 3
- Worst v18.6 11D regression error: 5.979e-17
- Max cocycle condition: 333.4120
- Largest old-source inverse-cocycle drift: 1.418e-15

- Root 1: xi=4.1351989579; material source s=694.168680622; angle old/total=90.000000 deg; signs at ±0.05 xi: 1.17598% / -1.18712%.
  - d(E·K2)/dxi: cocycle -19.6714751, new source -2.66276904, total -22.3342441; finite-difference -22.3342395.
  - Root difference with Hessian half step -6.7988e-09 xi; with coarse Simpson -9.5743e-08 xi.

- Root 2: xi=4.9973365558; material source s=694.284537893; angle old/total=90.000000 deg; signs at ±0.05 xi: -1.99915% / 2.16836%.
  - d(E·K2)/dxi: cocycle +38.303306, new source -0.262167591, total +38.0411384; finite-difference +38.0411341.
  - Root difference with Hessian half step +7.9108e-09 xi; with coarse Simpson +3.5520e-08 xi.

- Root 3: xi=7.1744358136; material source s=694.577104625; angle old/total=90.000000 deg; signs at ±0.05 xi: 0.970125% / -0.977194%.
  - d(E·K2)/dxi: cocycle -71.6762848, new source -0.780365225, total -72.45665; finite-difference -72.456644.
  - Root difference with Hessian half step -6.2131e-09 xi; with coarse Simpson -6.3242e-08 xi.

## Tolerance agreement

- Root 1: tight–standard = -2.695089e-11 xi
- Root 2: tight–standard = +2.229417e-11 xi
- Root 3: tight–standard = +3.213874e-11 xi

## Interpretation limit

The early vector is exactly fixed in the inverse fundamental-matrix coordinates up to floating-point error. The polar co-rotating frame removes one canonical Euclidean rotation of the complete cocycle, **not** a mechanically unique body rotation. The crossing of signed 11D projection is a geometric event of the chosen observable, and does not prove a physical torsional flip or hysteresis.

The source data are only sampled every 0.0125 xi. Root interpolation refines the location inside a sampled sign change but does not exclude unresolved pairs of crossings between adjacent nodes. The Simpson and Hessian/tolerance comparisons quantify discretization sensitivity; they are not interval enclosures.
