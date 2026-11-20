# Nonlinear constraint / metric provenance gate

At v18.11, no material C(z,theta), DC, D2C, or physical G has been established.
No Hessian of the kinematic phase leaf may be substituted for a material constraint Hessian.

A future nonlinear constraint audit must provide independent evidence of:

1. The original constraint C, state-coordinate definition, domain and its derivation.
2. DC and D2C validated by finite differences and constraint preservation DC f=0 (or explicitly time-dependent partial_t C + DC f=0).
3. Consistent second-order perturbation: DC w + D2C[u,u]=0 (with the declared coefficient convention).
4. Explicit physical provenance and nonsingularity of metric G, including whether its comparison is under pushed-forward metric for coordinate changes.
5. Interval or independently convergent bounds before using the word certificate.

A phase clock is NOT a conservation law or holonomic material constraint.
