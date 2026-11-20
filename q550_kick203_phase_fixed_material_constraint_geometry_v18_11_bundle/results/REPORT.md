# Q550 v18.11 — Phase-fixed material transport & constraint geometry

## Verified source

Nested v18.10 SHA256 c4cf85493a693755349f1b6d48e8fd8123739eb3e4ed8f5174c242fa8d37fb50. SHA provenance checked through historical v18.6.

## Mathematical structure

State y=(z,theta), theta_dot=2pi along original sampled trajectory, J=[[A,b],[0,0]]; the tangent bundle dtheta=0 is preserved by the linear cocycle P=[[Phi,q],[0,1]].
The 10D model is z_dot=f_z(z,theta0+2pi*(s-s0)), GENERALLY NONAUTONOMOUS. q_dot=A q+b is the effect of changing the initial phase. Do not identify the phase time-leaf with a material constraint.
The algebraic rank-2 plane span(E,K) is preserved as 10D-vectors, but its Euclidean projector is not GL-invariant. Projector speed decomposes into homogeneous transport and fresh-source birth. No metric-independent angular rate is claimed.

## Numbers

- standard: phase structural max defect 0.000e+00; maximal phase-to-z Jacobian column norm 15.2826.
  Median 10D Phi tangent rate finite-difference residual 1.306e-04; driven q rate median 9.020e-05.
  Median measured Grassmann total/transport/birth speed 1.80077, 1.91762, 0.340181 (Euclidean gauge).
  Median projector rate finite difference relative error for xi>=1.1: 4.714e-04.
- tight: phase structural max defect 0.000e+00; maximal phase-to-z Jacobian column norm 15.2826.
  Median 10D Phi tangent rate finite-difference residual 1.306e-04; driven q rate median 9.020e-05.
  Median measured Grassmann total/transport/birth speed 1.80077, 1.91762, 0.340181 (Euclidean gauge).
  Median projector rate finite difference relative error for xi>=1.1: 4.714e-04.

## Structural gate

Original-ODE server test re-integrates the 10D phase-driven variational block independently of the full 11D state and compares it to archived P; source modules are pinned with SHA256. This is a numerical model-consistency audit, not an invariant-manifold or physical-metric certificate.
No independent nonlinear material constraint C(y), D2C, or mechanical G(y) is supplied, so a second-order normal-curvature correction is neither guessed nor claimed. The ODE phase leaf C(s,y)=theta-theta0-2pi*(s-s0) is a kinematic clock identity only.
Every plot has its own declared 10D Euclidean gauge. Invariants: exact linear rank and the invariant phase-fixed variational subbundle relative to the supplied phase split (not arbitrary phase-mixing coordinate transforms).

