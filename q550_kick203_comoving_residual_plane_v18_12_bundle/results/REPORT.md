# Q550 v18.12 — Co-moving residual plane after common transport

## Result

- Source verified through v18.11 -> v18.10 -> v18.9 -> v18.7 -> v18.6, SHA256.
- Common homogeneous transport removed by the inverse 10x10 fixed-phase propagator.
- Both old and quadratic source identities reconstructed at 1e-8 relative tolerance.
- Source-frame residual plane turn rate median **1.11720608** per xi, min **0.610229298**, max **2.49843313** (Euclidean initial frame ONLY).
- Positive turn rates at **560/560** sampled nodes beyond xi=1.
- Integrated turn path **9.24311805** rad versus endpoint angle **0.875404879** rad. The discrepancy indicates bending/backtracking in the declared gauge, not physical hysteresis.
- Max condition number for Phi **330.958**; maximum 10D inverse transport errors E/K **1.89e-15/1.35e-15**.
- Independent finite differences of v(xi) agree with sourced analytical derivative at median relative error **0.000189**.
- GL(10) covariance of source 2-form: relative error **2.2e-16**.

## Mechanism

Let e=A(1), K=Phi A, E=Phi e, A_xi=s*g and define ehat=e/||e||, w=(I-ehat ehat^T)A, v=w/||w||.

    v_xi = s * (I-ehat ehat^T-v v^T) g / ||w||.

The homogeneous J transport is no longer present. A nonzero residual means fresh source drives the plane outside its previously occupied 2-plane, **within the chosen fixed initial-frame gauge**.

**Strict boundary**: a changing accumulated source A necessarily remembers earlier integration, but is not automatically hysteresis (multi-valued response to history, loop memory, irreversibility) and is not a certificate of physical body torsion. No source-derived physical metric G(y) or constraints C(y) were provided. The path length is metric dependent, although nonzero raw 3-blade e wedge A wedge g is algebraically invariant under invertible reparametrization.

## Next experimental gate

Compare histories that reach the same verified reduced material state/phase with different kicks, then measure unequal source-frame residual planes/response; alternatively construct a physical energy metric and constraint manifold from the real ODE. Do NOT infer hysteresis merely from cumulative source turning.
