# Q550 kick203 v18.1 — one-sided second-order observer convergence

B crossing: `xi_B = 210.861136603177698`. Same frozen A=hquarter and full11D M in all comparisons.

B observer derivatives (t is the 11D unit tangent):

- Forward first: `(t(s+h)-t(s))/h`
- Backward first: `(t(s)-t(s-h))/h`
- Forward second: `(-3t(s)+4t(s+h)-t(s+2h))/(2h)`
- Backward second: `(3t(s)-4t(s-h)+t(s-2h))/(2h)`

Only backward stencils are used before the event and only forward stencils after the event.
Four B steps: historical hhalf and hquarter, and extensions heighth and hsixteenth.

Tests: original central frames reproduce v17.9, and the first-order hquarter residuals reproduce v18.0. Shapley is algebraically checked in full11D.

Max regression error: 0.
Max Shapley closure: 5.74e-17.

No physical torsion/flip certification: off-shell allocation and observer sensitivity remain.
See second_order_full11D_arrays.npz, observer_orders.csv, residual_convergence.csv and figures.
