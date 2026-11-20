# Q550 kick203 centered/one-sided B-observer audit v18.0

B crossing `xi_B = 210.861136603177698`. Selected `47` of original `320`
points. No re-integration of the saved 11D variational cocycle M.

Observer differences (t is the unit 11D state tangent):

- centered: `d=(t(s+h)-t(s-h))/(2h)`;
- forward: `d=(t(s+h)-t(s))/h`;
- backward: `d=(t(s)-t(s-h))/h`;
- non-straddling: backward strictly before B, forward strictly after B.

Only the **B frame** is modified. A and M remain the original v17.9
arrays in all allocations. Same +/- nonlinear histories, same perturbations,
same Q550 parameters, same source provenance.

Expected centered support halfwidth in xi:
- h/2: `0.0069426298796`
- h/4: `0.0034713149398`

Strict reconstruction checks:
- max central B frame difference: `0.000e+00`
- max centered Shapley vector difference: `0.000e+00`
- max algebraic closure error: `6.262e-17`

## Boundaries

The archived v17.9 M and A frames are frozen **by design**, so any effect in
this audit is attributable to the B-observer construction, not to a changed
transport. The Shapley hybrid histories remain off-shell. This is not a
certification of physical flip or a force. The exact event node is undefined
for the non-straddling observer; one-sided offsets are compared instead.

See `observer_comparison.csv`, `symmetric_offset_comparison.csv`, figures,
`summary.json` and `central_regression.csv`.
