# Optional model-structure contract (NOT Q550 evidence)

Use **independently documented** formulas in the **same full 11D coordinates and units** as the frozen Q550 source state. Do not call ordinary diagonal coordinate rescalings "mechanical metrics". The script will not treat external callbacks as proven physical laws.

Example structure module **interface only** (do not put invented example physics into your Q550 results):

```python
# structure_candidate.py
import numpy as np

def metric(y: np.ndarray) -> np.ndarray:
    # REQUIRED only if testing a candidate material metric:
    # return verified 11x11 symmetric positive-definite tensor in original y-units
    raise NotImplementedError

def constraints(y: np.ndarray) -> np.ndarray:
    # OPTIONAL; return m verified constraint functions C(y) (not just a Jacobian)
    raise NotImplementedError

def constraint_jacobian(y: np.ndarray) -> np.ndarray:
    # OPTIONAL; otherwise centered numerical differentiation of C is used
    raise NotImplementedError
```

Only define the functions actually implemented. Functions that raise `NotImplementedError` should be removed.

Required provenance JSON keys: `name`, `physical_reference`, `coordinate_mapping`, `metric_or_constraint_formula`, `derivation_method`. Each must be a nonempty string pointing to a **verifiable** Q550 source. The audit computes algebraic consistency, checks `C(y)` drift and `DC f`, and considers

- first variation `DC u`, with original `u=P u0`;
- second variation `DC K2 + 0.5 D2C[u,u]`, appropriate when `K2` is the coefficient of epsilon squared in `y(epsilon)=y0+epsilon u+epsilon^2 K2+...`.

Because the historical kick was **off-constraint**, a nonzero residual is expected for many true constraints; no constrained physical conclusion follows. A successful diagnostic does not prove physical provenance. A true on-constraint experiment requires retracted initial conditions and comparison under actual ODE integration.

To run after base ODE sources are available:

```
python3 scripts/model_structure_check_v18_9.py \
  --base-root /root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle \
  --structure-module /absolute/path/to/structure_candidate.py \
  --structure-provenance /absolute/path/to/structure_provenance.json
```

The native Q550 audit module, support RHS and support Jacobian sources are checked against historical SHA256 before replay; a changed model is a **hard failure**. If these files are unavailable after a server crash, the core archived audit still runs and reports its limited conclusions.
