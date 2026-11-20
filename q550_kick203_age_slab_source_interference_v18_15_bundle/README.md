# Q550 v18.15 — Age-Slab Source Interference & Echo Localization Audit

This is the original Q550 source-history audit (kick 203), using the **full authenticated v18.14 server sendback** embedded unchanged as `inputs/v18_14_original_sendback.zip`. It recursively reuses original v18.7/6 numerical source inputs and the original reader scripts, validating SHA256 and ZIP CRC. It does **not** rerun the 400-orbit baseline.

## Exact server commands (one directory root inside ZIP)

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_age_slab_source_interference_v18_15_bundle.zip
cd q550_kick203_age_slab_source_interference_v18_15_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_15_server.log
```

You will find **`/root/elab1/hysteresis/q550_kick203_age_slab_source_interference_v18_15_sendback.zip`**. Send that ZIP back for review.

The bundle archive itself has exactly one top-level directory named `q550_kick203_age_slab_source_interference_v18_15_bundle/`. The absolute paths shown above have been tested from a fresh extraction in a clean directory; `run_on_server.sh` internally locates its root with a path derived from the script itself, not from the caller's working directory. Existing directories must not be overwritten accidentally; move or remove a previous extraction if rerunning in place.

## Mathematical test

Work in the **initial co-moving Euclidean 10D reference**, where `A(xi)` is the accumulated historical source and `s*g(xi)` its archived integrand. The early-memory direction is fixed at `e=A(1)`. For each of the four *previously observed* radial events (v18.13/v18.14), cut the history by source **age** a=xi-u, from most recent to oldest, into widths `.5`, `.25`, `.125`, `.0625` in xi units.

For each age slab `[a_i,a_{i+1}]`:

```
S_i = A(xi-a_i)-A(xi-a_{i+1})
Q_e = I - e e^T/(e^T e)
q_i = Q_e S_i
w = Q_e A(xi); v = w/||w||
L_i = q_i . v                           (signed axial support)
T_i = q_i - L_i v                       (lateral support)
I_ij = q_i . q_j                       (signed pairwise interference)
```

Exactly, `sum q_i=w`, `sum L_i=||w||`, `sum T_i=0`. Negative `L_i` identifies slabs opposing the **final direction**, not negative energy. The full source integral is unchanged; source terms are diagnostically decomposed, **not deleted from or fed back into the ODE**.

`destructive_fraction = 1-||w||/(sum |L_i|)` is a gauge- and partition-dependent indicator of cancellation: it is *not* a percentage of material energy or physical hysteresis. The age localization includes each slab's leave-one-out line-angle effect, the strongest negative pairwise products, old/recent h=2 interference, and the earliest age past which all larger age prefixes remain within 10°/25° of the full direction.

Independent consistency controls include:
- Original v18.14 event-location regression, plus `standard`/`tight` agreement.
- Closure of all slab vectors, signed projections, and pairwise norm identity.
- Age-width convergence across four resolutions, **including cancellation variation**, not merely algebraic closure.
- Independent integration of interpolated original source integrand `g` versus differences of interpolated A. Note these two interpolants are distinct and need not coincide at machine precision.
- Synthetic opposing/positive source slabs, coordinate-rank invariance, rejection of degenerate history, path/layout checks.

**Limits**: the physical material metric and mechanical constraint manifold have not been established. Angles, dot products, interference and event positions depend on the chosen Euclidean gauge. This is an observer-history diagnostic rather than a material hysteresis certificate. There are only four events and no universal alternation theorem.

## Contents

- `inputs/v18_14_original_sendback.zip` — authenticated original with full nested script genealogy.
- `scripts/checkpoint_io_v18_15.py` and original `checkpoint_io_v18_14/13/11.py` — authenticated input traversal.
- `scripts/age_slab_math.py` — mathematical decomposition and signed interference.
- `scripts/analyze_v18_15.py`, `scripts/plot_v18_15.py`, `scripts/pack_sendback.py` — executables that generate this audit.
- `tests/test_age_slabs.py` — independent synthetic and structural tests.
- `results/` — computed `summary.json`, `REPORT.md`, event and slab CSVs, refinement tables and images.
- `reference/` — original v18.9 reference routines.
- `requirements.txt`, `ENVIRONMENT.json`, `SHA256SUMS_OUTPUT.json` — environment and hashes.

Python 3.10+ and numpy/scipy/matplotlib.
