# Q550 kick203 — v18.13 Transverse Source Reorientation & Radial Turning Audit

**Scope.** Determine how the co-moving material-memory plane changes direction on either side of the four radial-rate turnings identified in v18.12. This is a numerical audit in the **declared initial Euclidean 10D frame**; it is not a physical torsion or hysteresis certificate.

## Frozen genealogical source

- `inputs/v18_12_original_sendback.zip` is the actual v18.12 sendback. It preserves the nested verified previous source chain (v18.11, v18.10, v18.9, v18.7 and v18.6). CRC and the v18.12, v18.11, v18.10, v18.9 SHA256 manifests are verified before analysis.
- `reference/` holds executable v18.9 lineage code; `scripts/checkpoint_io_v18_11.py` and `scripts/residual_plane_math_v18_12.py` preserve verified original methods.
- `scripts/checkpoint_io_v18_13.py` validates the v18.12 archive, extracts the actual nested v18.11 input, and invokes the historical checkpoint loader.
- `scripts/reorientation_math.py`, `scripts/analyze_v18_13.py`, and `scripts/plot_v18_13.py` are actual executable sources for the NEW audit. `tests/test_geometry.py` contains seven independent sanity/negative tests.
- `results/` holds the local frozen checkpoint audit output. `SHA256SUMS_OUTPUT.json` records hashes of the bundle contents. `scripts/pack_sendback.py` refreshes it during execution and writes a verified sendback ZIP to the parent directory.

## Execute on original server

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_transverse_source_reorientation_v18_13_bundle.zip
cd q550_kick203_transverse_source_reorientation_v18_13_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_13_server.log
```

Send back: `/root/elab1/hysteresis/q550_kick203_transverse_source_reorientation_v18_13_sendback.zip`.

No 400-orbit replay and no fresh full-ODE integration. The previously independent ODE replay results are carried by the nested genealogy.

## Mathematical protocol

Let `e` be the pre-existing source at `xi=1`, `A(xi)` the cumulative source in the co-moving material frame, `g(xi)` its source integrand, and `s = ds_ODE/dxi`. In the declared source-frame Euclidean geometry set:

    ehat = e / ||e||
    w = A - (ehat·A) ehat
    rho = ||w||
    v = w / rho
    gperp = g - (ehat·g) ehat
    rho' = s (v·gperp)
    v' = s (gperp - (v·gperp)v)/rho

The radial turnings are sign changes of `rho'`; locate them using linear, PCHIP and natural cubic interpolation, and verify downsampled parity grids. The tangent direction `v'/||v'||` on either side of an event is compared using shortest-geodesic Levi-Civita parallel transport on the unit sphere in the 10D Euclidean source frame. Test windows `±0.05`, `±0.10`, `±0.25`, `±0.50` xi. Also vary the Euclidean gauge by explicit fixed diagonal coordinate rescalings.

**Important identity:** at a radial turning, `rho' = 0`, the source has zero radial component by definition. A near-one transverse source fraction at that event is not independent evidence for torsion or hysteresis.

## Observed results (original v18.12 frozen checkpoints)

- Four radial events near xi `2.62565, 3.94961, 5.56011, 6.95829` in the initial Euclidean source frame.
- At `±0.05 xi`, parallel-transported directional cosines are `0.98875, 0.96467, 0.99153, 0.97576`.
- At `±0.50 xi`, they are `0.02665, -0.98998, 0.11877, -0.83942`. Thus immediate directional continuity can coexist with substantial longer-scale reorientation.
- Standard/tight maximum radial root difference is ~`2.9e-10 xi`; interpolation methods disagree by at most `~7.5e-5 xi`, indicating interpolation uncertainty exceeds the solver-tolerance spread.
- The radial-turning locations **shift under anisotropic coordinate reweighting**. They are not coordinate-free material events.

See `results/REPORT.md` and CSVs for exact numbers, test definitions, and the metric sensitivity table.

## Evidence boundary

- Linear dependence/independence and vanishing/nonvanishing wedge products are invariant under invertible changes of coordinates, within numerical rank resolution. Normalized volumes, radial turnings, unit-sphere geodesic coherence, lengths and cosines are **not** invariant under arbitrary GL(10) coordinate changes.
- A nonvanishing residual transverse source, even with local continuity, does **not** establish a physical body torsion, mechanical constraints, non-holonomic effects, irreversibility, or history-dependent hysteresis. Distinct histories at matched phase-fixed material states were not compared.
- Root refinements are numerical interpolation checks, **not** interval certificates or independent original ODE replays.
