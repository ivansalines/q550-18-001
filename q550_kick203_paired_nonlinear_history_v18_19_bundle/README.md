# Q550 v18.19 — Paired Nonlinear Trajectory & History Discrimination Audit

**Scientific aim.** Compare distinct histories and futures of *independently integrated nonlinear original Q550 11D ODE trajectories* that share an exactly matched **reduced** material observation at xi=4, with phase held equal. Some pairs additionally match the *first derivative of that observable to linear order* (jet matching). The initial full states intentionally differ. This is an **OFF-CONSTRAINT, model-reduction diagnostic**, not evidence of a constitutive memory law or physical torsional flip.

## Reproduction — single ZIP root verified

Upload `q550_kick203_paired_nonlinear_history_v18_19_bundle.zip` into `/root/elab1/hysteresis`:

```bash
cd /root/elab1/hysteresis
unzip q550_kick203_paired_nonlinear_history_v18_19_bundle.zip
cd q550_kick203_paired_nonlinear_history_v18_19_bundle
bash scripts/run_on_server.sh 2>&1 | tee v18_19_server.log
```

Send back `/root/elab1/hysteresis/q550_kick203_paired_nonlinear_history_v18_19_sendback.zip`.

The original Q550 model is **not replaced**: the script expects the pinned original files (sha256 checked) in `/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle`, and supports changing that root by appending `--base-root /correct/path` to the server run command. Requires scipy/numpy/matplotlib/pandas and Python 3.11 or later. The 400-orbit base is **not regenerated**: original v18.6 kick203 state from authenticated v18.18 nested sendback is reused.

## Procedure

1. Archive gate: authenticate uploaded **original server** v18.18 sendback, follow v18.17→v18.6 genealogy and read frozen J, P, source and kick203 state.
2. Run baseline original nonlinear ODE from xi=0 to xi=4. Check the independently differentiated original RHS Jacobian against archived J at xi=4; reject mismatched model.
3. Define observer H in material 10D with three families: `material3` (span of archived old/material/instantaneous source directions at xi=4), raw component `raw4` (first four), and `raw6` (first six). `material3` is a **diagnostic projection**, not a physically established observer.
4. Construct antisymmetric full-state pairs y±(4)=y_ref(4)±epsilon (h,0), with H h=0 exactly. `jet` seeds additionally satisfy H A(4) h=0 to first order (A=original material Jacobian). Each pair shares identical phase, reduced state, and (approximately to second order) reduced velocity at the meeting point.
5. Integrate the two members **separately** with independent adaptive solver calls, backward to xi=1 and forward to xi=8 with original nonlinear RHS, at eps=[1e-6,3e-6,1e-5] (six directions × three eps = 18 paired experiments). Test nonlinear vs cocycle-derived tangent sensitivity and phase invariance; evaluate history/future observational separation.
6. Include a control where the **entire state is identical**: deterministic future must coincide. No physical/constitutive memory claim is permitted. Any effect observed is compatible with omitted components of a deterministic 11D state.

### A precise interpretation gate

Matching only reduced observables (even their first temporal derivative) at xi=4 while permitting different hidden 10D state components can yield different observed futures. If found, this proves an **empirical failure of Markov closure for that chosen reduction** on the tested pairs. It does *not* prove the model has an additional, independently causal hereditary material variable. If a genuine model-specific constraint C(y) is supplied later, pairs must be projected/retracted on that constraint, and the off-constraint conclusion reassessed.

### Output and fallback

`results/frozen_preflight_summary.json` is generated locally and must always be labeled `FROZEN_PREFLIGHT_ONLY_NO_NONLINEAR_SERVER_RESULT`. Server creates `results/nonlinear_pair_summary.json`, `results/pair_summary.csv`, `results/paired_trajectories.csv`, `results/paired_future_discrimination.png`. Active `v18_19_server.log` and `.runtime/` are **excluded** from SHA manifest to prevent `tee` races; immutable outputs are hashed and packed only after nonlinear run. If it fails, send `v18_19_server.log` and do not label frozen preflight a completed nonlinear result.

The Jacobian is differentiated numerically from the original RHS, independently of the frozen checkpoint; `standard`/`tight` checkpoint tolerance cross-checks are available. Backward nonlinear integration is a numerical diagnostic; if it becomes unstable, result is inconclusive, not an invented successful history comparison.
