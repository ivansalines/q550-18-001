#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
PYTHON="${PYTHON:-python3.13}"
command -v "$PYTHON" >/dev/null 2>&1 || PYTHON=python3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
code=0
"$PYTHON" -m unittest discover -s tests -p 'test_*.py' -v
"$PYTHON" scripts/independent_frozen_audit.py
if [[ "${Q550_SKIP_NONLINEAR:-0}" != "1" ]]; then
    if ! "$PYTHON" scripts/nonlinear_phase_probe.py --base-root "${Q550_BASE_ROOT:-/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle}"; then
        echo "NONLINEAR ODE AUDIT FAILED; packaging diagnostic anyway"
        code=1
    fi
else
    echo 'NONLINEAR ORIGINAL ODE PROBE SKIPPED — frozen checkpoint only'
fi
"$PYTHON" scripts/pack_sendback.py
exit "$code"
