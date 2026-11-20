#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
PYTHON="${PYTHON:-python3.13}"
if ! command -v "$PYTHON" >/dev/null 2>&1; then PYTHON=python3;fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PYTHON" -m unittest discover -s tests -p 'test_*.py' -v
"$PYTHON" scripts/analyze_plane_v18_10.py
"$PYTHON" scripts/plot_plane_v18_10.py
# Run pinned original model replay automatically on Q550 server.
if [[ "${Q550_SKIP_REPLAY:-0}" != "1" ]];then
    "$PYTHON" scripts/inspect_original_ode_v18_10.py \
      --base-root "${Q550_BASE_ROOT:-/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle}"
else
    echo 'Original ODE replay explicitly skipped by Q550_SKIP_REPLAY=1 (results remain a frozen-source audit).'
fi
"$PYTHON" scripts/pack_sendback_v18_10.py
