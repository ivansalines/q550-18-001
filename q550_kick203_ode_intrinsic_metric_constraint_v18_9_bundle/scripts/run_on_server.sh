#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
PYTHON="${PYTHON:-python3.13}"
if ! command -v "$PYTHON" >/dev/null 2>&1;then PYTHON=python3;fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PYTHON" -m unittest discover -s tests -p 'test_*.py' -v
"$PYTHON" scripts/analyze_cached_v18_9.py
if [[ "${Q550_REPLAY:-0}" == "1" ]];then
  "$PYTHON" scripts/model_structure_check_v18_9.py \
    --base-root "${Q550_BASE_ROOT:-/root/elab1/hysteresis/q550_full11D_transverse_slow_state_dense_sampling_v1_bundle}"
else
  echo 'ODE replay optional: set Q550_REPLAY=1 if original pinned model sources survived the server rescue.'
fi
"$PYTHON" scripts/pack_sendback_v18_9.py
