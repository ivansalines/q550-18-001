#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
PYTHON=${PYTHON:-python3.13}
if ! command -v "$PYTHON" >/dev/null 2>&1; then PYTHON=python3;fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PYTHON" -m unittest discover -s tests -p 'test_*.py' -v
"$PYTHON" scripts/preflight_v18_21.py
"$PYTHON" scripts/nonlinear_server_v18_21.py "$@"
"$PYTHON" scripts/plot_v18_21.py
"$PYTHON" scripts/build_manifest.py
"$PYTHON" scripts/pack_sendback.py
