#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PYTHON="${PYTHON:-python3.13}"
"$PYTHON" -m unittest discover -s tests -v
"$PYTHON" scripts/preflight_v18_6.py
"$PYTHON" scripts/run_source_transport_v18_6.py "$@" 2>&1 | tee v18_6_solver.log
"$PYTHON" scripts/analyze_source_transport_v18_6.py 2>&1 | tee v18_6_analysis.log
"$PYTHON" scripts/pack_sendback_v18_6.py
