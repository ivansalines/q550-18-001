#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PYTHON="${PYTHON:-python3.13}"
"$PYTHON" -m unittest discover -s tests -v
"$PYTHON" scripts/preflight_v18_4.py
# v18.3 anchor levels are consumed as read-only inputs; only new amplitudes run.
"$PYTHON" scripts/run_scan_v18_4.py "$@" 2>&1 | tee v18_4_solver.log
"$PYTHON" scripts/analyze_scan_v18_4.py 2>&1 | tee v18_4_analysis.log
"$PYTHON" scripts/pack_sendback_v18_4.py
