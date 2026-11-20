#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PYTHON="${PYTHON:-python3.13}"
"$PYTHON" -m unittest discover -s tests -v
"$PYTHON" scripts/analyze_metric_v18_8.py "$@" 2>&1 | tee v18_8_analysis.log
"$PYTHON" scripts/pack_sendback_v18_8.py
