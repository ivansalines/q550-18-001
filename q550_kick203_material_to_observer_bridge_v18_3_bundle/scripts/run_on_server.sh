#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PYTHON="${PYTHON:-python3.13}"
echo "v18.3 postprocessing-only bridge; original Q550 orbit is NOT replayed"
"$PYTHON" -m unittest discover -s tests -v
"$PYTHON" scripts/run_material_to_observer_bridge_v18_3.py "$@" \
  2>&1 | tee results/kick203_v18_3/console.log
"$PYTHON" scripts/pack_sendback_v18_3.py
