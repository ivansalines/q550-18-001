#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PYTHON="${PYTHON:-python3.13}"
"$PYTHON" -m unittest discover -s tests -v
"$PYTHON" scripts/run_kick203_material_transport_v18_2.py \
  --out "$ROOT/results/kick203_v18_2" "$@" \
  2>&1 | tee "$ROOT/v18_2_console.log"
"$PYTHON" scripts/pack_sendback_v18_2.py \
  --out "$ROOT/q550_kick203_on_shell_full11d_material_transport_v18_2_sendback.zip"
