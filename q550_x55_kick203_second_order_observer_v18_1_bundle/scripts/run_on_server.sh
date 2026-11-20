#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
python3.13 -m unittest discover -s tests -v
python3.13 scripts/run_kick203_second_order_observer_v18_1.py \
  --out "$ROOT/results/kick203_v18_1" 2>&1 | tee "$ROOT/v18_1_console.log"
python3.13 scripts/pack_sendback_v18_1.py \
  --results "$ROOT/results/kick203_v18_1" \
  --out "$ROOT/q550_kick203_second_order_observer_v18_1_sendback.zip"
