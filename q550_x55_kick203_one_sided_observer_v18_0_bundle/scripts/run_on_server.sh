#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
python3.13 scripts/run_kick203_one_sided_observer_v18_0.py \
  --out results/kick203_v18_0 2>&1 | tee v18_0_console.log
python3.13 scripts/pack_sendback_v18_0.py \
  --results results/kick203_v18_0 \
  --out q550_kick203_one_sided_observer_v18_0_sendback.zip
