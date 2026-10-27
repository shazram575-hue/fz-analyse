#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_ROOT"
export PYTHONUTF8=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
python3 --version
df -h .
if [[ "${BIOINFO_INSTALL_PROGRAM_DEPS:-0}" == "1" ]]; then
  python3 -m venv .venv
  .venv/bin/python -m pip install -r requirements.txt
  .venv/bin/python scripts/project/doctor.py --evidence-only
else
  python3 scripts/project/doctor.py --evidence-only
fi
# This setup never downloads large datasets or runs R, NMF, docking, or association tests.
