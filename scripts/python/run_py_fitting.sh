#!/usr/bin/env bash
# Run Python fitting with the correct conda environments.
# MuSiCal -> conda env "musical"
# SigProfilerAssignment -> conda env "pyenv"
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
export MPLCONFIGDIR="${ROOT}/results/logs/mpl"
mkdir -p "$MPLCONFIGDIR"

TOOLS_ARG=""
COSMIC_ARG=""
DATASETS_ARG=""
PHASE_ARG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tools) TOOLS_ARG="$2"; shift 2 ;;
    --cosmic) COSMIC_ARG="$2"; shift 2 ;;
    --datasets) DATASETS_ARG="$2"; shift 2 ;;
    --phase) PHASE_ARG="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

run_one() {
  local env="$1"
  local tool="$2"
  local py="/opt/anaconda3/envs/${env}/bin/python"
  if [[ ! -x "$py" ]]; then
    echo "ERROR: missing python at $py" >&2
    return 1
  fi
  local args=(scripts/python/run_py_fitting.py --tools "$tool")
  [[ -n "$COSMIC_ARG" ]] && args+=(--cosmic "$COSMIC_ARG")
  [[ -n "$DATASETS_ARG" ]] && args+=(--datasets "$DATASETS_ARG")
  [[ -n "$PHASE_ARG" ]] && args+=(--phase "$PHASE_ARG")
  echo "==> Running $tool with conda env $env"
  "$py" "${args[@]}"
}

if [[ -z "$TOOLS_ARG" || "$TOOLS_ARG" == *MuSiCal* ]]; then
  run_one musical MuSiCal
fi
if [[ -z "$TOOLS_ARG" || "$TOOLS_ARG" == *SigProfilerAssignment* ]]; then
  run_one pyenv SigProfilerAssignment
fi
