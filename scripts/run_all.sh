#!/usr/bin/env bash
# Run the complete sdc-bench pipeline: build, campaign, verification, timing, analysis.
#
#   scripts/run_all.sh quick     reduced campaign (about a minute); used by CI and tests/check_quick.sh
#   scripts/run_all.sh full      the full campaign reported in the paper
#
# Environment overrides: OUT (results folder), THREADS (0 = all cores), K, CLEAN, STRIDE,
# PRECISIONS (comma-separated subset of fp64,fp32,fp16,bf16), SKIP_TIMING=1.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

MODE="${1:-full}"
case "$MODE" in
  quick) K="${K:-100}";  CLEAN="${CLEAN:-3}";  STRIDE="${STRIDE:-50}"; OUT="${OUT:-results-quick}" ;;
  full)  K="${K:-2000}"; CLEAN="${CLEAN:-50}"; STRIDE="${STRIDE:-10}"; OUT="${OUT:-results}" ;;
  *) echo "usage: $0 [quick|full]" >&2; exit 2 ;;
esac
PRECISIONS="${PRECISIONS:-fp64,fp32,fp16,bf16}"
THREADS="${THREADS:-0}"
mkdir -p "$OUT/logs"

step() { printf '\n=== %s ===\n' "$*"; }

step "0/6 environment"
command -v "${CXX:-g++}" >/dev/null || { echo "error: C++ compiler '${CXX:-g++}' not found" >&2; exit 1; }
command -v python3 >/dev/null     || { echo "error: python3 not found" >&2; exit 1; }
python3 -c "import numpy" 2>/dev/null || { echo "error: numpy not found (pip install -r analysis/requirements.txt)" >&2; exit 1; }
[[ "$(uname -m)" == "x86_64" ]] || echo "warning: not x86-64; results are valid but not bit-identical to the reference run"
{
  echo "date: $(date -Is)"
  echo "mode: $MODE  K=$K CLEAN=$CLEAN STRIDE=$STRIDE PRECISIONS=$PRECISIONS"
  uname -a
  (command -v lscpu >/dev/null && lscpu | grep -E "Model name|^CPU\(s\)|Thread|MHz|L1d|L2|L3") || grep -m1 "model name" /proc/cpuinfo || true
  "${CXX:-g++}" --version | head -1
  python3 --version
  python3 -c "import numpy; print('numpy', numpy.__version__)"
  echo "nproc: $(nproc)"
  (cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null | sed 's/^/governor: /') || true
} > "$OUT/system.txt"
cat "$OUT/system.txt"

step "1/6 build"
make 2>&1 | tee "$OUT/logs/build.log"

step "2/6 campaign (Algorithm 1)"
{ time ./bin/campaign --out "$OUT" --K "$K" --clean "$CLEAN" --precisions "$PRECISIONS" --threads "$THREADS"; } 2>&1 | tee "$OUT/logs/campaign.log"

step "3/6 independent verification (Algorithm 3)"
set +e
./bin/verify --dir "$OUT" --stride "$STRIDE" --threads "$THREADS" 2>&1 | tee "$OUT/logs/verify.log"
VSTATUS=${PIPESTATUS[0]}
set -e
if [[ $VSTATUS -eq 3 ]]; then
  echo "error: the verifier found bit-level disagreements with the campaign (see $OUT/verify_summary.csv)" >&2
  exit 3
fi

step "4/6 timing benchmark"
if [[ "${SKIP_TIMING:-0}" == "1" || "$MODE" == "quick" ]]; then
  echo "skipped"
else
  echo "Run on an idle machine on mains power; this step takes a few minutes."
  ./bin/bench_timing | tee "$OUT/timing.csv"
fi

step "5/6 analysis (Algorithm 2)"
python3 analysis/analyze.py "$OUT" 2>&1 | tee "$OUT/logs/analyze.log"

step "6/6 checksums"
( cd "$OUT" && sha256sum injections.bin clean_checks.bin problems.csv > SHA256SUMS && cat SHA256SUMS )

printf '\nDone. Report: %s/report.md   All numbers: %s/summary.json   Figure data: %s/fig/\n' "$OUT" "$OUT" "$OUT"
