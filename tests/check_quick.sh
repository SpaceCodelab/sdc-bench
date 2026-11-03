#!/usr/bin/env bash
# Determinism regression test.
#
# Runs the quick campaign into a temporary folder and checks that the records are
# bit-identical to the reference digests in tests/quick_reference.sha256, then checks that
# the independent verifier reports no implementation disagreement and that the analysis
# runs. A mismatch on x86-64 means the arithmetic or the generators have changed.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="$(mktemp -d)"; trap 'rm -rf "$OUT"' EXIT

OUT="$OUT" SKIP_TIMING=1 scripts/run_all.sh quick > "$OUT/run.log" 2>&1 || {
  echo "FAIL: pipeline exited with an error"; tail -30 "$OUT/run.log"; exit 1; }

if [[ "$(uname -m)" != "x86_64" ]]; then
  echo "SKIP: digest comparison is defined for x86-64 only"; exit 0
fi
( cd "$OUT" && sha256sum -c "$ROOT/tests/quick_reference.sha256" ) || {
  echo "FAIL: records differ from the reference digests"; exit 1; }
grep -q " 0 output mismatches, 0 check mismatches" "$OUT/logs/verify.log" || {
  echo "FAIL: verifier reported disagreements"; exit 1; }
[[ -s "$OUT/summary.json" ]] || { echo "FAIL: analysis produced no summary"; exit 1; }
echo "PASS: quick campaign is bit-identical to the reference, verifier clean, analysis ran"
