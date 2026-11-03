#!/usr/bin/env bash
# Package the experimental records of a run for archiving (e.g. on Zenodo).
#
#   scripts/make_release.sh results        -> dist/sdc-bench-records-<date>.tar.gz
#
# The archive holds the raw records, their checksums, the run metadata, the verifier and
# timing outputs, the logs, the data licence and the record-format description. The code
# itself is archived from the Git repository (see docs/reproducing.md).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RES="$(cd "${1:?usage: $0 <results-dir>}" && pwd)"
mkdir -p "$ROOT/dist"
( cd "$RES" && sha256sum -c SHA256SUMS >/dev/null ) && echo "record checksums: OK"
DATE=$(date -r "$RES/injections.bin" +%Y%m%d)
STAGE="$(mktemp -d)"; trap 'rm -rf "$STAGE"' EXIT
NAME="sdc-bench-records-$DATE"; mkdir "$STAGE/$NAME"
for f in injections.bin clean_checks.bin problems.csv campaign_meta.json SHA256SUMS \
         timing.csv verify_summary.csv verify_examples.csv system.txt logs; do
  [[ -e "$RES/$f" ]] && cp -r "$RES/$f" "$STAGE/$NAME/"
done
cp "$ROOT/LICENSE-DATA" "$STAGE/$NAME/LICENSE"
cp "$ROOT/docs/data-format.md" "$STAGE/$NAME/DATA_FORMAT.md"
tar -C "$STAGE" -czf "$ROOT/dist/$NAME.tar.gz" "$NAME"
ls -l "$ROOT/dist/$NAME.tar.gz"
