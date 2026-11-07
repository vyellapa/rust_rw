#!/usr/bin/env bash
# Byte-for-byte comparison of delly-rs against a reference C++ delly build.
#   usage: compare.sh <reference-delly> <rust-delly> <outdir> <name> [delly sr args...]
# Compares the BCF, its .csi index and (if -d is given) the gzipped read dump.
set -euo pipefail
if [ "$#" -lt 4 ]; then
  echo "usage: $0 <reference-delly> <rust-delly> <outdir> <name> [delly sr args...]" >&2
  exit 2
fi
REF=$1; RS=$2; OUT=$3; NAME=$4; shift 4
case "$NAME" in ""|*/*) echo "name must be a nonempty filename component" >&2; exit 2;; esac
command -v bcftools >/dev/null || { echo "bcftools is required" >&2; exit 2; }
mkdir -p "$OUT"
args=(); dump=""
while [ $# -gt 0 ]; do
  case "$1" in
    -d|--dump)
      [ "$#" -ge 2 ] || { echo "$1 requires a filename" >&2; exit 2; }
      dump=$2; shift 2;;
    --dump=*) dump=${1#--dump=}; shift;;
    -d?*) dump=${1#-d}; shift;;
    -o*|--outfile*) echo "compare.sh chooses the output paths; omit -o/--outfile" >&2; exit 2;;
    --) args+=("$@"); break;;
    *) args+=("$1"); shift;;
  esac
done
# Never compare files left behind by a previous run. Keep final logs even on failure.
RUN=$(mktemp -d "$OUT/.${NAME}.XXXXXX")
trap 'rm -rf -- "$RUN"' EXIT
failed=0
for impl in ref rs; do
  bin=$REF; [ $impl = rs ] && bin=$RS
  extra=(); [ -n "$dump" ] && extra=(-d "$RUN/$NAME.$impl.dump.gz")
  if "$bin" sr -o "$RUN/$NAME.$impl.bcf" ${extra[@]+"${extra[@]}"} ${args[@]+"${args[@]}"} 2> "$RUN/$NAME.$impl.log"; then
    :
  else
    code=$?
    echo "$NAME: $impl binary failed (exit $code)" >&2
    failed=1
  fi
  cp "$RUN/$NAME.$impl.log" "$OUT/$NAME.$impl.log"
done
[ "$failed" -eq 0 ] || exit 1
suffixes=(bcf bcf.csi)
[ -z "$dump" ] || suffixes+=(dump.gz)
status=IDENTICAL
for suffix in "${suffixes[@]}"; do
  for impl in ref rs; do
    if [ ! -s "$RUN/$NAME.$impl.$suffix" ]; then
      echo "$NAME: missing or empty $impl $suffix output" >&2
      exit 1
    fi
  done
  if ! cmp -s "$RUN/$NAME.ref.$suffix" "$RUN/$NAME.rs.$suffix"; then
    status="DIFFER ($suffix)"
    failed=1
  fi
done
# pipefail also rejects unreadable/truncated BCF instead of reporting zero SVs.
n=$(bcftools view -H "$RUN/$NAME.ref.bcf" | wc -l | tr -d '[:space:]')
bcftools view -H "$RUN/$NAME.rs.bcf" >/dev/null
for impl in ref rs; do
  for suffix in "${suffixes[@]}"; do
    mv "$RUN/$NAME.$impl.$suffix" "$OUT/$NAME.$impl.$suffix"
  done
done
echo "$NAME: $status ($n SVs)"
exit "$failed"
