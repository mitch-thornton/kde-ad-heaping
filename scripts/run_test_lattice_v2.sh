#!/usr/bin/env bash
# run_test_lattice_v2.sh -- Decision 1: validate and adopt the repaired mixed-grain reader.
# You use zsh, so invoke as:   bash run_test_lattice_v2.sh
# Needs heap_lattice_v2.R beside it. About a minute with the NHANES files attached.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
NHANES_DIR="${NHANES_DIR:-/Users/mitch/src/KDE-AD-HEAPING/data/NHANES}"
OUTDIR="${OUTDIR:-$PKG_DIR/results}"
RS="${RS:-$HERE/test_lattice_v2.R}"
export LATTICE_V2="${LATTICE_V2:-$HERE/heap_lattice_v2.R}"
command -v Rscript >/dev/null 2>&1 || { echo "error: Rscript not on PATH" >&2; exit 1; }
[ -d "$PKG_DIR/R" ]      || { echo "error: no R/ under PKG_DIR=$PKG_DIR" >&2; exit 1; }
[ -f "$RS" ]             || { echo "error: test_lattice_v2.R not found at $RS" >&2; exit 1; }
[ -f "$LATTICE_V2" ]     || { echo "error: heap_lattice_v2.R not found at $LATTICE_V2" >&2; exit 1; }
[ -f "$NHANES_DIR/SMQ_J.xpt" ] && echo "NHANES found, parts 3 and 4 will run." \
  || echo "note: SMQ_J.xpt not found in $NHANES_DIR, parts 3 and 4 will skip."
mkdir -p "$OUTDIR"
PKG_DIR="$PKG_DIR" NHANES_DIR="$NHANES_DIR" OUTDIR="$OUTDIR" \
  NSIM="${NSIM:-200}" NBOOT="${NBOOT:-1000}" Rscript "$RS" 2>&1 | tee "$OUTDIR/test_lattice_v2.log"
echo
echo "Paste the log back. Part 3 is the one that decides whether the published numbers stand."
