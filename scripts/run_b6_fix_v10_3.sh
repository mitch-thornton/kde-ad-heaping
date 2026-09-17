#!/usr/bin/env bash
# run_b6_fix_v10_3.sh -- corrected B6, replacing the double-rounded version.
# You use zsh, so invoke as:   bash run_b6_fix_v10_3.sh
# Needs SMQ_J.xpt. Takes about a minute.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
NHANES_DIR="${NHANES_DIR:-/Users/mitch/src/KDE-AD-HEAPING/data/NHANES}"
OUTDIR="${OUTDIR:-$PKG_DIR/results}"
RS="${RS:-$HERE/b6_fix_v10_3.R}"
command -v Rscript >/dev/null 2>&1 || { echo "error: Rscript not on PATH" >&2; exit 1; }
[ -d "$PKG_DIR/R" ] || { echo "error: no R/ under PKG_DIR=$PKG_DIR" >&2; exit 1; }
[ -f "$NHANES_DIR/SMQ_J.xpt" ] || { echo "error: SMQ_J.xpt not found in $NHANES_DIR" >&2; exit 1; }
mkdir -p "$OUTDIR"
PKG_DIR="$PKG_DIR" NHANES_DIR="$NHANES_DIR" OUTDIR="$OUTDIR" NSIM="${NSIM:-200}" \
  Rscript "$RS" 2>&1 | tee "$OUTDIR/b6_fix_v10_3.log"
echo
echo "Paste the log back. The gate line under each base is the one that matters."
