#!/usr/bin/env bash
# run_verify_readers_v10_3.sh -- Step 2 of the v10.3 plan.
# You use zsh, so invoke this as:   bash run_verify_readers_v10_3.sh
# Takes about a minute. Nothing is modified; the shipped package is only called.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
NSEED="${NSEED:-10}"
OUTDIR="${OUTDIR:-$PKG_DIR/results}"
RS="${RS:-$HERE/verify_readers_v10_3.R}"
command -v Rscript >/dev/null 2>&1 || { echo "error: Rscript not on PATH" >&2; exit 1; }
[ -d "$PKG_DIR/R" ] || { echo "error: no R/ under PKG_DIR=$PKG_DIR" >&2; exit 1; }
[ -f "$RS" ]        || { echo "error: verify_readers_v10_3.R not found at $RS" >&2; exit 1; }
mkdir -p "$OUTDIR"
PKG_DIR="$PKG_DIR" NSEED="$NSEED" OUTDIR="$OUTDIR" Rscript "$RS" 2>&1 \
  | tee "$OUTDIR/verify_readers_v10_3.log"
echo
echo "Paste the whole output back. The two verdict lines at the bottom are what I need."
