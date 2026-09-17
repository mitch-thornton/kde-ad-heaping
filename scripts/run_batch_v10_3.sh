#!/usr/bin/env bash
# run_batch_v10_3.sh -- Step 3 of the v10.3 plan.
# You use zsh, so invoke this as:   bash run_batch_v10_3.sh
#
#   bash run_batch_v10_3.sh --check        fast smoke path, a few minutes
#   bash run_batch_v10_3.sh                full run
#   BLOCKS=B7,B9 bash run_batch_v10_3.sh   a subset
#
# B5 and B11 need the NHANES files. Point NHANES_DIR at the directory holding
# DEMO_J.xpt, BMX_J.xpt and SMQ_J.xpt. Without it those two blocks skip with a message
# and everything else still runs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
OUTDIR="${OUTDIR:-$PKG_DIR/results}"
RS="${RS:-$HERE/batch_v10_3.R}"
NHANES_DIR="${NHANES_DIR:-}"

command -v Rscript >/dev/null 2>&1 || { echo "error: Rscript not on PATH" >&2; exit 1; }
[ -d "$PKG_DIR/R" ] || { echo "error: no R/ under PKG_DIR=$PKG_DIR" >&2; exit 1; }
[ -f "$RS" ]        || { echo "error: batch_v10_3.R not found at $RS" >&2; exit 1; }

if [ "${1:-}" = "--check" ]; then
  echo "[check] smoke path: 3 replicates, 20 bootstrap resamples, 10 simulations"
  BLOCKS="${BLOCKS:-B6,B7,B9}" NSEED=3 NBOOT=20 NSIM=10 \
    PKG_DIR="$PKG_DIR" OUTDIR="${OUTDIR}-smoke" NHANES_DIR="$NHANES_DIR" Rscript "$RS"
  echo "[check] passed"; exit 0
fi

if [ -z "$NHANES_DIR" ]; then
  echo "note: NHANES_DIR is not set, so B5 and B11 will skip."
  echo "      set it to the directory holding DEMO_J.xpt, BMX_J.xpt and SMQ_J.xpt."
  echo
fi

mkdir -p "$OUTDIR"
PKG_DIR="$PKG_DIR" OUTDIR="$OUTDIR" NHANES_DIR="$NHANES_DIR" \
  NSEED="${NSEED:-20}" NBOOT="${NBOOT:-1000}" NSIM="${NSIM:-200}" BLOCKS="${BLOCKS:-B5,B6,B7,B8,B9,B11}" \
  Rscript "$RS" 2>&1 | tee "$OUTDIR/batch_v10_3.log"

echo
echo "wrote $OUTDIR/batch_v10_3.rds and $OUTDIR/batch_v10_3.log"
echo "Paste the log back and I will fold the numbers into v10.3."
