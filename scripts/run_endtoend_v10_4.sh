#!/usr/bin/env bash
# run_endtoend_v10_4.sh
#
# The two substantive reviewer items still open after v10.4, in one pass:
#   A  end-to-end grid detection, Reviewer 1 item 3 and the editor's automatic-detection item
#   B  known-truth test of the heaped-fraction reader, the other half of Reviewer 1 item 5
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash run_endtoend_v10_4.sh
#
# Defaults match your layout:
#   package   /Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping
#   mw15.R    beside this script, else /Users/mitch/src/KDE-AD-HEAPING/mw15.R
#   results   $PKG_DIR/results
#
# Fast check first (about half a minute), then the full run (a few minutes):
#
#     SMOKE=1 bash run_endtoend_v10_4.sh
#     bash run_endtoend_v10_4.sh
#
# Override NSEED, NSEEDB, PART, OUTDIR from the environment if you want a narrower run.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
OUTDIR="${OUTDIR:-$PKG_DIR/results}"
RS="${RS:-$HERE/endtoend_v10_4.R}"

if [ -z "${MW15_R:-}" ]; then
  if   [ -f "$HERE/mw15.R" ];                       then MW15_R="$HERE/mw15.R"
  elif [ -f "$PKG_DIR/../mw15.R" ];                 then MW15_R="$(cd "$PKG_DIR/.." && pwd)/mw15.R"
  elif [ -f "$PKG_DIR/scripts/mw15.R" ];            then MW15_R="$PKG_DIR/scripts/mw15.R"
  fi
fi

die() { echo "error: $*" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || die "Rscript is not on PATH"
[ -f "$RS" ] || die "endtoend_v10_4.R not found at $RS"
[ -d "$PKG_DIR/R" ] || die "no R/ directory under PKG_DIR=$PKG_DIR"
[ -n "${MW15_R:-}" ] && [ -f "$MW15_R" ] || die "mw15.R not found; set MW15_R"

mkdir -p "$OUTDIR"
echo "package : $PKG_DIR"
echo "mw15.R  : $MW15_R"
echo "results : $OUTDIR"
echo

PKG_DIR="$PKG_DIR" MW15_R="$MW15_R" OUTDIR="$OUTDIR" \
  Rscript "$RS" 2>&1 | tee "$OUTDIR/endtoend_v10_4.log"

echo
echo "log: $OUTDIR/endtoend_v10_4.log"
echo "rds: $OUTDIR/endtoend_v10_4.rds"
echo "Send both back and I will fold them into v10.5."
