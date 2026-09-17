#!/usr/bin/env bash
# run_benchmark_v10_2.sh
#
# Supersedes run_benchmark_v10.sh. Produces Table 1 and the full fifteen-density
# Marron-Wand supplementary table in one pass, at 50 seeds, with paired differences
# and paired standard errors.
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash run_benchmark_v10_2.sh
#
# not as ./run_benchmark_v10_2.sh and not by sourcing it.
#
# Files this script expects beside itself: benchmark_v10_2.R and mw15.R.
#
# Set PKG_DIR to your checkout of github.com/mitch-thornton/kde-ad-heaping. The default
# follows the convention of keeping clones under /Users/mitch/src/<PROJECT>/. Tell me the
# real layout and I will change the default.
#
#   bash run_benchmark_v10_2.sh --check      fast reduced run, about 25 seconds
#   bash run_benchmark_v10_2.sh --time-sem   time one SEM cell before committing to a full run
#   bash run_benchmark_v10_2.sh              full run
#
# Environment overrides: PKG_DIR NSEED SET OUTDIR NO_SEM DELTA KAPPA

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
NSEED="${NSEED:-50}"
SET="${SET:-both}"
OUTDIR="${OUTDIR:-$PKG_DIR/results}"
RS="${RS:-$HERE/benchmark_v10_2.R}"
export MW15_R="${MW15_R:-$HERE/mw15.R}"

die() { echo "error: $*" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || die "Rscript is not on PATH"
[ -d "$PKG_DIR/R" ] || die "no R/ directory under PKG_DIR=$PKG_DIR ; set PKG_DIR to your kde-ad-heaping checkout"
[ -f "$RS" ]        || die "benchmark_v10_2.R not found at $RS"
[ -f "$MW15_R" ]    || die "mw15.R not found at $MW15_R"

case "${1:-}" in
  --check)
    echo "[check] reduced run, 2 seeds, 2 grids, 3 densities per set"
    SMOKE=1 NSEED=2 SET=both PKG_DIR="$PKG_DIR" OUTDIR="${OUTDIR}-smoke" Rscript "$RS"
    echo "[check] verifying the Marron-Wand definitions"
    Rscript "$HERE/verify_mw15.R" | sed -n '1,20p'
    echo "[check] passed"
    exit 0 ;;
  --time-sem)
    echo "[time-sem] timing one cell with and without the Kernelheaping SEM baseline."
    echo "           The full run is 76 cells at NSEED replicates each."
    Rscript -e 'if (!requireNamespace("Kernelheaping", quietly=TRUE))
      stop("Kernelheaping is not installed. install.packages(\"Kernelheaping\")")'
    echo "--- without SEM ---"
    time ( SMOKE=1 NSEED=2 SET=table1 NO_SEM=1 PKG_DIR="$PKG_DIR" OUTDIR=/tmp/t1 Rscript "$RS" >/dev/null )
    echo "--- with SEM ---"
    time ( SMOKE=1 NSEED=2 SET=table1 PKG_DIR="$PKG_DIR" OUTDIR=/tmp/t2 Rscript "$RS" >/dev/null )
    echo
    echo "Multiply the difference by 76 cells / 6 smoke cells and by NSEED/2 to project"
    echo "the full run. If that is too long, run SET=table1 with SEM and SET=mw15 with"
    echo "NO_SEM=1, and say in the supplementary caption that SEM is omitted there."
    exit 0 ;;
esac

echo "=== benchmark_v10_2 ==="
echo "  package dir : $PKG_DIR"
echo "  mw15.R      : $MW15_R"
echo "  seeds       : $NSEED"
echo "  set         : $SET"
echo "  output      : $OUTDIR"
echo "  SEM         : ${NO_SEM:+disabled}${NO_SEM:-enabled if Kernelheaping is installed}"
echo

Rscript -e 'if (!requireNamespace("Kernelheaping", quietly=TRUE))
  cat("note: Kernelheaping is not installed, so the SEM column will be empty.\n      install.packages(\"Kernelheaping\")\n\n")'

mkdir -p "$OUTDIR"
PKG_DIR="$PKG_DIR" NSEED="$NSEED" SET="$SET" OUTDIR="$OUTDIR" \
  Rscript "$RS" 2>&1 | tee "$OUTDIR/benchmark_v10_2.log"

echo
echo "wrote:"
echo "  $OUTDIR/benchmark_v10_2.rds"
echo "  $OUTDIR/benchmark_v10_2.json"
echo "  $OUTDIR/benchmark_v10_2.log"
echo
echo "Paste the log back and I will rebuild Table 1, the supplementary table, and the"
echo "narrative from the JSON."
