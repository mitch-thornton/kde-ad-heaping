#!/usr/bin/env bash
# run_make_figures_v12.sh
#
# Regenerates Figures 1, 2 and 4 in R and stores the NHANES study at full precision,
# closing pending markers B15, B2, B3 and B4.
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash run_make_figures_v12.sh
#
# Fast check first (about a minute), then the full run:
#
#     SMOKE=1 bash run_make_figures_v12.sh
#     bash run_make_figures_v12.sh
#
# Defaults match your layout. Figures go to ../figures beside this script, so running it
# from a bundle's scripts/ directory writes into that bundle's figures/ directory.
# Result files go to the package results/ store alongside everything else.
#
# Parts: 1 is Figure 1 and takes a second, 2 is Figure 2 and takes a few minutes at 50
# seeds, 3 is Figure 4 and Table 2 and needs the NHANES files. Run one at a time with
# PARTS=2 if you prefer.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
NHANES_DIR="${NHANES_DIR:-/Users/mitch/src/KDE-AD-HEAPING/data/NHANES}"
OUTDIR="${OUTDIR:-$PKG_DIR/results}"
FIGDIR="${FIGDIR:-$(cd "$HERE/.." && pwd)/figures}"
RS="${RS:-$HERE/make_figures_v12.R}"

if [ -z "${MW15_R:-}" ]; then
  if   [ -f "$HERE/mw15.R" ];            then MW15_R="$HERE/mw15.R"
  elif [ -f "$PKG_DIR/../mw15.R" ];      then MW15_R="$(cd "$PKG_DIR/.." && pwd)/mw15.R"
  elif [ -f "$PKG_DIR/scripts/mw15.R" ]; then MW15_R="$PKG_DIR/scripts/mw15.R"
  fi
fi

die() { echo "error: $*" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || die "Rscript is not on PATH"
[ -f "$RS" ] || die "make_figures_v12.R not found at $RS"
[ -d "$PKG_DIR/R" ] || die "no R/ directory under PKG_DIR=$PKG_DIR"
[ -n "${MW15_R:-}" ] && [ -f "$MW15_R" ] || die "mw15.R not found; set MW15_R"

mkdir -p "$OUTDIR" "$FIGDIR"
echo "package : $PKG_DIR"
echo "mw15.R  : $MW15_R"
echo "NHANES  : $NHANES_DIR"
echo "figures : $FIGDIR"
echo "results : $OUTDIR"
echo

PKG_DIR="$PKG_DIR" MW15_R="$MW15_R" NHANES_DIR="$NHANES_DIR" \
OUTDIR="$OUTDIR" FIGDIR="$FIGDIR" \
  Rscript "$RS" 2>&1 | tee "$OUTDIR/make_figures_v10_7.log"

echo
echo "log: $OUTDIR/make_figures_v10_7.log"
echo "Three PDFs were overwritten in $FIGDIR. Rebuild the paper to see them:"
echo "    cd $(cd "$HERE/.." && pwd) && bash build.sh"
echo "Send back the log and fig1b_v10_7.rds, fig2_ise_v10_7.rds and nhanes_v10_7.rds"
echo "and I will wire the captions and Table 2 to them."
