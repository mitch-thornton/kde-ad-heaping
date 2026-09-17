#!/usr/bin/env bash
# run_make_tables_v10_2.sh
#
# Rebuilds Table 1, Supplementary Table S2 and Supplementary Table S3 directly from
# benchmark_v10_2.rds, so no experimental value in the paper is typed by hand.
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash run_make_tables_v10_2.sh
#
# Defaults match the layout your v10.2 run reported:
#   package   /Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping
#   mw15.R    /Users/mitch/src/KDE-AD-HEAPING/mw15.R
#   results   /Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping/results
#
# Override with RESULTS and OUTDIR if you move things.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
RESULTS="${RESULTS:-$PKG_DIR/results}"
OUTDIR="${OUTDIR:-$RESULTS/tables}"
RS="${RS:-$HERE/make_tables_v10_2.R}"

die() { echo "error: $*" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || die "Rscript is not on PATH"
[ -f "$RS" ] || die "make_tables_v10_2.R not found at $RS"
[ -f "$RESULTS/benchmark_v10_2.rds" ] || \
  die "benchmark_v10_2.rds not found in $RESULTS ; run run_benchmark_v10_2.sh first"

mkdir -p "$OUTDIR"
RESULTS="$RESULTS" OUTDIR="$OUTDIR" Rscript "$RS" | tee "$OUTDIR/tables_v10_build.log"

echo
echo "LaTeX written to $OUTDIR :"
echo "  table1_v10.tex          replaces the Table 1 block in the paper"
echo "  supp_tableS2_mw15.tex   new supplementary table, all fifteen Marron-Wand densities"
echo "  supp_tableS3_paired.tex new supplementary table, paired differences and paired SEs"
echo "  tables_v10_summary.txt  counts and provenance, for the response letter"
echo
echo "Add \\input{} lines for these rather than pasting, so a rerun updates the PDF."
