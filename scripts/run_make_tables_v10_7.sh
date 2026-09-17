#!/usr/bin/env bash
# run_make_tables_v10_7.sh
#
# Emits Supplementary Tables S4 to S9 and the numbers_v10_7.tex macro file from the
# stored result files, so no v10.3 number in the paper is typed by hand. Run
# run_make_tables_v10_2.sh first for Table 1 and Tables S2 and S3.
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash run_make_tables_v10_7.sh
#
# Needs, under RESULTS (default $PKG_DIR/results):
#   batch_v10_7.rds  test_lattice_v2.rds  verify_readers_v10_7.rds  benchmark_v10_2.rds

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
RESULTS="${RESULTS:-$PKG_DIR/results}"
OUTDIR="${OUTDIR:-$RESULTS/tables}"
RS="${RS:-$HERE/make_tables_v10_7.R}"

die() { echo "error: $*" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || die "Rscript is not on PATH"
[ -f "$RS" ] || die "make_tables_v10_7.R not found at $RS"
for f in batch_v10_7.rds test_lattice_v2.rds verify_readers_v10_7.rds benchmark_v10_2.rds; do
  [ -f "$RESULTS/$f" ] || die "$f not found in $RESULTS"
done

mkdir -p "$OUTDIR"
RESULTS="$RESULTS" OUTDIR="$OUTDIR" Rscript "$RS" | tee "$OUTDIR/tables_v10_7_build.log"

echo
echo "LaTeX written to $OUTDIR :"
echo "  supp_tableS4_thresholds.tex  supp_tableS5_cutoff.tex  supp_tableS6_deconv.tex"
echo "  supp_tableS7_nhanes_bw.tex   supp_tableS8_lattice.tex supp_tableS9_detector.tex"
echo "  numbers_v10_7.tex            macros the main text reads; numbers_v10_7.txt for the letter"
