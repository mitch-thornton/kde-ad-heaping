#!/usr/bin/env bash
# run_make_tableS12_v13.sh
#
# Emits Supplementary Table S12 and the numbers_v13 macros from
# results/sensitivity_n_v13.rds. build.sh calls the same script, so this wrapper exists only
# for running the emission on its own.
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash scripts/run_make_tableS12_v13.sh
#
# TARGET selects the density by its Table 1 index. The default is T2, the bimodal mixture,
# which is the one representative scenario the editor asked for. TARGET=T1,T2 emits both
# targets; read the v13.0 notes before doing that, since the Gaussian row carries two
# concessions nobody asked for.
#
# Environment overrides: RESULTS OUTDIR TARGET

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE="$(cd "$HERE/.." && pwd)"
RESULTS="${RESULTS:-$BUNDLE/results}"
OUTDIR="${OUTDIR:-$BUNDLE/tables}"

die() { echo "error: $*" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || die "Rscript is not on PATH"
[ -f "$RESULTS/sensitivity_n_v13.rds" ] || die "results/sensitivity_n_v13.rds is missing.
       Run  bash scripts/run_sensitivity_n_v13.sh  first."
[ -f "$RESULTS/benchmark_v10_2.rds" ] || die "results/benchmark_v10_2.rds is missing, and the
       emitter re-checks the n = 4000 row against it before writing anything."

echo "=== make_tableS12_v13 ==="
echo "  results : $RESULTS"
echo "  output  : $OUTDIR"
echo "  target  : ${TARGET:-T2}   (T2 is the bimodal mixture)"
echo

export RESULTS OUTDIR
if [ -n "${TARGET:-}" ]; then export TARGET; fi
Rscript "$HERE/make_tableS12_v13.R"

echo
echo "wrote:"
echo "  $OUTDIR/supp_tableS12_n.tex"
echo "  $OUTDIR/numbers_v13.tex"
echo "  $OUTDIR/numbers_v13.txt"
