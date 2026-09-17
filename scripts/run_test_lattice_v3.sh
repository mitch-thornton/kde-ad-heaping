#!/usr/bin/env bash
# run_test_lattice_v3.sh
# Regime map for the mixed-grain reader, the third reader on the real counts, and the
# model-based refinement. You use zsh; invoke as
#     NHANES_DIR=/Users/mitch/src/KDE-AD-HEAPING/data/NHANES bash run_test_lattice_v3.sh
# Expects heap_lattice_v3.R and test_lattice_v3.R beside this script (Downloads is fine).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
export PKG_DIR
export LATTICE_V3="${LATTICE_V3:-$HERE/heap_lattice_v3.R}"
export NHANES_DIR="${NHANES_DIR:-/Users/mitch/src/KDE-AD-HEAPING/data/NHANES}"
export OUTDIR="${OUTDIR:-$PKG_DIR/results}"
export NSIM="${NSIM:-200}" NBOOT="${NBOOT:-1000}"
[ -f "$LATTICE_V3" ] || { echo "heap_lattice_v3.R not found at $LATTICE_V3" >&2; exit 1; }
mkdir -p "$OUTDIR"
Rscript "$HERE/test_lattice_v3.R" 2>&1 | tee "$OUTDIR/test_lattice_v3.log"
echo; echo "log: $OUTDIR/test_lattice_v3.log   rds: $OUTDIR/test_lattice_v3.rds"
