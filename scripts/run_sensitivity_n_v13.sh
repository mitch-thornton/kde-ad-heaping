#!/usr/bin/env bash
# run_sensitivity_n_v13.sh
#
# The sample-size sensitivity run asked for by Reviewer 1 and by the editor in round five.
# One grid width, five sample sizes, three estimators, fifty seeds, nested draws.
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash scripts/run_sensitivity_n_v13.sh
#
# not as ./run_sensitivity_n_v13.sh and not by sourcing it.
#
#   bash scripts/run_sensitivity_n_v13.sh --check    three seeds, two sizes, a few seconds
#   bash scripts/run_sensitivity_n_v13.sh            the full run
#
# Two things this script will not let you get wrong.
#
# The output goes into THIS BUNDLE'S results/ directory by default, not into the package
# checkout's results/, because the table emitter reads the bundle's results/. Override with
# OUTDIR if you want it elsewhere.
#
# The run is checked against the published Table 1 cell before it is believed. The n = 4000
# row has to reproduce, per seed and in the fifty-seed mean, the bimodal and Gaussian
# entries at D = 0.5 already stored in results/benchmark_v10_2.rds. The script exits
# nonzero if it does not, after writing the result files so the failure can be diagnosed.
#
# Environment overrides: PKG_DIR OUTDIR BENCH_RDS TARGETS NS GRID_D NSEED CHECK_TOL

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE="$(cd "$HERE/.." && pwd)"

PKG_DIR="${PKG_DIR:-/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping}"
OUTDIR="${OUTDIR:-$BUNDLE/results}"
RS="${RS:-$HERE/sensitivity_n_v13.R}"
export MW15_R="${MW15_R:-$HERE/mw15.R}"
export BENCH_R="${BENCH_R:-$HERE/benchmark_v10_2.R}"

die() { echo "error: $*" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || die "Rscript is not on PATH"
[ -f "$RS" ]        || die "sensitivity_n_v13.R not found at $RS"
[ -f "$MW15_R" ]    || die "mw15.R not found at $MW15_R"
[ -f "$BENCH_R" ]   || die "benchmark_v10_2.R not found at $BENCH_R (the target densities and the grid are read out of it)"

# The package checkout is the source of the estimators. If the usual path is not there,
# fall back to this bundle's own repository image, which is the same code at adheaping
# 1.1.0, and say so rather than failing on a path guess.
if [ ! -d "$PKG_DIR/R" ]; then
  if [ -d "$BUNDLE/github/R" ]; then
    echo "note: no R/ under PKG_DIR=$PKG_DIR"
    echo "      falling back to this bundle's repository image, $BUNDLE/github"
    echo "      Set PKG_DIR explicitly if you want your working checkout instead."
    echo
    PKG_DIR="$BUNDLE/github"
  else
    die "no R/ directory under PKG_DIR=$PKG_DIR ; set PKG_DIR to your kde-ad-heaping checkout"
  fi
fi

# The stored benchmark lives with the bundle, not with wherever this run is told to write,
# so it is anchored on $BUNDLE and not on $OUTDIR. Overriding OUTDIR must not quietly
# switch the consistency check off.
BENCH_RDS="${BENCH_RDS:-$BUNDLE/results/benchmark_v10_2.rds}"
if [ ! -f "$BENCH_RDS" ]; then
  if [ -n "${ALLOW_NO_ANCHOR:-}" ]; then
    echo "warning: $BENCH_RDS is missing and ALLOW_NO_ANCHOR is set, so this run will"
    echo "         not be checked against the published Table 1 cell."
    echo
  else
    die "the stored benchmark is missing at $BENCH_RDS, so the n = 4000 row cannot be
       checked against the published Table 1 cell and the run would be unverifiable.
       Point BENCH_RDS at results/benchmark_v10_2.rds, or set ALLOW_NO_ANCHOR=1 if you
       really want an unchecked run."
  fi
fi

if [ "${1:-}" = "--check" ]; then
  echo "[check] three seeds at n = 250 and n = 4000, both targets"
  echo "[check] the n = 4000 values must match the stored per-seed values exactly"
  mkdir -p "${OUTDIR}-smoke"
  SMOKE=1 PKG_DIR="$PKG_DIR" OUTDIR="${OUTDIR}-smoke" BENCH_RDS="$BENCH_RDS" Rscript "$RS"
  echo "[check] passed"
  exit 0
fi

echo "=== sensitivity_n_v13 ==="
echo "  package dir : $PKG_DIR"
echo "  benchmark R : $BENCH_R"
echo "  stored cell : $BENCH_RDS"
echo "  output      : $OUTDIR"
echo "  targets     : ${TARGETS:-T1,T2}   (Gaussian, Bimodal)"
echo "  sizes       : ${NS:-250,500,1000,2000,4000}"
echo "  grid width  : ${GRID_D:-0.5}"
echo "  seeds       : ${NSEED:-50}"
echo

mkdir -p "$OUTDIR"
# Pass only what was actually set. A word that expands to NAME=value is not an assignment
# to bash, so these have to be exports rather than a command prefix.
# Written as if blocks rather than with &&. Under set -e a bare `test && export` whose test
# fails ends the script, which is the failure mode that silently skipped an edit in v12.3.
export PKG_DIR OUTDIR BENCH_RDS
if [ -n "${TARGETS:-}" ];   then export TARGETS;   fi
if [ -n "${NS:-}" ];        then export NS;        fi
if [ -n "${GRID_D:-}" ];    then export GRID_D;    fi
if [ -n "${NSEED:-}" ];     then export NSEED;     fi
if [ -n "${CHECK_TOL:-}" ]; then export CHECK_TOL; fi
set +e
Rscript "$RS" 2>&1 | tee "$OUTDIR/sensitivity_n_v13.log"
rc=${PIPESTATUS[0]}
set -e

echo
echo "wrote:"
echo "  $OUTDIR/sensitivity_n_v13.rds"
echo "  $OUTDIR/sensitivity_n_v13.json"
echo "  $OUTDIR/sensitivity_n_v13.log"
echo

if [ "$rc" -ne 0 ]; then
  echo "The Table 1 anchor failed. Send back the log. Do not read the numbers." >&2
  exit "$rc"
fi

echo "Send back sensitivity_n_v13.log and sensitivity_n_v13.rds and the supplementary"
echo "table and its macros will be emitted from the rds, never typed."
