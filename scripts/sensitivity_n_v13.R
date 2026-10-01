#!/usr/bin/env Rscript
## sensitivity_n_v13.R
## Sample-size sensitivity, requested by Reviewer 1 and by the editor in round five.
##
## One grid width, several sample sizes, three estimators, fifty seeds. The question is
## what happens to the de-heaping advantage as the sample shrinks, which is the empirical
## counterpart of the Methods paragraph stating that the empirical characteristic function
## carries a variance falling as 1/n, that a smaller sample raises the noise floor, and
## that the identifiability limit pi/D does not move with n.
##
## Three properties of the construction matter and are not negotiable.
##
##   1. Nothing is retyped. The four target densities, the evaluation grid, the grid
##      spacing, the integrated-squared-error functional and the benchmark sample size are
##      read out of scripts/benchmark_v10_2.R by parsing that file and evaluating only the
##      assignments named below. The two scripts therefore cannot disagree about what the
##      bimodal mixture is or where the grid sits. If the extraction fails, or if what
##      comes back is not what Table 1 was built from, this script stops before computing
##      anything.
##
##   2. The draws are nested. Each seed draws the benchmark's full 4000 observations from
##      the same seed under the same rule, rounds them to the grid, and the smaller
##      samples are the leading n of that rounded draw. So the seed rule of record is
##      unchanged, nothing in the reproducibility section needs extending, and the
##      comparison across n is paired rather than independent.
##
##   3. The largest sample is a published cell, so it is an assertion rather than a
##      coincidence. At n = 4000 every per-seed value must reproduce the corresponding
##      per-seed value stored in results/benchmark_v10_2.rds, and the fifty-seed means
##      must reproduce the published Table 1 entries. The check runs after the results
##      have been written, so a failure costs the diagnosis but not the run. If it fails,
##      the run is wrong and nothing else in the output should be read.
##
## Usage:
##   PKG_DIR=/path/to/kde-ad-heaping Rscript sensitivity_n_v13.R
##   PKG_DIR=... SMOKE=1 Rscript sensitivity_n_v13.R
##
## Environment:
##   PKG_DIR   checkout of github.com/mitch-thornton/kde-ad-heaping            (required)
##   BENCH_R   path to benchmark_v10_2.R            (default: beside this script)
##   MW15_R    path to mw15.R                       (default: beside this script)
##   OUTDIR    results directory                    (default: $PKG_DIR/results)
##   BENCH_RDS stored benchmark result for the check (default: $OUTDIR/benchmark_v10_2.rds)
##   TARGETS   comma-separated Table 1 indices      (default T1,T2)
##   NS        comma-separated sample sizes         (default 250,500,1000,2000,4000)
##   GRID_D    the single grid width                (default 0.5)
##   NSEED     replicates per cell                  (default 50)
##   CHECK_TOL per-seed tolerance on ISE x 1e3      (default 1e-4)
##   SMOKE     set to 1 for a fast reduced run

t_start <- Sys.time()

.selfdir <- function() {
  a <- commandArgs(trailingOnly = FALSE); fa <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(fa)) dirname(normalizePath(fa)) else "."
}
HERE <- .selfdir()

pkg <- Sys.getenv("PKG_DIR", unset = NA)
if (is.na(pkg)) stop("set PKG_DIR to your kde-ad-heaping checkout")
rdir <- file.path(pkg, "R")
if (!dir.exists(rdir)) stop("no R/ directory under PKG_DIR=", pkg)
for (f in list.files(rdir, pattern = "\\.R$", full.names = TRUE)) source(f)

mwfile <- Sys.getenv("MW15_R", unset = file.path(HERE, "mw15.R"))
if (!file.exists(mwfile)) stop("mw15.R not found at ", mwfile, " (set MW15_R)")
source(mwfile)

benchR <- Sys.getenv("BENCH_R", unset = file.path(HERE, "benchmark_v10_2.R"))
if (!file.exists(benchR))
  stop("benchmark_v10_2.R not found at ", benchR, " (set BENCH_R). ",
       "This script reads the target densities and the grid out of it rather than ",
       "carrying its own copies.")

SMOKE     <- nzchar(Sys.getenv("SMOKE", ""))
NSEED     <- as.integer(Sys.getenv("NSEED", if (SMOKE) "3" else "50"))
GRID_D    <- as.numeric(Sys.getenv("GRID_D", "0.5"))
CHECK_TOL <- as.numeric(Sys.getenv("CHECK_TOL", "1e-4"))
TARGETS   <- strsplit(Sys.getenv("TARGETS", "T1,T2"), ",")[[1]]
NS        <- as.integer(strsplit(Sys.getenv("NS",
               if (SMOKE) "250,4000" else "250,500,1000,2000,4000"), ",")[[1]])
OUT       <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BENCH_RDS <- Sys.getenv("BENCH_RDS", unset = file.path(OUT, "benchmark_v10_2.rds"))

## ---- read the definitions out of the benchmark script ----------------------------------
## Only top-level assignments to the names below are evaluated, in file order, in a private
## environment. Everything else in that file, including its PKG_DIR check and its loops, is
## parsed and skipped.
.lift <- function(path, wanted) {
  exprs <- parse(path)
  env <- new.env(parent = globalenv())
  got <- character(0)
  for (e in exprs) {
    if (!is.call(e) || length(e) < 3L) next
    op <- tryCatch(as.character(e[[1L]]), error = function(x) "")
    if (!length(op) || !op[1L] %in% c("<-", "=", "<<-")) next
    if (!is.name(e[[2L]])) next
    nm <- as.character(e[[2L]])
    if (!nm %in% wanted) next
    eval(e, envir = env)
    got <- c(got, nm)
  }
  miss <- setdiff(wanted, got)
  if (length(miss))
    stop("could not lift ", paste(miss, collapse = ", "), " out of ", path,
         ". That file has been restructured, so this script's claim that it shares the ",
         "benchmark's definitions no longer holds. Stopping rather than guessing.")
  env
}

B <- .lift(benchR, c("TABLE1", "M", "grid", "dx", "n", "ise"))
TABLE1 <- B$TABLE1; M <- B$M; grid <- B$grid; dx <- B$dx; n_bench <- B$n; ise <- B$ise

## What was lifted has to be what Table 1 was built from. These are cheap and they are the
## only thing standing between a silent restructuring of the benchmark file and a
## supplementary table that disagrees with the main text.
if (!identical(as.integer(M), 2048L)) stop("lifted M is ", M, ", expected 2048")
if (length(grid) != 2048L) stop("lifted grid has length ", length(grid), ", expected 2048")
if (abs(grid[1] + 10) > 1e-12 || abs(dx - 20 / 2048) > 1e-12)
  stop("lifted grid does not run from -10 at spacing 20/2048")
if (!identical(as.integer(n_bench), 4000L)) stop("lifted benchmark n is ", n_bench, ", expected 4000")
if (!is.function(ise)) stop("lifted ise is not a function")
if (length(TABLE1) != 4L) stop("lifted TABLE1 has ", length(TABLE1), " entries, expected 4")
names(TABLE1) <- sapply(TABLE1, function(d) as.character(d$index))
if (!identical(names(TABLE1), c("T1", "T2", "T3", "T4")))
  stop("lifted TABLE1 indices are ", paste(names(TABLE1), collapse = ","), ", expected T1..T4")
if (!identical(TABLE1$T1$name, "Gaussian") || !identical(TABLE1$T2$name, "Bimodal"))
  stop("lifted TABLE1 names T1,T2 are ", TABLE1$T1$name, ",", TABLE1$T2$name,
       ", expected Gaussian,Bimodal")
if (!isTRUE(all.equal(TABLE1$T2$w, c(.5, .5))) ||
    !isTRUE(all.equal(TABLE1$T2$mu, c(-1.2, 1.2))) ||
    !isTRUE(all.equal(TABLE1$T2$sd, c(.5, .5))))
  stop("lifted bimodal mixture parameters are not the published ones")
if (!isTRUE(all.equal(ise(c(1, 0), c(0, 0)), 1 * dx)))
  stop("lifted ise does not integrate against the lifted dx")

bad <- setdiff(TARGETS, names(TABLE1))
if (length(bad)) stop("TARGETS names no such Table 1 index: ", paste(bad, collapse = ","))
if (any(NS < 50L)) stop("NS below 50 is not a meaningful sample for this estimator")
if (max(NS) > n_bench)
  stop("NS asks for ", max(NS), " but the nested draw is the benchmark's ", n_bench,
       ". Sizes above the benchmark's n would need an independent draw and would break ",
       "both the seed rule of record and the Table 1 anchor.")
NS <- sort(unique(NS))

METHODS <- c("naive", "deconv", "combined")
N_DRAW  <- n_bench          # always the benchmark's draw; the smaller samples are prefixes

cat(sprintf(paste0("sensitivity_n_v13  D=%.2f seeds=%d targets=%s n=%s  draw=%d  M=%d ",
                   "dx=%.6f%s\n"),
            GRID_D, NSEED, paste(TARGETS, collapse = "+"), paste(NS, collapse = ","),
            N_DRAW, M, dx, if (SMOKE) "  [SMOKE]" else ""))
cat(sprintf("definitions lifted from %s\n\n", benchR))
cat(sprintf("%-10s %6s %10s %10s %10s %8s %10s %7s %6s %6s\n",
            "density", "n", "naive", "deconv", "combined", "c/naive",
            "naive-c", "se", "wins", "D/h"))

rows <- list()
for (ti in TARGETS) {
  d  <- TABLE1[[ti]]
  ft <- mw_d(grid, d)

  ## per-seed ISE, indexed [seed, method, n]
  E <- array(NA_real_, dim = c(NSEED, length(METHODS), length(NS)),
             dimnames = list(NULL, METHODS, as.character(NS)))
  RHO  <- matrix(NA_real_, NSEED, length(NS), dimnames = list(NULL, as.character(NS)))
  DOH  <- matrix(NA_real_, NSEED, length(NS), dimnames = list(NULL, as.character(NS)))
  PICK <- matrix(NA_character_, NSEED, length(NS), dimnames = list(NULL, as.character(NS)))

  for (s in seq_len(NSEED)) {
    ## the seed rule of record, unchanged, and the benchmark's own draw
    set.seed(20260627 + 1000 * s + round(GRID_D * 100))
    y_full <- GRID_D * round(mw_r(N_DRAW, d) / GRID_D)

    for (j in seq_along(NS)) {
      y <- y_full[seq_len(NS[j])]                       # nested: the leading n
      fc <- adkde(y, GRID_D, grid)
      RHO[s, j]  <- attr(fc, "rho")
      PICK[s, j] <- attr(fc, "pick")
      ## the benchmark's own diagnostic: the grid against the reference bandwidth. It is
      ## recorded because it is the only quantity in this experiment that moves with n by
      ## construction, through the n^(-1/5) in the bandwidth, and it is what decides
      ## whether the uncorrected kernel is wide enough to smooth the lattice away.
      DOH[s, j] <- GRID_D / (1.06 * stats::sd(y) * length(y)^(-1 / 5))
      E[s, "naive",    j] <- ise(naive_kde(y, grid), ft)
      E[s, "deconv",   j] <- ise(deconv_kde(y, GRID_D, grid), ft)
      E[s, "combined", j] <- ise(as.numeric(fc), ft)
    }
  }

  for (j in seq_along(NS)) {
    mn  <- apply(E[, , j, drop = FALSE], 2, function(v) mean(v, na.rm = TRUE)) * 1e3
    sdv <- apply(E[, , j, drop = FALSE], 2, function(v) stats::sd(v, na.rm = TRUE)) * 1e3
    adv <- (E[, "naive", j] - E[, "combined", j]) * 1e3       # paired, positive is a gain
    advd <- (E[, "deconv", j] - E[, "combined", j]) * 1e3
    pk  <- table(PICK[, j])

    rows[[sprintf("%s %s D%.2f n%d", d$index, d$name, GRID_D, NS[j])]] <- list(
      index = as.character(d$index), name = d$name, D = GRID_D,
      n = NS[j], n_draw = N_DRAW, nested = TRUE, nseed = NSEED, M = M, dx = dx,
      mean_x1e3 = round(mn, 6), sd_x1e3 = round(sdv, 6),
      perseed_x1e3 = setNames(lapply(METHODS, function(k) round(E[, k, j] * 1e3, 6)), METHODS),
      ratio_combined_naive = mn["combined"] / mn["naive"],
      ratio_deconv_naive   = mn["deconv"]   / mn["naive"],
      advantage_vs_naive_x1e3  = c(mean = mean(adv),  se = stats::sd(adv)  / sqrt(length(adv))),
      advantage_vs_deconv_x1e3 = c(mean = mean(advd), se = stats::sd(advd) / sqrt(length(advd))),
      wins_vs_naive  = sum(adv  > 0), wins_vs_deconv = sum(advd > 0),
      rho_mean = mean(RHO[, j], na.rm = TRUE),
      D_over_h_mean = mean(DOH[, j], na.rm = TRUE),
      pick_counts = setNames(as.integer(pk), names(pk)))

    cat(sprintf("%-10s %6d %10.4f %10.4f %10.4f %8.4f %10.4f %7.4f %4d/%-3d %5.2f\n",
                substr(d$name, 1, 10), NS[j], mn["naive"], mn["deconv"], mn["combined"],
                mn["combined"] / mn["naive"], mean(adv),
                stats::sd(adv) / sqrt(length(adv)), sum(adv > 0), NSEED,
                mean(DOH[, j], na.rm = TRUE)))
  }
  cat("\n")
}

## ---- provenance, written before the check so a failure costs the diagnosis only --------
prov <- list(
  script = "sensitivity_n_v13.R",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  elapsed_sec = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
  r_version = R.version.string, platform = R.version$platform,
  pkg_dir = normalizePath(pkg), mw15_file = normalizePath(mwfile),
  benchmark_script = normalizePath(benchR),
  definitions_lifted = "TABLE1,M,grid,dx,n,ise",
  targets = paste(TARGETS, collapse = ","),
  ns = paste(NS, collapse = ","), n_draw = N_DRAW, nested_draws = TRUE,
  D = GRID_D, nseed = NSEED, M = M, dx = dx, grid_lo = -10, grid_hi = 10,
  methods = paste(METHODS, collapse = ","),
  seed_rule = "20260627 + 1000*replicate + round(D*100)",
  smoke = SMOKE)

## ---- the Table 1 anchor -----------------------------------------------------------------
## Reproduce the published cell. The per-seed comparison is the strong form and runs at any
## seed count; the mean comparison is the published form and runs only at the stored count.
chk <- list(ran = FALSE, pass = NA, note = "", detail = NULL)
if (!(n_bench %in% NS)) {
  chk$note <- sprintf("skipped: NS does not include the benchmark's n=%d", n_bench)
} else if (!file.exists(BENCH_RDS)) {
  chk$note <- paste("skipped: no stored benchmark at", BENCH_RDS)
} else if (abs(GRID_D - 0.5) > 1e-12) {
  chk$note <- sprintf("skipped: the stored cell used for the anchor is D=0.50, this run is D=%.2f", GRID_D)
} else {
  bch <- readRDS(BENCH_RDS)
  det <- list(); ok <- TRUE
  for (ti in TARGETS) {
    d <- TABLE1[[ti]]
    key <- sprintf("table1 %s %s %s", d$index, d$name, GRID_D)
    if (is.null(bch$rows[[key]])) {
      ok <- FALSE
      det[[ti]] <- list(status = "missing", key = key)
      next
    }
    br <- bch$rows[[key]]
    me <- rows[[sprintf("%s %s D%.2f n%d", d$index, d$name, GRID_D, n_bench)]]
    ## per-seed, over the seeds this run actually computed
    ps <- sapply(METHODS, function(k) {
      a <- me$perseed_x1e3[[k]]
      b <- br$perseed_x1e3[[k]][seq_len(NSEED)]
      if (length(b) < NSEED || all(is.na(b))) return(NA_real_)
      max(abs(a - b), na.rm = TRUE)
    })
    ## means, only when this run used the stored replicate count
    mm <- if (identical(as.integer(NSEED), as.integer(br$nseed)))
      sapply(METHODS, function(k) abs(me$mean_x1e3[[k]] - br$mean_x1e3[[k]])) else
      setNames(rep(NA_real_, length(METHODS)), METHODS)
    bad_ps <- any(!is.na(ps) & ps > CHECK_TOL)
    bad_mm <- any(!is.na(mm) & mm > 5e-4)
    if (bad_ps || bad_mm) ok <- FALSE
    det[[ti]] <- list(status = if (bad_ps || bad_mm) "FAIL" else "pass",
                      key = key, nseed_stored = br$nseed, nseed_run = NSEED,
                      max_abs_perseed_x1e3 = round(ps, 9),
                      abs_mean_diff_x1e3 = round(mm, 9),
                      stored_mean_x1e3 = sapply(METHODS, function(k) br$mean_x1e3[[k]]),
                      run_mean_x1e3 = sapply(METHODS, function(k) me$mean_x1e3[[k]]))
  }
  chk <- list(ran = TRUE, pass = ok, tol_perseed_x1e3 = CHECK_TOL,
              tol_mean_x1e3 = 5e-4, note = "", detail = det)
}
prov$table1_anchor <- chk

saveRDS(list(provenance = prov, rows = rows), file.path(OUT, "sensitivity_n_v13.rds"))

## minimal JSON writer, base R only, same shape as benchmark_v10_2.json
esc  <- function(s) gsub('"', '\\\\"', as.character(s))
jnum <- function(v) if (length(v) == 0 || is.na(v)) "null" else format(v, scientific = FALSE, trim = TRUE)
jarr <- function(v) paste0("[", paste(sapply(v, jnum), collapse = ","), "]")
jstr <- function(s) paste0('"', esc(s), '"')
kv   <- function(k, v) paste0('"', esc(k), '":', v)
jobj <- function(nms, vals) paste0("{", paste(mapply(kv, nms, vals), collapse = ","), "}")
provj <- jobj(names(prov)[names(prov) != "table1_anchor"],
              sapply(prov[names(prov) != "table1_anchor"], function(v)
                if (is.numeric(v)) jnum(v) else if (is.logical(v)) tolower(as.character(v))
                else jstr(v)))
rowj <- sapply(names(rows), function(rn) {
  r <- rows[[rn]]
  paste0(jstr(rn), ":", jobj(
    c("index", "name", "D", "n", "n_draw", "nseed", "mean_x1e3", "sd_x1e3",
      "perseed_x1e3", "ratio_combined_naive", "ratio_deconv_naive",
      "advantage_vs_naive_x1e3", "advantage_vs_deconv_x1e3",
      "wins_vs_naive", "wins_vs_deconv", "rho_mean", "D_over_h_mean", "pick_counts"),
    c(jstr(r$index), jstr(r$name), jnum(r$D), jnum(r$n), jnum(r$n_draw), jnum(r$nseed),
      jobj(names(r$mean_x1e3), sapply(r$mean_x1e3, jnum)),
      jobj(names(r$sd_x1e3),   sapply(r$sd_x1e3, jnum)),
      jobj(names(r$perseed_x1e3), sapply(r$perseed_x1e3, jarr)),
      jnum(r$ratio_combined_naive), jnum(r$ratio_deconv_naive),
      paste0('{"mean":', jnum(r$advantage_vs_naive_x1e3["mean"]),
             ',"se":',   jnum(r$advantage_vs_naive_x1e3["se"]), "}"),
      paste0('{"mean":', jnum(r$advantage_vs_deconv_x1e3["mean"]),
             ',"se":',   jnum(r$advantage_vs_deconv_x1e3["se"]), "}"),
      jnum(r$wins_vs_naive), jnum(r$wins_vs_deconv), jnum(r$rho_mean),
      jnum(r$D_over_h_mean),
      jobj(names(r$pick_counts), sapply(r$pick_counts, jnum)))))
})
writeLines(paste0('{"provenance":', provj, ',"rows":{', paste(rowj, collapse = ","), "}}"),
           file.path(OUT, "sensitivity_n_v13.json"))

cat(sprintf("-> %s\n-> %s\nelapsed %.1f s over %d cells\n\n",
            file.path(OUT, "sensitivity_n_v13.rds"),
            file.path(OUT, "sensitivity_n_v13.json"),
            prov$elapsed_sec, length(rows)))

## ---- report the anchor, and fail loudly ------------------------------------------------
cat("== Table 1 anchor, n = 4000 against results/benchmark_v10_2.rds ==\n")
if (!chk$ran) {
  cat(sprintf("   NOT RUN: %s\n", chk$note))
} else {
  for (ti in names(chk$detail)) {
    dd <- chk$detail[[ti]]
    if (identical(dd$status, "missing")) {
      cat(sprintf("   %-3s FAIL  stored benchmark has no row keyed '%s'\n", ti, dd$key))
      next
    }
    cat(sprintf("   %-3s %-4s  seeds run %d, stored %d\n", ti, dd$status,
                dd$nseed_run, dd$nseed_stored))
    for (k in METHODS)
      cat(sprintf("        %-9s per-seed max |diff| %-12s  mean stored %-11s run %-11s |diff| %s\n",
                  k,
                  if (is.na(dd$max_abs_perseed_x1e3[[k]])) "n/a" else
                    format(dd$max_abs_perseed_x1e3[[k]], scientific = TRUE, digits = 3),
                  format(round(dd$stored_mean_x1e3[[k]], 4), nsmall = 4),
                  format(round(dd$run_mean_x1e3[[k]], 4), nsmall = 4),
                  if (is.na(dd$abs_mean_diff_x1e3[[k]])) "n/a (seed counts differ)" else
                    format(dd$abs_mean_diff_x1e3[[k]], scientific = TRUE, digits = 3)))
  }
  cat(sprintf("   tolerance: per-seed %g, mean %g, both on ISE x 1e3\n",
              chk$tol_perseed_x1e3, chk$tol_mean_x1e3))
}
cat("\n")

if (isTRUE(chk$ran) && !isTRUE(chk$pass)) {
  cat("********************************************************************************\n")
  cat("* The n = 4000 row does not reproduce the published Table 1 cell.              *\n")
  cat("* The result files above were written so this can be diagnosed, but no number  *\n")
  cat("* in them should be used. Send back this log.                                  *\n")
  cat("********************************************************************************\n")
  quit(status = 1L)
}
if (!isTRUE(chk$ran))
  cat("note: the anchor did not run, so nothing here has been checked against Table 1.\n")
