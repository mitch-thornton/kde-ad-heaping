#!/usr/bin/env Rscript
## benchmark_v10_2.R
## Supersedes benchmark_v10.R. Produces Table 1 and the full Marron-Wand supplementary
## table from one pass, so the two cannot disagree on the densities they share.
##
## Changes from the shipped benchmark.R:
##   * seed count raised from 8 to NSEED (default 50), Reviewer 2 item E5
##   * paired per-seed differences and paired standard errors emitted, so the tie
##     decisions in the table can be checked from the shipped artifact
##   * the full fifteen-density Marron-Wand battery added for supplementary, C11
##   * three candidate abstention rules evaluated alongside the unmodified estimator,
##     so the rule for Reviewer 2 item 6 is chosen from evidence
##   * a provenance block written into the results file
##
## Seeds are generated, not listed. The seed for replicate s at grid D is
##     20260627 + 1000*s + round(D*100)
## which is the rule benchmark.R already uses, anchored on the seed of record 20260627.
##
## Usage:
##   PKG_DIR=/path/to/kde-ad-heaping NSEED=50 SET=both Rscript benchmark_v10_2.R
##   PKG_DIR=... NSEED=2 SMOKE=1 Rscript benchmark_v10_2.R
##
## Environment:
##   PKG_DIR  checkout of github.com/mitch-thornton/kde-ad-heaping   (required)
##   MW15_R   path to mw15.R                       (default: beside this script)
##   SET      table1 | mw15 | both                 (default both)
##   NSEED    replicates per cell                  (default 50)
##   OUTDIR   results directory                    (default $PKG_DIR/results)
##   NO_SEM   set to 1 to skip the Kernelheaping SEM baseline, which dominates runtime
##   DELTA    taper-rule abstention constant       (default 0.05)
##   KAPPA    scale-rule abstention constant       (default 0.5)
##   SMOKE    set to 1 for a fast reduced run

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

NSEED  <- as.integer(Sys.getenv("NSEED", "50"))
SET    <- Sys.getenv("SET", "both")
SMOKE  <- nzchar(Sys.getenv("SMOKE", ""))
NO_SEM <- nzchar(Sys.getenv("NO_SEM", ""))
DELTA  <- as.numeric(Sys.getenv("DELTA", "0.05"))
KAPPA  <- as.numeric(Sys.getenv("KAPPA", "0.5"))
OUT    <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

has_km <- (!NO_SEM) && requireNamespace("Kernelheaping", quietly = TRUE)

## ---- the four densities currently in Table 1 -------------------------------
## Verified against the Marron-Wand definitions by verify_mw15.R. Numbers 1 and 3 below
## are Marron-Wand members 1 and 4 exactly. Numbers 2 and 4 are not Marron-Wand
## densities, so the paper text that attributes them to that benchmark needs correcting.
TABLE1 <- list(
  list(index = "T1", name = "Gaussian",       mw = 1L,
       w = 1, mu = 0, sd = 1),
  list(index = "T2", name = "Bimodal",        mw = NA_integer_,
       w = c(.5, .5), mu = c(-1.2, 1.2), sd = c(.5, .5)),
  list(index = "T3", name = "Kurtotic",       mw = 4L,
       w = c(2/3, 1/3), mu = c(0, 0), sd = c(1, .1)),
  list(index = "T4", name = "Strongly skewed (paper variant)", mw = NA_integer_,
       w = rep(.2, 5), mu = c(0, .5, 1.0833, 1.4167, 1.6875),
       sd = c(1, 2/3, 4/9, 8/27, 16/81)))

SETS <- list()
if (SET %in% c("table1", "both")) SETS$table1 <- TABLE1
if (SET %in% c("mw15",   "both")) SETS$mw15   <- MW15
if (!length(SETS)) stop("SET must be table1, mw15 or both")

M <- 2048; grid <- seq(-10, 10, length.out = M + 1)[-(M + 1)]; dx <- grid[2] - grid[1]
n <- 4000; Ds <- c(0.25, 0.5, 1.0, 1.5)
if (SMOKE) { Ds <- c(0.25, 1.0); SETS <- lapply(SETS, function(s) s[seq_len(min(3, length(s)))]) }
ise <- function(fh, ft) sum((fh - ft)^2) * dx
sincf <- function(u) ifelse(abs(u) < 1e-12, 1, sin(u)/u)

## ---- candidate abstention statistics ---------------------------------------
## Rule C, the taper rule. Let w* be the frequency at which the empirical power spectrum
## descends to the residue floor. The correction divides by sinc(wD/2). When that taper
## stays near unity out to w*, the correction has nothing to do and only adds variance.
.taper_at_wstar <- function(y, D, grid) {
  Mg <- length(grid); b <- bin_prob(y, grid)
  S <- Mod(fft(b$p))^2; w <- .fftw(Mg, b$dx)
  fl <- .residue_floor(S, w)
  half <- 2:(Mg %/% 2); wp <- w[half]; Sp <- S[half]
  hit <- which(Sp <= fl)
  wstar <- if (length(hit)) wp[hit[1]] else max(wp)
  sincf(wstar * D / 2)
}

METHODS <- c("naive", "deconv", "sem", "imput", "deheap", "super",
             "combined", "abst_det", "abst_scale", "abst_taper")

cat(sprintf("benchmark_v10_2  n=%d seeds=%d M=%d dx=%.6f sets=%s SEM=%s delta=%.2f kappa=%.2f\n",
            n, NSEED, M, dx, paste(names(SETS), collapse = "+"), has_km, DELTA, KAPPA))
cat(sprintf("%-8s %2s %-26s %5s %9s %9s %9s %9s %9s %6s\n",
            "set", "#", "density", "D", "naive", "combined", "abst-det", "abst-scl", "abst-tpr", "fire"))

rows <- list()
for (sn in names(SETS)) {
  for (d in SETS[[sn]]) {
    ft <- mw_d(grid, d)
    for (D in Ds) {
      e <- setNames(lapply(METHODS, function(k) rep(NA_real_, NSEED)), METHODS)
      fired <- rep(NA, NSEED); tap <- rep(NA_real_, NSEED)
      rho <- rep(NA_real_, NSEED); dovh <- rep(NA_real_, NSEED)

      for (s in seq_len(NSEED)) {
        set.seed(20260627 + 1000 * s + round(D * 100))
        y <- D * round(mw_r(n, d) / D)

        f_naive <- naive_kde(y, grid)
        f_comb  <- adkde(y, D, grid)
        rho[s]  <- attr(f_comb, "rho"); f_comb <- as.numeric(f_comb)

        det <- tryCatch(heap_detect(y), error = function(x) list(detected = FALSE))
        fired[s] <- isTRUE(det$detected)
        h <- 1.06 * stats::sd(y) * n^(-1/5); dovh[s] <- D / h
        tap[s] <- tryCatch(.taper_at_wstar(y, D, grid), error = function(x) NA_real_)

        e$naive[s]    <- ise(f_naive, ft)
        e$deconv[s]   <- ise(deconv_kde(y, D, grid), ft)
        e$imput[s]    <- ise(heitjan_mi(y, D, grid, M = 8), ft)
        e$deheap[s]   <- ise(deheap_kde(y, D, grid), ft)
        e$super[s]    <- ise(superpose_kde(y, D, grid), ft)
        e$combined[s] <- ise(f_comb, ft)
        ## rule A, Reviewer 2 as written: abstain when the comb detector does not fire
        e$abst_det[s]   <- if (fired[s]) e$combined[s] else e$naive[s]
        ## rule B, scale: abstain when the grid is small against the kernel bandwidth
        e$abst_scale[s] <- if (D >= KAPPA * h) e$combined[s] else e$naive[s]
        ## rule C, taper: abstain when the box taper is nearly flat over the informative band
        e$abst_taper[s] <- if (!is.na(tap[s]) && tap[s] > 1 - DELTA) e$naive[s] else e$combined[s]
        if (has_km) { fs <- sem_kde(y, D, grid); if (!is.null(fs)) e$sem[s] <- ise(fs, ft) }
      }

      mn  <- sapply(e, function(v) mean(v, na.rm = TRUE))
      sdv <- sapply(e, function(v) stats::sd(v, na.rm = TRUE))
      pair <- lapply(METHODS, function(k) {
        dd <- e[[k]] - e$combined
        if (all(is.na(dd))) return(c(mean = NA_real_, se = NA_real_))
        c(mean = mean(dd, na.rm = TRUE),
          se = stats::sd(dd, na.rm = TRUE) / sqrt(sum(!is.na(dd))))
      })
      names(pair) <- METHODS

      key <- paste(sn, d$index, d$name, D)
      rows[[key]] <- list(
        set = sn, index = as.character(d$index), name = d$name,
        mw_member = if (is.null(d$mw)) d$index else d$mw,
        D = D, n = n, nseed = NSEED, M = M, dx = dx,
        mean_x1e3 = round(mn * 1e3, 5), sd_x1e3 = round(sdv * 1e3, 5),
        paired_vs_combined_x1e3 = lapply(pair, function(p) round(p * 1e3, 5)),
        perseed_x1e3 = lapply(e, function(v) round(v * 1e3, 6)),
        detector_fire_rate = mean(fired, na.rm = TRUE),
        taper_mean = mean(tap, na.rm = TRUE),
        D_over_h_mean = mean(dovh, na.rm = TRUE),
        rho_mean = mean(rho, na.rm = TRUE))

      cat(sprintf("%-8s %2s %-26s %5.2f %9.3f %9.3f %9.3f %9.3f %9.3f %6.2f\n",
                  sn, as.character(d$index), substr(d$name, 1, 26), D,
                  mn["naive"]*1e3, mn["combined"]*1e3, mn["abst_det"]*1e3,
                  mn["abst_scale"]*1e3, mn["abst_taper"]*1e3, mean(fired, na.rm = TRUE)))
    }
  }
}

prov <- list(
  script = "benchmark_v10_2.R",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  elapsed_sec = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
  r_version = R.version.string, platform = R.version$platform,
  pkg_dir = normalizePath(pkg), mw15_file = normalizePath(mwfile),
  set = SET, nseed = NSEED, n = n, M = M, dx = dx, grid_lo = -10, grid_hi = 10,
  Ds = paste(Ds, collapse = ","),
  seed_rule = "20260627 + 1000*replicate + round(D*100)",
  delta = DELTA, kappa = KAPPA, kernelheaping = has_km, smoke = SMOKE)

saveRDS(list(provenance = prov, rows = rows), file.path(OUT, "benchmark_v10_2.rds"))

## minimal JSON writer, base R only
esc <- function(s) gsub('"', '\\\\"', as.character(s))
jnum <- function(v) if (length(v) == 0 || is.na(v)) "null" else format(v, scientific = FALSE, trim = TRUE)
jarr <- function(v) paste0("[", paste(sapply(v, jnum), collapse = ","), "]")
kv <- function(k, v) paste0('"', esc(k), '":', v)
jobj <- function(nms, vals) paste0("{", paste(mapply(kv, nms, vals), collapse = ","), "}")
provj <- jobj(names(prov), sapply(prov, function(v)
  if (is.numeric(v)) jnum(v) else if (is.logical(v)) tolower(as.character(v))
  else paste0('"', esc(v), '"')))
rowj <- sapply(names(rows), function(rn) {
  r <- rows[[rn]]
  ms <- jobj(names(r$mean_x1e3), sapply(r$mean_x1e3, jnum))
  ss <- jobj(names(r$sd_x1e3), sapply(r$sd_x1e3, jnum))
  pp <- jobj(names(r$paired_vs_combined_x1e3),
             sapply(r$paired_vs_combined_x1e3, function(p)
               paste0('{"mean":', jnum(p["mean"]), ',"se":', jnum(p["se"]), "}")))
  ps <- jobj(names(r$perseed_x1e3), sapply(r$perseed_x1e3, jarr))
  paste0('"', esc(rn), '":', jobj(
    c("set","index","name","mw_member","D","mean_x1e3","sd_x1e3",
      "paired_vs_combined_x1e3","perseed_x1e3","detector_fire_rate",
      "taper_mean","D_over_h_mean","rho_mean"),
    c(paste0('"', esc(r$set), '"'), paste0('"', esc(r$index), '"'),
      paste0('"', esc(r$name), '"'), jnum(r$mw_member), jnum(r$D),
      ms, ss, pp, ps, jnum(r$detector_fire_rate), jnum(r$taper_mean),
      jnum(r$D_over_h_mean), jnum(r$rho_mean))))
})
writeLines(paste0('{"provenance":', provj, ',"rows":{', paste(rowj, collapse = ","), "}}"),
           file.path(OUT, "benchmark_v10_2.json"))

cat(sprintf("\n-> %s\n-> %s\nelapsed %.1f s over %d cells\n",
            file.path(OUT, "benchmark_v10_2.rds"), file.path(OUT, "benchmark_v10_2.json"),
            prov$elapsed_sec, length(rows)))
