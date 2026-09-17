#!/usr/bin/env Rscript
## endtoend_v10_4.R
## The two substantive reviewer items still open after v10.4, in one pass. They share the
## targets, the grid, the seed rule and the package, so they share a script.
##
##   PART A  End-to-end grid detection (work item B12; Reviewer 1 item 3; editor item
##           "end-to-end automatic grid detection"). The grid width is estimated from the
##           data with no knowledge of the true value, the estimate is fed to the combined
##           estimator, and the integrated squared error of that pipeline is compared
##           against the same estimator handed the true grid. Four pipelines are run in
##           every cell so the reader's contribution is separable from the estimator's:
##             oracle     the true D, which is what Table 1 reports
##             verifying  D from heap_grid with the true grid supplied as the hint
##             blind      D from heap_grid with no hint
##             detector   D from heap_detect, in the replicates where it fires
##           The naive kernel is carried as the do-nothing reference, so the question
##           "is an automatic pipeline better than leaving the data alone" is answered
##           directly rather than by inference.
##
##   PART B  Known-truth test of the heaped-fraction reader (work item B13, the half of
##           Reviewer 1 item 5 that the mixed-grain sweep does not cover). The manuscript
##           states that fractions of 0.2, 0.4 and 0.6 are recovered as 0.175, 0.334 and
##           0.506, "so p is identifiable up to a stable calibratable attenuation". That
##           sentence carries three claims and this part tests all three: that recovery is
##           monotone in p, that it is near-linear, and that the attenuation is STABLE,
##           meaning one calibration serves across targets and grids. The third is the one
##           that makes "calibratable" mean anything, because a calibration that depends on
##           the unknown density is not available to a practitioner. Also reported is what
##           the reader returns at p = 0, where there is no comb at all.
##
## Both parts are indicative until run on your machine. Nothing here is typed into the
## paper; the run writes endtoend_v10_4.rds and the table emitter reads it.
##
## Usage (bash, not zsh):
##     PKG_DIR=/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping bash run_endtoend_v10_4.sh
##
## Environment:
##   PKG_DIR  your kde-ad-heaping checkout                      (required)
##   MW15_R   path to mw15.R                                    (default: beside this script)
##   NSEED    replicates per cell, part A                       (default 20)
##   NSEEDB   replicates per cell, part B                       (default 20)
##   PART     A | B | both                                      (default both)
##   OUTDIR   results directory                                 (default $PKG_DIR/results)
##   SMOKE    set to 1 for a fast reduced run
##
## Seed rules, both generated rather than listed:
##   part A   20260627 + 1000*s + round(100*D)                  (the rule Table 1 uses)
##   part B   20260627 + 1000*s + round(100*D) + 10000*round(10*p)

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

NSEED  <- as.integer(Sys.getenv("NSEED",  "20"))
NSEEDB <- as.integer(Sys.getenv("NSEEDB", "20"))
PART   <- Sys.getenv("PART", "both")
SMOKE  <- nzchar(Sys.getenv("SMOKE", ""))
OUT    <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- the four Table 1 targets, identical to benchmark_v10_2.R --------------
TABLE1 <- list(
  list(index = "T1", name = "Gaussian",       mw = 1L,  w = 1, mu = 0, sd = 1),
  list(index = "T2", name = "Bimodal",        mw = NA_integer_,
       w = c(.5, .5), mu = c(-1.2, 1.2), sd = c(.5, .5)),
  list(index = "T3", name = "Kurtotic",       mw = 4L,  w = c(2/3, 1/3), mu = c(0, 0), sd = c(1, .1)),
  list(index = "T4", name = "Skewed",         mw = NA_integer_,
       w = rep(.2, 5), mu = c(0, .5, 1.0833, 1.4167, 1.6875),
       sd = c(1, 2/3, 4/9, 8/27, 16/81)))

M <- 2048; grid <- seq(-10, 10, length.out = M + 1)[-(M + 1)]; dx <- grid[2] - grid[1]
n <- 4000
DS_A <- c(0.25, 0.5, 1.0, 1.5)                 # part A, the Table 1 grids
DS_B <- c(0.5, 1.0, 2.0)                       # part B, the grids of the partial-heaping claim
PS   <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
if (SMOKE) { NSEED <- 2; NSEEDB <- 3; DS_A <- c(0.5, 1.0); DS_B <- c(1.0); TABLE1 <- TABLE1[1:2] }
ise <- function(fh, ft) sum((fh - ft)^2) * dx
hr  <- function(t) cat("\n", strrep("=", 78), "\n ", t, "\n", strrep("=", 78), "\n", sep = "")

## A grid estimate is usable when it is finite, resolvable on the analysis grid, and not
## the span artifact the blind reader returns when it finds no tooth. Estimates outside
## this range are counted, not quietly replaced, so the failure rate is visible.
D_LO <- 2 * dx; D_HI <- (max(grid) - min(grid)) / 5
usable <- function(D) is.finite(D) && D > D_LO && D <= D_HI

res <- list()
cat(sprintf("endtoend_v10_4  n=%d M=%d dx=%.6f partA_seeds=%d partB_seeds=%d part=%s\n",
            n, M, dx, NSEED, NSEEDB, PART))
cat(sprintf("usable grid estimate: %.4f < Dhat <= %.4f\n", D_LO, D_HI))

## ============================================================================
## PART A. End-to-end grid detection
## ============================================================================
if (PART %in% c("A", "both")) {
  hr("A  end-to-end: estimate D from the data, then reconstruct with the estimate")
  cat(sprintf("  %-9s %5s | %9s %9s %9s %9s %9s | %s\n",
              "target", "D", "naive", "oracle", "verifying", "blind", "detector", "blind Dhat"))
  cells <- list()
  for (d in TABLE1) {
    ft <- mw_d(grid, d)
    for (D in DS_A) {
      e <- list(naive = rep(NA_real_, NSEED), oracle = rep(NA_real_, NSEED),
                verifying = rep(NA_real_, NSEED), blind = rep(NA_real_, NSEED),
                detector = rep(NA_real_, NSEED))
      Dh <- list(verifying = rep(NA_real_, NSEED), blind = rep(NA_real_, NSEED),
                 detector = rep(NA_real_, NSEED))
      ok <- list(verifying = rep(FALSE, NSEED), blind = rep(FALSE, NSEED),
                 detector = rep(FALSE, NSEED))
      fired <- rep(FALSE, NSEED)

      for (s in seq_len(NSEED)) {
        set.seed(20260627 + 1000 * s + round(D * 100))
        y <- D * round(mw_r(n, d) / D)

        e$naive[s]  <- ise(naive_kde(y, grid), ft)
        e$oracle[s] <- ise(as.numeric(adkde(y, D, grid)), ft)

        Dh$verifying[s] <- tryCatch(heap_grid(y, grid, near = D), error = function(x) NA_real_)
        Dh$blind[s]     <- tryCatch(heap_grid(y, grid),           error = function(x) NA_real_)
        det <- tryCatch(heap_detect(y), error = function(x) list(detected = FALSE, D_hat = NA_real_))
        fired[s] <- isTRUE(det$detected)
        Dh$detector[s] <- if (fired[s]) det$D_hat else NA_real_

        ## reconstruct from each estimate. An unusable estimate is recorded as a failure
        ## of the pipeline rather than replaced by the true value.
        for (k in c("verifying", "blind", "detector")) {
          Dk <- Dh[[k]][s]
          if (usable(Dk)) {
            v <- tryCatch(ise(as.numeric(adkde(y, Dk, grid)), ft), error = function(x) NA_real_)
            if (is.finite(v)) { e[[k]][s] <- v; ok[[k]][s] <- TRUE }
          }
        }
      }

      mn <- sapply(e, function(v) mean(v, na.rm = TRUE))
      key <- paste(d$name, D)
      cells[[key]] <- list(
        target = d$name, D = D, nseed = NSEED, ise = e, Dhat = Dh, usable = ok,
        fire_rate = mean(fired),
        mean_ise = mn,
        n_ok = sapply(ok, sum),
        blind_Dhat_mean = mean(Dh$blind, na.rm = TRUE),
        blind_relerr = mean(Dh$blind / D - 1, na.rm = TRUE),
        ## per-seed paired comparisons, so the claims are checkable
        paired_blind_vs_oracle = { dd <- e$blind - e$oracle; dd <- dd[is.finite(dd)]
          if (length(dd) > 1) c(mean = mean(dd), se = stats::sd(dd)/sqrt(length(dd))) else c(mean = NA_real_, se = NA_real_) },
        paired_blind_vs_naive  = { dd <- e$blind - e$naive;  dd <- dd[is.finite(dd)]
          if (length(dd) > 1) c(mean = mean(dd), se = stats::sd(dd)/sqrt(length(dd))) else c(mean = NA_real_, se = NA_real_) },
        paired_verif_vs_oracle = { dd <- e$verifying - e$oracle; dd <- dd[is.finite(dd)]
          if (length(dd) > 1) c(mean = mean(dd), se = stats::sd(dd)/sqrt(length(dd))) else c(mean = NA_real_, se = NA_real_) })

      f1 <- function(v) if (!is.finite(v)) "      --" else if (v >= 100) sprintf("%9.1f", v) else sprintf("%9.3f", v)
      cat(sprintf("  %-9s %5.2f | %s %s %s %s %s | %6.3f (%d/%d ok)\n",
                  d$name, D, f1(1e3*mn["naive"]), f1(1e3*mn["oracle"]), f1(1e3*mn["verifying"]),
                  f1(1e3*mn["blind"]), f1(1e3*mn["detector"]),
                  mean(Dh$blind, na.rm = TRUE), sum(ok$blind), NSEED))
    }
  }

  ## ---- summary that answers the reviewer question directly -----------------
  nm <- names(cells)
  gm <- function(k) sapply(cells, function(c) c$mean_ise[[k]])
  nok <- function(k) sapply(cells, function(c) c$n_ok[[k]])
  tot <- length(cells) * NSEED
  blind_ok    <- sum(nok("blind"));  verif_ok <- sum(nok("verifying")); det_ok <- sum(nok("detector"))
  ## a cell counts as recovered when the blind pipeline is within 10 percent of the oracle
  within10 <- sum(is.finite(gm("blind")) & gm("blind") <= 1.10 * gm("oracle"))
  beat_naive_blind  <- sum(is.finite(gm("blind"))  & gm("blind")  < gm("naive"))
  beat_naive_oracle <- sum(gm("oracle") < gm("naive"))
  cat("\n")
  cat(sprintf("  usable grid estimates, of %d replicate-cells : verifying %d, blind %d, detector %d\n",
              tot, verif_ok, blind_ok, det_ok))
  cat(sprintf("  cells where the verifying pipeline matches the oracle to 3 decimals : %d of %d\n",
              sum(abs(gm("verifying") - gm("oracle")) < 5e-4, na.rm = TRUE), length(cells)))
  cat(sprintf("  cells where the blind pipeline is within 10%% of the oracle          : %d of %d\n",
              within10, length(cells)))
  cat(sprintf("  cells where the blind pipeline beats the naive kernel               : %d of %d\n",
              beat_naive_blind, length(cells)))
  cat(sprintf("  cells where the oracle pipeline beats the naive kernel              : %d of %d\n",
              beat_naive_oracle, length(cells)))
  cat("\n  READ THIS AS: the difference between the oracle and verifying columns is the cost\n")
  cat("  of reading a grid you already know. The difference between oracle and blind is the\n")
  cat("  cost of not knowing it, which is what Reviewer 1 item 3 asks for. If the blind\n")
  cat("  column is worse than naive in most cells, the end-to-end pipeline is not usable\n")
  cat("  as an automatic procedure and the paper should say so with these numbers.\n")

  res$partA <- list(cells = cells, n = n, M = M, dx = dx, nseed = NSEED,
                    D_lo = D_LO, D_hi = D_HI,
                    summary = list(total_replicate_cells = tot, usable_verifying = verif_ok,
                                   usable_blind = blind_ok, usable_detector = det_ok,
                                   cells = length(cells), blind_within10 = within10,
                                   blind_beats_naive = beat_naive_blind,
                                   oracle_beats_naive = beat_naive_oracle))
}

## ============================================================================
## PART B. Heaped-fraction reader against known truth
## ============================================================================
if (PART %in% c("B", "both")) {
  hr("B  heaped-fraction reader against known truth")
  cat(sprintf("  n = %d, %d replicates per cell, p in {%s}\n", n, NSEEDB, paste(PS, collapse = ", ")))
  cat("  the reader is given the true D, so this tests recovery of p alone\n\n")
  cat(sprintf("  %-9s %5s |%s| %6s %6s %6s\n", "target", "D",
              paste(sprintf("%7.1f", PS), collapse = ""), "slope", "intcpt", "R2"))
  pcells <- list()
  for (d in TABLE1) {
    for (D in DS_B) {
      P <- matrix(NA_real_, NSEEDB, length(PS), dimnames = list(NULL, as.character(PS)))
      for (j in seq_along(PS)) {
        p <- PS[j]
        for (s in seq_len(NSEEDB)) {
          set.seed(20260627 + 1000 * s + round(D * 100) + 10000 * round(10 * p))
          x <- mw_r(n, d)
          hp <- stats::runif(n) < p
          y <- x; y[hp] <- D * round(x[hp] / D)
          P[s, j] <- tryCatch(heap_fraction(y, D, grid), error = function(e) NA_real_)
        }
      }
      m <- colMeans(P, na.rm = TRUE); se <- apply(P, 2, function(v) stats::sd(v, na.rm = TRUE)) / sqrt(NSEEDB)
      fit <- stats::lm(m ~ PS)
      r2  <- summary(fit)$r.squared
      mono <- all(diff(m) > -1e-9)
      key <- paste(d$name, D)
      pcells[[key]] <- list(target = d$name, D = D, p = PS, phat_mean = m, phat_se = se,
                            phat = P, slope = unname(coef(fit)[2]), intercept = unname(coef(fit)[1]),
                            r2 = r2, monotone = mono, at_zero = m[1], nseed = NSEEDB)
      cat(sprintf("  %-9s %5.2f |%s| %6.3f %6.3f %6.3f%s\n", d$name, D,
                  paste(sprintf("%7.3f", m), collapse = ""),
                  coef(fit)[2], coef(fit)[1], r2, if (mono) "" else "  NOT MONOTONE"))
    }
  }

  sl <- sapply(pcells, function(c) c$slope); ic <- sapply(pcells, function(c) c$intercept)
  z0 <- sapply(pcells, function(c) c$at_zero); r2 <- sapply(pcells, function(c) c$r2)
  mono_all <- all(sapply(pcells, function(c) c$monotone))
  ## one calibration for all cells, then the worst error it leaves behind
  pooled <- stats::lm(as.vector(sapply(pcells, function(c) c$phat_mean)) ~ rep(PS, length(pcells)))
  a <- unname(coef(pooled)[2]); b <- unname(coef(pooled)[1])
  inv_err <- sapply(pcells, function(c) max(abs((c$phat_mean - b) / a - PS)))
  cat("\n")
  cat(sprintf("  monotone in every cell              : %s\n", if (mono_all) "yes" else "NO"))
  cat(sprintf("  slope across the %2d cells           : %.3f to %.3f (spread %.0f%% of the mean)\n",
              length(pcells), min(sl), max(sl), 100 * (max(sl) - min(sl)) / mean(sl)))
  cat(sprintf("  intercept across cells              : %.3f to %.3f\n", min(ic), max(ic)))
  cat(sprintf("  reader output at p = 0              : %.3f to %.3f\n", min(z0), max(z0)))
  cat(sprintf("  linear fit R2                       : %.3f to %.3f\n", min(r2), max(r2)))
  cat(sprintf("  one pooled calibration phat = %.3f p + %.3f\n", a, b))
  cat(sprintf("  worst error in p after that single calibration, by cell : %.3f to %.3f\n",
              min(inv_err), max(inv_err)))
  cat("\n  READ THIS AS: the paper claims p is identifiable up to a STABLE calibratable\n")
  cat("  attenuation. Monotone recovery in every cell supports the first half. The second\n")
  cat("  half needs the slope to be about the same across targets and grids, since a\n")
  cat("  calibration that depends on the unknown density is not available in practice. If\n")
  cat("  the slope spread is large, or a single pooled calibration leaves errors in p of\n")
  cat("  more than a few hundredths, the sentence should be narrowed to monotonicity.\n")

  res$partB <- list(cells = pcells, p = PS, Ds = DS_B, n = n, nseed = NSEEDB,
                    pooled = c(slope = a, intercept = b),
                    summary = list(monotone_all = mono_all, slope_min = min(sl), slope_max = max(sl),
                                   intercept_min = min(ic), intercept_max = max(ic),
                                   at_zero_min = min(z0), at_zero_max = max(z0),
                                   r2_min = min(r2), r2_max = max(r2),
                                   inv_err_min = min(inv_err), inv_err_max = max(inv_err)),
                    paper_claim = c(p0.2 = 0.175, p0.4 = 0.334, p0.6 = 0.506))
}

## ---- save ------------------------------------------------------------------
saveRDS(list(provenance = list(
  script = "endtoend_v10_4.R",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = R.version.string, platform = R.version$platform,
  pkg_dir = normalizePath(pkg), n = n, M = M, dx = dx,
  nseed_A = NSEED, nseed_B = NSEEDB, part = PART, smoke = SMOKE,
  Ds_A = DS_A, Ds_B = DS_B, ps = PS,
  seed_rule_A = "20260627 + 1000*replicate + round(100*D)",
  seed_rule_B = "20260627 + 1000*replicate + round(100*D) + 10000*round(10*p)",
  elapsed_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs"))),
  results = res),
  file.path(OUT, "endtoend_v10_4.rds"))

cat(sprintf("\n-> %s\nelapsed %.1f s\n", file.path(OUT, "endtoend_v10_4.rds"),
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
