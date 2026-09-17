#!/usr/bin/env Rscript
## b6_fix_v10_3.R
## Corrected replacement for the B6 block of batch_v10_3.R.
##
## The first B6 was mine and it was wrong. It resampled the base counts from the real
## NHANES reports and then applied a known grain mixture on top. Those reports are
## already heaped, so the simulation rounded twice. The symptom is in the output you
## returned: the simulated share of multiples of five came out at 0.768 against 0.570
## in the real data, ten at 0.598 against 0.430, and twenty at 0.355 against 0.220,
## inflated by factors of 1.35 to 1.61. A simulation that is more heaped than the data
## it is meant to mimic cannot test how the reader behaves on that data.
##
## The fix is to build a base that is NOT already heaped, then impose the known mixture.
## Three candidate bases are offered and all three are run, because which one is the
## fairest stand-in is a judgment call and it is better to see the spread.
##
##   dejitter  real counts, uniformly dequantized within their apparent grain, then
##             resampled. Keeps the shape of the real distribution, removes the comb.
##   nonround  only the real counts that are NOT multiples of five. These are the
##             reports least likely to have been rounded, so they are the closest thing
##             to an unheaped sample the data contains. Resampled with jitter.
##   gamma     a smooth parametric stand-in fitted to the mean and variance of the
##             dejittered counts. Fully smooth by construction.
##
## The block also reports, for each base, whether the simulated round-number shares
## match the real ones. That check is the gate. If the simulated shares do not land
## near 0.570 / 0.430 / 0.220, the simulation is not a stand-in for the real data and
## its recovery numbers should not be used.
##
## Usage (bash, not zsh):
##     NHANES_DIR=/Users/mitch/src/KDE-AD-HEAPING/data/NHANES bash run_b6_fix_v10_3.sh
##
## Environment: PKG_DIR, NHANES_DIR, NSIM (default 200), OUTDIR.

t0 <- Sys.time()
pkg <- Sys.getenv("PKG_DIR", unset = "/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping")
if (!dir.exists(file.path(pkg, "R"))) stop("no R/ under PKG_DIR=", pkg)
for (f in list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
NH   <- Sys.getenv("NHANES_DIR", unset = "")
NSIM <- as.integer(Sys.getenv("NSIM", "200"))
OUT  <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GRAINS <- c(1, 5, 10, 20)
TRUEW  <- c(0.41, 0.18, 0.12, 0.16); TRUEU <- 0.13
wfull  <- c(TRUEW, TRUEU) / sum(c(TRUEW, TRUEU))
REAL   <- c(m5 = 0.570, m10 = 0.430, m20 = 0.220)   # the paper's observed shares

## ---- real counts -----------------------------------------------------------
cig <- NULL
if (nzchar(NH) && file.exists(file.path(NH, "SMQ_J.xpt"))) {
  suppressWarnings(suppressMessages(library(foreign)))
  smq <- read.xport(file.path(NH, "SMQ_J.xpt"))
  cn <- intersect(c("SMD650","SMD641","SMQ020"), names(smq))
  if (length(cn)) {
    v <- suppressWarnings(as.numeric(smq[[cn[1]]]))
    cig <- v[is.finite(v) & v > 0 & v < 200]
  }
}
if (is.null(cig) || length(cig) < 200)
  stop("B6 needs the real counts. Set NHANES_DIR to the directory holding SMQ_J.xpt.")
nn <- length(cig)
cat(sprintf("b6_fix_v10_3   n = %d real reports, %d simulations per base\n\n", nn, NSIM))
cat(sprintf("real data shares      5/10/20 : %.3f %.3f %.3f\n",
            mean(cig %% 5 == 0), mean(cig %% 10 == 0), mean(cig %% 20 == 0)))
cat(sprintf("paper's stated shares 5/10/20 : %.3f %.3f %.3f\n\n", REAL[1], REAL[2], REAL[3]))

## ---- three candidate unheaped bases ----------------------------------------
apparent_grain <- function(x) ifelse(x %% 20 == 0, 20, ifelse(x %% 10 == 0, 10,
                              ifelse(x %% 5 == 0, 5, 1)))
mk_base <- function(kind, m) {
  if (kind == "dejitter") {
    g <- apparent_grain(cig)
    pmax(0.5, cig + runif(nn, -g/2, g/2))            # dequantize within apparent grain
  } else if (kind == "nonround") {
    pool <- cig[cig %% 5 != 0]
    if (length(pool) < 50) return(NULL)
    pmax(0.5, sample(pool, m, TRUE) + runif(m, -0.5, 0.5))
  } else {
    g <- apparent_grain(cig); dj <- pmax(0.5, cig + runif(nn, -g/2, g/2))
    mu <- mean(dj); v <- stats::var(dj)
    pmax(0.5, rgamma(m, shape = mu^2/v, scale = v/mu))
  }
}

BASES <- c("dejitter", "nonround", "gamma")
res <- list()
for (kind in BASES) {
  rec <- matrix(NA_real_, NSIM, length(GRAINS)); unr <- numeric(NSIM)
  sh <- matrix(NA_real_, NSIM, 3)
  ok <- TRUE
  for (s in seq_len(NSIM)) {
    set.seed(20260627 + s)
    base <- mk_base(kind, nn)
    if (is.null(base)) { ok <- FALSE; break }
    comp <- sample(seq_along(wfull), nn, TRUE, wfull)
    y <- pmax(1, round(base))                        # grain one, and the residual class
    for (j in seq_along(GRAINS)) {
      g <- GRAINS[j]; sel <- comp == j
      if (g > 1) y[sel] <- pmax(g, g * round(base[sel] / g))
    }
    r <- tryCatch(heap_lattice(y, grains = GRAINS), error = function(e) NULL)
    if (!is.null(r)) { rec[s, ] <- r$weights; unr[s] <- r$unrounded }
    sh[s, ] <- c(mean(y %% 5 == 0), mean(y %% 10 == 0), mean(y %% 20 == 0))
  }
  if (!ok) { cat(sprintf("base %-9s unavailable\n", kind)); next }
  ms <- colMeans(sh, na.rm = TRUE)
  gate <- all(abs(ms - REAL) < 0.06)
  cat(sprintf("== base: %s ==\n", kind))
  cat(sprintf("  simulated shares 5/10/20 : %.3f %.3f %.3f   against real %.3f %.3f %.3f   %s\n",
              ms[1], ms[2], ms[3], REAL[1], REAL[2], REAL[3],
              if (gate) "MATCHES, usable" else "does NOT match, do not use"))
  cat(sprintf("  %-10s %8s %10s %8s\n", "quantity", "true", "recovered", "ratio"))
  for (j in seq_along(GRAINS))
    cat(sprintf("  grain %-4d %8.3f %10.3f %8.3f\n", GRAINS[j], TRUEW[j],
                mean(rec[,j], na.rm=TRUE), mean(rec[,j], na.rm=TRUE)/TRUEW[j]))
  cat(sprintf("  %-10s %8.3f %10.3f %8.3f\n", "residual", TRUEU, mean(unr), mean(unr)/TRUEU))
  rr <- mean(rec[,1], na.rm=TRUE)/TRUEW[1]
  dir <- sapply(seq_along(GRAINS), function(j) mean(rec[,j], na.rm=TRUE)/TRUEW[j])
  cat(sprintf("  ratios across grains: %s\n", paste(sprintf("%.2f", dir), collapse = " ")))
  cat(sprintf("  uniform attenuation would put every ratio below 1 and close together. %s\n\n",
              if (all(dir < 1)) "All are below 1." else "They are not all below 1."))
  res[[kind]] <- list(shares = ms, gate = gate, recovered = rec, unrounded = unr, ratios = dir)
}

cat("INTERPRETATION\n")
cat("Use only the bases whose simulated round-number shares match the real ones. On those,\n")
cat("if the recovery ratios are all below one and close together, the attenuation account\n")
cat("in the paper is supported. If they straddle one, it is not, and the sentence should be\n")
cat("withdrawn rather than defended. Note also that the shipped R reader defines its\n")
cat("residual share as one minus the sum of the non-unit weights, which with weights that\n")
cat("sum to one makes it identically the grain-one weight, so the residual column carries\n")
cat("no information beyond the first row.\n")

saveRDS(list(provenance = list(script = "b6_fix_v10_3.R",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = R.version.string, platform = R.version$platform,
  pkg_dir = normalizePath(pkg), n = nn, nsim = NSIM,
  seed_rule = "20260627 + simulation index"), results = res),
  file.path(OUT, "b6_fix_v10_3.rds"))
cat(sprintf("\n-> %s\nelapsed %.1f s\n", file.path(OUT, "b6_fix_v10_3.rds"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
