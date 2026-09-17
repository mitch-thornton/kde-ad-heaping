#!/usr/bin/env Rscript
## test_lattice_v3.R
## Maps the regime of the mixed-grain reader and applies the third reader to the real
## cigarette counts. Four parts.
##
##   1  The theory check. A known-truth mixture at the cigarette sample size, with the
##      base scaled from the cigarette scale (mean 12) up by factors of two. The reader is
##      exact when the base is smooth at the coarsest grain and biased when it is not;
##      kappa = sd(base) / coarsest grain is the regime statistic.
##   2  The cigarette-scale known-truth test with the third reader, the same test the
##      second reader failed, for a like-for-like comparison.
##   3  The third reader on the real 1019 counts, with bootstrap standard errors, and
##      kappa estimated for the real data.
##   4  The model-based refinement (NOT adopted), on known truth and on the real counts,
##      to show that a base model recovers most of what the comb alone cannot.
##
## Usage (bash, not zsh):
##     NHANES_DIR=/Users/mitch/src/KDE-AD-HEAPING/data/NHANES bash run_test_lattice_v3.sh
## Part 3 and the real-data half of part 4 need NHANES_DIR; the rest runs without it.
## Environment: PKG_DIR (only for results dir default), NHANES_DIR, LATTICE_V3, NSIM, NBOOT, OUTDIR.

t0 <- Sys.time()
pkg <- Sys.getenv("PKG_DIR", unset = "/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping")
v3  <- Sys.getenv("LATTICE_V3", unset = file.path(dirname(pkg), "heap_lattice_v3.R"))
if (!file.exists(v3)) stop("heap_lattice_v3.R not found at ", v3, "; set LATTICE_V3")
source(v3)
NH    <- Sys.getenv("NHANES_DIR", unset = "")
NSIM  <- as.integer(Sys.getenv("NSIM", "200"))
NBOOT <- as.integer(Sys.getenv("NBOOT", "1000"))
OUT   <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
GR <- c(1, 5, 10, 20); TRUEW <- c(0.41, 0.18, 0.12, 0.16); TRUEU <- 0.13
TRUE4 <- c(unit = TRUEW[1] + TRUEU, g5 = TRUEW[2], g10 = TRUEW[3], g20 = TRUEW[4])
wfull <- c(TRUEW, TRUEU) / sum(c(TRUEW, TRUEU))
n0 <- 1019
hr <- function(t) cat("\n", strrep("=", 78), "\n ", t, "\n", strrep("=", 78), "\n", sep = "")

simulate <- function(mu, sd, n = n0) {
  base <- pmax(0.5, rgamma(n, shape = mu^2 / sd^2, scale = sd^2 / mu))
  comp <- sample(5, n, TRUE, wfull); g <- c(1, 5, 10, 20, 1)[comp]
  list(y = pmax(g, g * round(base / g)), base = base)
}
w4 <- function(r) { w <- r$weights; c(unit = unname(w["1"]), g5 = unname(w["5"]), g10 = unname(w["10"]), g20 = unname(w["20"])) }
res <- list()

## ---- 1 regime sweep --------------------------------------------------------
hr("1  known-truth recovery against the base scale (third reader)")
cat(sprintf("  n = %d, %d simulations per scale, base gamma with sd = 0.75 mean, grains 1 5 10 20\n", n0, NSIM))
cat(sprintf("  %-9s %6s | %7s %7s %7s %7s | %s\n", "base mean", "kappa", "unit", "g5", "g10", "g20", "ratios to truth"))
sweep <- list()
for (mu in c(12, 24, 48, 96, 192)) {
  sd <- 0.75 * mu; W <- matrix(NA, NSIM, 4); kap <- numeric(NSIM)
  for (s in seq_len(NSIM)) { set.seed(20260627 + s); si <- simulate(mu, sd); W[s, ] <- w4(heap_lattice3(si$y, GR)); kap[s] <- sd(si$base) / max(GR) }
  m <- colMeans(W); sweep[[as.character(mu)]] <- list(mu = mu, sd = sd, kappa = mean(kap), mean = m, se = apply(W, 2, sd) / sqrt(NSIM), W = W)
  cat(sprintf("  %-9d %6.2f | %7.3f %7.3f %7.3f %7.3f | %s\n", mu, mean(kap), m[1], m[2], m[3], m[4], paste(sprintf("%.2f", m / TRUE4), collapse = " ")))
}
cat(sprintf("  %-9s %6s | %7.3f %7.3f %7.3f %7.3f |\n", "truth", "", TRUE4[1], TRUE4[2], TRUE4[3], TRUE4[4]))
cat("  READ THIS AS: the reader is unbiased once the base is smooth at the coarsest grain\n",
    " (kappa above about 3) and biased below, worst on the coarsest grain. The cigarette\n",
    " data sit at the first row.\n")
res$sweep <- sweep

## ---- 2 cigarette-scale test, like for like ---------------------------------
hr("2  cigarette-scale known truth, third reader (the test the second reader failed)")
W <- matrix(NA, NSIM, 4); sh <- matrix(NA, NSIM, 3)
for (s in seq_len(NSIM)) { set.seed(20260627 + s); si <- simulate(12, 9); W[s, ] <- w4(heap_lattice3(si$y, GR)); sh[s, ] <- c(mean(si$y %% 5 == 0), mean(si$y %% 10 == 0), mean(si$y %% 20 == 0)) }
m <- colMeans(W)
cat(sprintf("  simulated shares 5/10/20 : %.3f %.3f %.3f\n", colMeans(sh)[1], colMeans(sh)[2], colMeans(sh)[3]))
cat(sprintf("  %-8s %7s %9s %7s\n", "quantity", "true", "recovered", "ratio"))
for (i in 1:4) cat(sprintf("  %-8s %7.3f %9.3f %7.2f\n", names(TRUE4)[i], TRUE4[i], m[i], m[i] / TRUE4[i]))
cat("  Second reader on the same construction gave ratios 1.15 (grain 1) 0.11 1.42 0.95;\n",
    " the unit share here is grain one plus the residual, 0.54.\n")
res$cig_scale <- list(mean = m, W = W, shares = colMeans(sh))

## ---- 3 real counts ---------------------------------------------------------
cig <- NULL
if (nzchar(NH) && file.exists(file.path(NH, "SMQ_J.xpt"))) {
  suppressWarnings(suppressMessages(library(foreign)))
  smq <- read.xport(file.path(NH, "SMQ_J.xpt"))
  cn <- intersect(c("SMD650", "SMD641", "SMQ020"), names(smq))
  if (length(cn)) { v <- suppressWarnings(as.numeric(smq[[cn[1]]])); cig <- v[is.finite(v) & v > 0 & v < 200] }
}
if (!is.null(cig) && length(cig) >= 200) {
  hr("3  third reader on the real NHANES counts")
  r <- heap_lattice3(cig, GR, nboot = NBOOT)
  w <- w4(r); se <- c(r$boot_se["1"], r$boot_se["5"], r$boot_se["10"], r$boot_se["20"])
  cat(sprintf("  n = %d\n", length(cig)))
  cat(sprintf("  exclusive-center amplitudes A20 %.3f  A10 %.3f  A5 %.3f\n", r$amplitudes["20"], r$amplitudes["10"], r$amplitudes["5"]))
  cat(sprintf("  %-8s %9s %8s %10s %10s\n", "quantity", "estimate", "boot SE", "previous", "repaired"))
  prev <- c(0.41 + 0.13, 0.18, 0.12, 0.16); rep2 <- c(0.442 + 0.117, 0.081, 0.190, 0.169)
  for (i in 1:4) cat(sprintf("  %-8s %9.3f %8.3f %10.3f %10.3f\n", names(TRUE4)[i], w[i], se[i], prev[i], rep2[i]))
  g0 <- ifelse(cig %% 20 == 0, 20, ifelse(cig %% 10 == 0, 10, ifelse(cig %% 5 == 0, 5, 1)))
  dj <- pmax(0.5, cig + runif(length(cig), -g0 / 2, g0 / 2))
  kap <- sd(dj) / max(GR)
  cat(sprintf("  kappa for the real data (sd of dequantized counts %.2f over grain 20): %.2f\n", sd(dj), kap))
  cat("  READ THIS AS: at this kappa part 1 says the reader is biased, so these weights are\n",
      " not a measurement either. The unit share and the heaped fraction 1 - unit are the\n",
      " quantities part 1 shows are recovered at every scale.\n")
  res$real <- list(n = length(cig), weights = w, se = se, amplitudes = r$amplitudes, kappa = kap, boot_w = r$boot_w)
} else cat("\n  (part 3 skipped: set NHANES_DIR to the directory holding SMQ_J.xpt)\n")

## ---- 4 model-based refinement, not adopted ---------------------------------
hr("4  model-based refinement (NOT adopted): known truth at the cigarette scale")
NREF <- min(20L, NSIM); W <- matrix(NA, NREF, 4)
for (s in seq_len(NREF)) { set.seed(20260627 + s); si <- simulate(12, 9); W[s, ] <- heap_lattice_refine(si$y, GR, iters = 6, seed = 20260627 + s)$weights }
m <- colMeans(W)
cat(sprintf("  %d simulations, 6 iterations each\n", NREF))
cat(sprintf("  %-8s %7s %9s %7s %7s\n", "quantity", "true", "recovered", "ratio", "sd"))
for (i in 1:4) cat(sprintf("  %-8s %7.3f %9.3f %7.2f %7.3f\n", names(TRUE4)[i], TRUE4[i], m[i], m[i] / TRUE4[i], sd(W[, i])))
cat("  READ THIS AS: with a base estimated from the data by posterior dequantization, most\n",
    " of the mixture is recovered at the cigarette scale, at the price of a rounding model\n",
    " and an imputation loop, which is the lineage the paper contrasts itself with.\n")
res$refine_sim <- list(mean = m, W = W)
if (!is.null(cig) && length(cig) >= 200) {
  rr <- heap_lattice_refine(cig, GR, iters = 8)
  cat("\n  refinement on the real counts, weights by iteration (unit, g5, g10, g20):\n")
  for (it in seq_len(nrow(rr$trajectory))) cat(sprintf("    iter %d  %s\n", it, paste(sprintf("%.3f", rr$trajectory[it, ]), collapse = "  ")))
  res$refine_real <- rr
}

saveRDS(list(provenance = list(script = "test_lattice_v3.R", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = R.version.string, platform = R.version$platform, lattice_v3 = normalizePath(v3),
  nsim = NSIM, nboot = NBOOT, grains = GR, truth = TRUE4, seed_rule = "20260627 + index"), results = res),
  file.path(OUT, "test_lattice_v3.rds"))
cat(sprintf("\n-> %s\nelapsed %.1f s\n", file.path(OUT, "test_lattice_v3.rds"), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
