#!/usr/bin/env Rscript
## verify_mw15.R -- self-checks on the Marron-Wand definitions, and two things the
## v10 revision needs to know before the supplementary benchmark is run:
##   (a) which of the four densities currently in Table 1 are actually Marron-Wand members
##   (b) how fine the evaluation grid has to be for each density before the grid itself,
##       rather than the estimator, dominates the integrated squared error
source(file.path(dirname(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])), "mw15.R"))

cat("=== 1. weights sum to one, and the density integrates to one ===\n")
for (d in MW15) {
  xg <- seq(-15, 15, length.out = 2^20); dxg <- xg[2] - xg[1]
  I <- sum(mw_d(xg, d)) * dxg
  cat(sprintf("%2d %-24s ncomp=%2d  sum(w)=%.12f  integral=%.8f  min(sd)=%.4f\n",
              d$index, d$name, length(d$w), sum(d$w), I, mw_min_sd(d)))
}

cat("\n=== 2. mean and variance against numerical quadrature ===\n")
for (d in MW15) {
  xg <- seq(-15, 15, length.out = 2^20); dxg <- xg[2] - xg[1]; f <- mw_d(xg, d)
  m_num <- sum(xg * f) * dxg
  v_num <- sum((xg - m_num)^2 * f) * dxg
  m_an <- sum(d$w * d$mu)
  v_an <- sum(d$w * (d$sd^2 + d$mu^2)) - m_an^2
  cat(sprintf("%2d %-24s mean %9.6f vs %9.6f   var %9.6f vs %9.6f   %s\n",
              d$index, d$name, m_num, m_an, v_num, v_an,
              if (abs(m_num - m_an) < 1e-6 && abs(v_num - v_an) < 1e-6) "ok" else "MISMATCH"))
}

cat("\n=== 3. do the four densities in Table 1 match Marron-Wand members? ===\n")
PAPER <- list(
  gaussian = list(w = 1, mu = 0, sd = 1),
  bimodal  = list(w = c(.5, .5), mu = c(-1.2, 1.2), sd = c(.5, .5)),
  kurtotic = list(w = c(2/3, 1/3), mu = c(0, 0), sd = c(1, .1)),
  skewed   = list(w = rep(.2, 5), mu = c(0, .5, 1.0833, 1.4167, 1.6875),
                  sd = c(1, 2/3, 4/9, 8/27, 16/81)))
xg <- seq(-8, 8, length.out = 2^18); dxg <- xg[2] - xg[1]
for (pn in names(PAPER)) {
  fp <- mw_d(xg, PAPER[[pn]])
  best <- NULL
  for (d in MW15) {
    l1 <- sum(abs(fp - mw_d(xg, d))) * dxg
    if (is.null(best) || l1 < best$l1) best <- list(l1 = l1, idx = d$index, nm = d$name)
  }
  exact <- best$l1 < 1e-9
  cat(sprintf("paper '%-8s' : closest MW member is #%2d %-24s  L1 distance %8.5f  %s\n",
              pn, best$idx, best$nm, best$l1,
              if (exact) "IDENTICAL" else "NOT a Marron-Wand density"))
}

cat("\n=== 4. grid resolution floor per density ===\n")
cat("ISE x1e3 of the closed-form density against its own binned-and-reconstructed version.\n")
cat("This is the error no estimator on that grid can go below. Table 1 uses M=2048 on\n")
cat("[-10,10), so dx = 0.009766.\n\n")
floor_ise <- function(d, M, lo, hi) {
  grid <- seq(lo, hi, length.out = M + 1)[-(M + 1)]; dx <- grid[2] - grid[1]
  ft <- mw_d(grid, d)
  ## exact mass per cell, then the piecewise-constant density that a binning recovers
  edges <- c(grid - dx/2, grid[M] + dx/2)
  p <- numeric(M)
  for (j in seq_along(d$w)) p <- p + d$w[j] * diff(pnorm(edges, d$mu[j], d$sd[j]))
  fb <- p / dx
  sum((fb - ft)^2) * dx
}
cfgs <- list(c(2048, -10, 10), c(8192, -10, 10), c(8192, -5, 5), c(32768, -10, 10))
hdr <- sapply(cfgs, function(c) sprintf("M=%-5d[%g,%g]", c[1], c[2], c[3]))
cat(sprintf("%2s %-24s %8s %14s %14s %14s %14s\n", "#", "density", "min sd", hdr[1], hdr[2], hdr[3], hdr[4]))
for (d in MW15) {
  vals <- sapply(cfgs, function(c) floor_ise(d, c[1], c[2], c[3]) * 1e3)
  cat(sprintf("%2d %-24s %8.4f %14.4f %14.4f %14.4f %14.4f\n",
              d$index, d$name, mw_min_sd(d), vals[1], vals[2], vals[3], vals[4]))
}

cat("\n=== 5. sd per bin at each configuration ===\n")
cat(sprintf("%2s %-24s %10s %10s %10s %10s\n", "#", "density", hdr[1], hdr[2], hdr[3], hdr[4]))
for (d in MW15) {
  bins <- sapply(cfgs, function(c) mw_min_sd(d) / ((c[3] - c[2]) / c[1]))
  cat(sprintf("%2d %-24s %10.2f %10.2f %10.2f %10.2f\n", d$index, d$name,
              bins[1], bins[2], bins[3], bins[4]))
}
