#!/usr/bin/env Rscript
## benchmark.R -- head-to-head of the combined AD heaping-KDE against three mainstream
## external methods (deconvolution KDE; real Kernelheaping SEM; Heitjan-Rubin MI replica)
## on the four Gaussian-mixture targets, ISE against the closed-form density. Seed 20260627.
## locate the package R/ directory relative to this script (or PKG_DIR env)
args <- commandArgs(trailingOnly = FALSE)
fa <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(fa)) dirname(normalizePath(fa)) else Sys.getenv("PKG_DIR", ".")
for (f in list.files(file.path(here, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
has_km <- requireNamespace("Kernelheaping", quietly = TRUE)

## ---- targets: Gaussian mixtures with closed-form density ----
TARGETS <- list(
  gaussian = list(w = 1, mu = 0, sd = 1),
  bimodal  = list(w = c(.5, .5), mu = c(-1.2, 1.2), sd = c(.5, .5)),
  kurtotic = list(w = c(2/3, 1/3), mu = c(0, 0), sd = c(1, .1)),
  skewed   = list(w = rep(.2, 5), mu = c(0, .5, 1.0833, 1.4167, 1.6875),
                  sd = c(1, 2/3, 4/9, 8/27, 16/81)))
dmix <- function(x, t) rowSums(sapply(seq_along(t$w), function(j) t$w[j] * dnorm(x, t$mu[j], t$sd[j])))
rmix <- function(n, t) { k <- sample(seq_along(t$w), n, TRUE, t$w); rnorm(n, t$mu[k], t$sd[k]) }

ise <- function(fh, ft, dx) sum((fh - ft)^2) * dx
M <- 2048; grid <- seq(-10, 10, length.out = M + 1)[- (M + 1)]; dx <- grid[2] - grid[1]
n <- 4000; Ds <- c(0.25, 0.5, 1.0, 1.5); nseed <- 8

cat(sprintf("adheaping benchmark  n=%d  seeds=%d  Kernelheaping=%s\n", n, nseed, has_km))
cat(sprintf("%-9s %5s | %7s %7s %7s %7s | %7s %7s | %8s\n",
            "target","D","naive","deconv","SEM","Heitjan","de-heap","super","COMBINED"))
rows <- list()
for (nm in names(TARGETS)) {
  t <- TARGETS[[nm]]; ft <- dmix(grid, t)
  for (D in Ds) {
    e <- list(naive=c(), deconv=c(), sem=c(), heitjan=c(), deheap=c(), super=c(), comb=c())
    for (s in seq_len(nseed)) {
      set.seed(20260627 + 1000*s + round(D*100))
      y <- D * round(rmix(n, t) / D)
      e$naive   <- c(e$naive,   ise(naive_kde(y, grid), ft, dx))
      e$deconv  <- c(e$deconv,  ise(deconv_kde(y, D, grid), ft, dx))
      e$heitjan <- c(e$heitjan, ise(heitjan_mi(y, D, grid, M = 8), ft, dx))
      e$deheap  <- c(e$deheap,  ise(deheap_kde(y, D, grid), ft, dx))
      e$super   <- c(e$super,   ise(superpose_kde(y, D, grid), ft, dx))
      e$comb    <- c(e$comb,    ise(as.numeric(adkde(y, D, grid)), ft, dx))
      if (has_km) { fs <- sem_kde(y, D, grid); if (!is.null(fs)) e$sem <- c(e$sem, ise(fs, ft, dx)) }
    }
    mn <- lapply(e, function(v) if (length(v)) mean(v) else NA)
    rows[[paste(nm, D)]] <- c(target = nm, D = D, lapply(mn, function(v) round(v * 1e3, 3)))
    cat(sprintf("%-9s %5.2f | %7.2f %7.2f %7s %7.2f | %7.2f %7.2f | %8.2f\n",
        nm, D, mn$naive*1e3, mn$deconv*1e3,
        if (is.na(mn$sem)) "   -  " else sprintf("%.2f", mn$sem*1e3),
        mn$heitjan*1e3, mn$deheap*1e3, mn$super*1e3, mn$comb*1e3))
  }
}
saveRDS(rows, file.path(here, "benchmark_results.rds"))
cat("\n-> benchmark_results.rds\n")
