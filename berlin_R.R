#!/usr/bin/env Rscript
## berlin_R.R -- Berlin resident-register age controlled-grouping study in R, all methods,
## and the age-density figure (with an SEM curve). Usage: Rscript berlin_R.R <dir with CSV>
args <- commandArgs(trailingOnly = FALSE)
fa <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(fa)) dirname(normalizePath(fa)) else "."
for (f in list.files(file.path(here, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
ta <- commandArgs(trailingOnly = TRUE); datadir <- if (length(ta)) ta[1] else "../berlin-data"
has_km <- requireNamespace("Kernelheaping", quietly = TRUE)

csv <- read.csv2(file.path(datadir, "EWR201512E_Matrix.csv"), stringsAsFactors = FALSE, check.names = FALSE)
bandcols <- grep("^E_E[0-9][0-9]_[0-9]+$", names(csv), value = TRUE)
edges <- t(sapply(bandcols, function(nm) as.integer(strsplit(sub("E_E", "", nm), "_")[[1]])))
counts <- sapply(bandcols, function(nm) sum(as.numeric(gsub(",", ".", gsub("\\.", "", csv[[nm]]))), na.rm = TRUE))
lo <- edges[,1]; hi <- edges[,2]
set.seed(20260627)
N <- 40000; wts <- counts / sum(counts)
k <- sample(seq_along(counts), N, TRUE, wts)
x <- lo[k] + runif(N) * (hi[k] - lo[k])                      # dequantized native-grain age sample

M <- 2048; grid <- seq(-25, 135, length.out = M + 1)[-(M + 1)]; dx <- grid[2] - grid[1]; w <- .fftw(M, dx)
href <- 1.06 * sd(x) * N^(-1/5)
fref <- .renorm(.reconstruct(fft(bin_prob(x, grid)$p), exp(-0.5 * (href * w)^2), dx), dx)
ise <- function(fh) sum((fh - fref)^2) * dx

cat(sprintf("Berlin register: %d residents, %d bands; ISE x1e3 vs native reference (Kernelheaping=%s)\n",
            round(sum(counts)), length(counts), has_km))
cat(sprintf("%3s %8s %8s %8s %8s %8s %8s %8s\n","D","naive","deconv","SEM","Heitjan","de-heap","super","COMBINED"))
rows <- list()
for (D in c(5,10,15,20)) {
  y <- D * round(x / D)
  r <- list(naive=ise(naive_kde(y,grid)), deconv=ise(deconv_kde(y,D,grid)),
            heitjan=ise(heitjan_mi(y,D,grid,M=8)), deheap=ise(deheap_kde(y,D,grid)),
            super=ise(superpose_kde(y,D,grid)), comb=ise(as.numeric(adkde(y,D,grid))))
  fs <- if (has_km) sem_kde(y,D,grid) else NULL
  r$sem <- if (!is.null(fs)) ise(fs) else NA; r$Dhat <- heap_grid(y, grid, near = D)
  rows[[as.character(D)]] <- r
  cat(sprintf("%3d %8.3f %8.3f %8s %8.3f %8.3f %8.3f %8.3f\n", D, r$naive*1e3, r$deconv*1e3,
      if (is.na(r$sem)) "  -  " else sprintf("%.3f", r$sem*1e3), r$heitjan*1e3, r$deheap*1e3,
      r$super*1e3, r$comb*1e3))
}
saveRDS(rows, file.path(here, "berlin_results.rds"))

## figure: age density at 20-year grain, legend outside (horizontal), with an SEM curve
Dp <- 20; y <- Dp * round(x / Dp)
f_naive <- naive_kde(y, grid); f_comb <- as.numeric(adkde(y, Dp, grid))
f_sem <- if (has_km) sem_kde(y, Dp, grid) else NULL
m <- grid >= 0 & grid <= 100
pdf(file.path(here, "fig_berlin_age.pdf"), width = 4.6, height = 3.2)
par(mar = c(4, 4, 3, 1), xpd = FALSE)
plot(grid[m], fref[m], type = "l", lwd = 2.2, col = "black", xlab = "age (years)", ylab = "density",
     ylim = c(0, max(fref[m], f_naive[m]) * 1.05))
lines(grid[m], f_naive[m], lwd = 1.6, lty = 3, col = "grey45")
if (!is.null(f_sem)) lines(grid[m], f_sem[m], lwd = 1.6, col = "#1b7837")
lines(grid[m], f_comb[m], lwd = 1.8, col = "#b2182b")
par(xpd = NA)
leg <- c("native-grain reference", "naive KDE (20 yr)", if (!is.null(f_sem)) "SEM (Kernelheaping)", "combined de-heaping")
cols <- c("black", "grey45", if (!is.null(f_sem)) "#1b7837", "#b2182b")
ltys <- c(1, 3, if (!is.null(f_sem)) 1, 1)
legend(x = 50, y = max(fref[m], f_naive[m]) * 1.26, legend = leg, col = cols, lty = ltys, lwd = 1.8,
       horiz = TRUE, bty = "n", cex = 0.60, xjust = 0.5, seg.len = 1.4)
dev.off()
cat("-> berlin_results.rds ; fig_berlin_age.pdf\n")
