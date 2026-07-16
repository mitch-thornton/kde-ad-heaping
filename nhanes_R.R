#!/usr/bin/env Rscript
## nhanes_R.R -- NHANES weight controlled-coarsening study in R with all methods.
## Usage: Rscript nhanes_R.R <data-dir with DEMO_J.xpt BMX_J.xpt>
args <- commandArgs(trailingOnly = FALSE)
fa <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(fa)) dirname(normalizePath(fa)) else "."
for (f in list.files(file.path(here, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
ta <- commandArgs(trailingOnly = TRUE); datadir <- if (length(ta)) ta[1] else "../data"
has_km <- requireNamespace("Kernelheaping", quietly = TRUE)
suppressWarnings(suppressMessages(library(foreign)))

KG2LB <- 2.2046226218
demo <- read.xport(file.path(datadir, "DEMO_J.xpt"))[, c("SEQN", "RIDAGEYR")]
bmx  <- read.xport(file.path(datadir, "BMX_J.xpt"))[, c("SEQN", "BMXWT")]
d <- merge(demo, bmx, by = "SEQN"); d <- d[d$RIDAGEYR >= 20 & !is.na(d$BMXWT), ]
wt <- d$BMXWT * KG2LB; n <- length(wt)

M <- 4096; grid <- seq(min(wt) - 40, max(wt) + 40, length.out = M + 1)[-(M + 1)]; dx <- grid[2] - grid[1]
w <- .fftw(M, dx)
href <- 1.06 * sd(wt) * n^(-1/5)
fref <- .renorm(.reconstruct(fft(bin_prob(wt, grid)$p), exp(-0.5 * (href * w)^2), dx), dx)
ise <- function(fh) sum((fh - fref)^2) * dx

cat(sprintf("NHANES weights n=%d ; ISE x1e3 vs no-heaping reference (Kernelheaping=%s)\n", n, has_km))
cat(sprintf("%3s %8s %8s %8s %8s %8s %8s %8s\n", "D","naive","deconv","SEM","Heitjan","de-heap","super","COMBINED"))
rows <- list()
for (D in c(5,10,20,30,40)) {
  y <- D * round(wt / D)
  r <- list(naive=ise(naive_kde(y,grid)), deconv=ise(deconv_kde(y,D,grid)),
            heitjan=ise(heitjan_mi(y,D,grid,M=8)), deheap=ise(deheap_kde(y,D,grid)),
            super=ise(superpose_kde(y,D,grid)), comb=ise(as.numeric(adkde(y,D,grid))))
  fs <- if (has_km) sem_kde(y,D,grid) else NULL
  r$sem <- if (!is.null(fs)) ise(fs) else NA
  r$Dhat <- heap_grid(y, grid, near = D)
  rows[[as.character(D)]] <- r
  cat(sprintf("%3d %8.3f %8.3f %8s %8.3f %8.3f %8.3f %8.3f\n", D, r$naive*1e3, r$deconv*1e3,
      if (is.na(r$sem)) "  -  " else sprintf("%.3f", r$sem*1e3), r$heitjan*1e3, r$deheap*1e3,
      r$super*1e3, r$comb*1e3))
}
saveRDS(rows, file.path(here, "nhanes_results.rds"))

## colored two-panel figure consistent with the paper's methods
Ds <- c(5,10,20,30,40)
curves <- sapply(c("naive","deconv","sem","comb"), function(k) sapply(as.character(Ds), function(d) rows[[d]][[k]]))
pdf(file.path(here, "fig_nhanes_weight.pdf"), width = 7.2, height = 3.0)
par(mfrow = c(1,2), mar = c(4,4,2,1))
cols <- c(naive="grey45", deconv="#2166ac", sem="#1b7837", comb="#b2182b")
lts  <- c(naive=3, deconv=2, sem=1, comb=1)
matplot(Ds, log10(pmax(curves,1e-6)*1e3), type="n", xlab="imposed grid D (lb)",
        ylab=expression(log[10]~ISE~x10^3), main="(A) robustness to coarsening")
for (k in names(cols)) lines(Ds, log10(pmax(curves[,k],1e-6)*1e3), col=cols[k], lty=lts[k], lwd=2)
legend("topleft", c("naive KDE","deconvolution","SEM (Kernelheaping)","combined de-heaping"),
       col=cols, lty=lts, lwd=2, bty="o", box.lwd=0.6, bg="white", cex=0.6)
D <- 30; y <- D*round(wt/D)
plot(grid, fref, type="l", lwd=2.2, xlim=c(90,320), xlab="weight (lb)", ylab="density",
     main="(B) weight rounded to 30 lb")
lines(grid, naive_kde(y,grid), lty=3, lwd=1.6, col="grey45")
lines(grid, as.numeric(adkde(y,D,grid)), lwd=1.8, col="#b2182b")
legend("topright", c("no-heaping reference","naive KDE","combined de-heaping"),
       col=c("black","grey45","#b2182b"), lty=c(1,3,1), lwd=1.8, bty="o", box.lwd=0.6, bg="white", cex=0.6)
dev.off()
cat("-> nhanes_results.rds ; fig_nhanes_weight.pdf\n")
