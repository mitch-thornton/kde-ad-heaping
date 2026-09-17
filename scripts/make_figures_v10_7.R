#!/usr/bin/env Rscript
## make_figures_v10_7.R
## Regenerates the three figures that still carry pending markers, in R, from the same
## library that computes every table. Closes B15, B2, B3 and, as a by-product, B4.
##
##   PART 1  Figure 1, the identifiability band (B15). The panel (b) curve in the
##           submitted figure was normalized by the supremum of |phi| over the band, and
##           its title still read "exact in-band", which is the claim Reviewer 2 had
##           withdrawn. Both normalizations are now drawn and named, so the pointwise
##           error of order one at the band edge, which is what Proposition 1 bounds, is
##           visible rather than hidden by the normalization. No sampling, so this part
##           is exact and takes a second.
##
##   PART 2  Figure 2, integrated squared error against the heaping grid (B2, C8). The
##           submitted figure carried four methods at eight seeds and omitted the
##           combined estimator, which the caption nonetheless named. It is redrawn with
##           the comparison set of Table 1 at the same seed count, from this library.
##
##   PART 3  Figure 4 and Table 2, the NHANES controlled study (B3, B4). Panel (A) plotted
##           log10 of the error while Table 2 printed the error, so a reader could not move
##           between them; it now uses a logarithmic axis whose tick labels are the table's
##           own numbers. The run also stores every value at full precision, so Table 2 can
##           be emitted rather than typed, which is what B4 asks for and which removes the
##           last hand-typed table from the paper.
##
## Figure 3 (spectral detector) and Figure 5 (Berlin) carry no pending marker and are not
## touched. Legends are placed below the axes rather than inside the plot box.
##
## Usage (bash, not zsh):
##     PKG_DIR=/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping \
##     NHANES_DIR=/Users/mitch/src/KDE-AD-HEAPING/data/NHANES \
##     bash run_make_figures_v10_7.sh
##
## Environment:
##   PKG_DIR     the kde-ad-heaping checkout, the R library          (required)
##   NHANES_DIR  directory holding DEMO_J.xpt and BMX_J.xpt          (part 3 only)
##   MW15_R      path to mw15.R                                      (default: beside this script)
##   FIGDIR      where the PDFs go       (default: ../figures beside this script)
##   OUTDIR      where the rds files go  (default: $PKG_DIR/results)
##   NSEED       replicates for part 2                               (default 50)
##   NHSEED      seed anchoring the drawing baselines in part 3       (default 20260627)
##   PARTS       comma list from 1,2,3                               (default 1,2,3)
##   NO_SEM      set to 1 to skip the Kernelheaping baseline
##   SMOKE       set to 1 for a fast reduced run
##
## Seed rule for part 2, the rule of record: 20260627 + 1000*replicate + round(100*D).

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

FIG    <- Sys.getenv("FIGDIR", unset = normalizePath(file.path(HERE, "..", "figures"), mustWork = FALSE))
OUT    <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
NHDIR  <- Sys.getenv("NHANES_DIR", unset = "")
NSEED  <- as.integer(Sys.getenv("NSEED", "50"))
PARTS  <- strsplit(Sys.getenv("PARTS", "1,2,3"), ",")[[1]]
NO_SEM <- nzchar(Sys.getenv("NO_SEM", ""))
NHSEED <- as.integer(Sys.getenv("NHSEED", "20260627"))   # anchors the drawing baselines in part 3
SMOKE  <- nzchar(Sys.getenv("SMOKE", ""))
if (SMOKE) NSEED <- 3L
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
has_km <- (!NO_SEM) && requireNamespace("Kernelheaping", quietly = TRUE)

hr <- function(t) cat("\n", strrep("=", 78), "\n ", t, "\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("make_figures_v10_7  parts=%s  seeds=%d  SEM=%s\n  figures -> %s\n  results -> %s\n",
            paste(PARTS, collapse = ","), NSEED, has_km, FIG, OUT))

## the bimodal target of Table 1, defined exactly as benchmark_v10_2.R defines it
BIMODAL <- list(name = "Bimodal", w = c(.5, .5), mu = c(-1.2, 1.2), sd = c(.5, .5))
sincf <- function(u) ifelse(abs(u) < 1e-12, 1, sin(u) / u)

## Legends live in their own layout row, below the axes and above the caption, never
## inside the plot box. legend_cell() is called once the layout has advanced to that row.
legend_cell <- function(labels, cols, ltys, lwds = 2, ncol = 1, cex = 0.62) {
  op <- par(mar = c(0, 0, 0, 0)); plot.new()
  legend("top", legend = labels, col = cols, lty = ltys, lwd = lwds, ncol = ncol,
         bty = "n", cex = cex, seg.len = 2.4, x.intersp = 0.8, y.intersp = 1.1)
  par(op)
}

## ============================================================================
## PART 1. Figure 1, the identifiability band
## ============================================================================
if ("1" %in% PARTS) {
  hr("1  Figure 1, the identifiability band, and the normalization of panel (b)")
  D <- 0.7
  wt <- seq(0.01, 1.6 * pi / D, length.out = 600)
  ms <- -8:8
  phi_true <- mw_phi(wt, BIMODAL)
  phiY <- sapply(wt, function(x) sum(mw_phi(x + 2 * pi * ms / D, BIMODAL) *
                                     sincf((x + 2 * pi * ms / D) * D / 2)))
  phi_desh <- phiY / sincf(wt * D / 2)
  err_abs  <- Mod(phi_desh - phi_true)
  err_point <- err_abs / Mod(phi_true)          # what Proposition 1 bounds
  err_sup   <- err_abs / max(Mod(phi_true))     # the normalization the previous figure used
  frac <- wt / (pi / D)
  at <- sapply(c(0.6, 0.9, 1.0), function(f) which.min(abs(frac - f)))

  pdf(file.path(FIG, "fig_identifiability.pdf"), width = 7.2, height = 3.8)
  layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE), heights = c(4.0, 1.25))
  par(mar = c(4.0, 4.2, 2.4, 1.0), mgp = c(2.3, 0.8, 0), cex.main = 0.92)
  ## panel (a)
  plot(wt, Mod(phi_true), type = "l", lwd = 2.2, col = "black", ylim = c(0, 1.6),
       xlab = "frequency w", ylab = "modulus", main = "(a) de-Sheppard and the band")
  lines(wt, Mod(phiY), lty = 2, lwd = 1.6, col = "#b2182b")
  lines(wt, Mod(phi_desh), lty = 4, lwd = 1.6, col = "#2166ac")
  abline(v = pi / D, col = "grey40", lty = 3)
  text(pi / D * 1.04, 1.5, expression(pi/D), col = "grey30", cex = 0.8)
  ## panel (b), both normalizations
  ylo <- max(min(c(err_point, err_sup)), 1e-11)
  plot(frac, err_point, type = "l", log = "y", lwd = 1.8, col = "#b2182b",
       ylim = c(ylo, 3),
       xlab = "w / (pi/D)", ylab = "relative error", main = "(b) error after dividing by sinc")
  lines(frac, err_sup, lwd = 1.8, col = "#2166ac", lty = 2)
  abline(v = 1, col = "grey40", lty = 3); abline(h = 1, col = "grey60", lty = 2)
  ## legend row
  legend_cell(c("|phi| true", "|phi_Y| heaped", "|phi_Y / sinc| de-Sheppard"),
              c("black", "#b2182b", "#2166ac"), c(1, 2, 4), cex = 0.62)
  legend_cell(c("pointwise, divided by |phi(w)|", "divided by sup |phi| over the band"),
              c("#b2182b", "#2166ac"), c(1, 2), cex = 0.62)
  dev.off()

  cat(sprintf("  bimodal target, D = %.2f, replica sum truncated at |m| = %d\n", D, max(ms)))
  cat(sprintf("  %-14s %12s %12s\n", "band fraction", "pointwise", "sup-norm"))
  for (i in seq_along(at))
    cat(sprintf("  %-14.2f %12.3e %12.3e\n", frac[at[i]], err_point[at[i]], err_sup[at[i]]))
  cat("  READ THIS AS: the pointwise error reaches order one at the band edge, which is what\n")
  cat("  the text states and Proposition 1 bounds. The sup-norm curve is the one the\n")
  cat("  submitted figure drew, and it is about two orders of magnitude lower there. The\n")
  cat("  pointwise curve rises without bound at the zeros of |phi|, where a relative error\n")
  cat("  has no meaning; the panel is capped so the band-edge value of one stays legible.\n")
  saveRDS(list(provenance = list(script = "make_figures_v10_7.R", part = 1,
                 generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                 r_version = R.version.string, platform = R.version$platform,
                 D = D, m_trunc = max(ms), target = "bimodal"),
               band_fraction = frac[at], pointwise = err_point[at], supnorm = err_sup[at],
               edge_pointwise = err_point[length(err_point)], curves = list(
                 frac = frac, pointwise = err_point, supnorm = err_sup)),
          file.path(OUT, "fig1b_v10_7.rds"))
  cat(sprintf("  -> %s\n  -> %s\n", file.path(FIG, "fig_identifiability.pdf"),
              file.path(OUT, "fig1b_v10_7.rds")))
}

## ============================================================================
## PART 2. Figure 2, integrated squared error against the grid
## ============================================================================
if ("2" %in% PARTS) {
  hr("2  Figure 2, error against the heaping grid, with the Table 1 comparison set")
  M <- 2048; grid <- seq(-10, 10, length.out = M + 1)[-(M + 1)]; dx <- grid[2] - grid[1]
  n <- 4000
  Ds <- if (SMOKE) c(0.25, 1.0) else c(0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
  ft <- mw_d(grid, BIMODAL)
  ise <- function(fh) sum((fh - ft)^2) * dx
  METH <- c("naive", "deconv", "sem", "imput", "deheap", "super", "combined")
  curves <- matrix(NA_real_, length(Ds), length(METH), dimnames = list(as.character(Ds), METH))
  sds    <- curves
  cat(sprintf("  n = %d, %d seeds, M = %d, dx = %.6f, bimodal target\n", n, NSEED, M, dx))
  cat(sprintf("  %5s %9s %9s %9s %9s %9s\n", "D", "naive", "deconv", "SEM", "imput", "combined"))
  for (i in seq_along(Ds)) {
    D <- Ds[i]
    e <- setNames(lapply(METH, function(k) rep(NA_real_, NSEED)), METH)
    for (s in seq_len(NSEED)) {
      set.seed(20260627 + 1000 * s + round(D * 100))
      y <- D * round(mw_r(n, BIMODAL) / D)
      e$naive[s]    <- ise(naive_kde(y, grid))
      e$deconv[s]   <- ise(deconv_kde(y, D, grid))
      e$imput[s]    <- ise(heitjan_mi(y, D, grid, M = 8))
      e$deheap[s]   <- ise(deheap_kde(y, D, grid))
      e$super[s]    <- ise(superpose_kde(y, D, grid))
      e$combined[s] <- ise(as.numeric(adkde(y, D, grid)))
      if (has_km) { fs <- sem_kde(y, D, grid); if (!is.null(fs)) e$sem[s] <- ise(fs) }
    }
    curves[i, ] <- sapply(e, function(v) mean(v, na.rm = TRUE))
    sds[i, ]    <- sapply(e, function(v) stats::sd(v, na.rm = TRUE))
    cat(sprintf("  %5.2f %9.3f %9.3f %9s %9.3f %9.3f\n", D, 1e3 * curves[i, "naive"],
                1e3 * curves[i, "deconv"],
                if (is.na(curves[i, "sem"])) "  -  " else sprintf("%.3f", 1e3 * curves[i, "sem"]),
                1e3 * curves[i, "imput"], 1e3 * curves[i, "combined"]))
  }

  ## Table 1 compares five; the two components are drawn lighter so the figure carries the
  ## same comparison set without becoming unreadable.
  main5 <- c("naive", "deconv", "sem", "imput", "combined")
  comp2 <- c("deheap", "super")
  col5 <- c(naive = "grey45", deconv = "#2166ac", sem = "#1b7837", imput = "#762a83", combined = "#b2182b")
  lty5 <- c(naive = 3, deconv = 2, sem = 4, imput = 5, combined = 1)
  colc <- c(deheap = "#f4a582", super = "#92c5de"); ltyc <- c(deheap = 1, super = 1)
  pdf(file.path(FIG, "fig_ise_grid.pdf"), width = 4.3, height = 4.4)
  layout(matrix(c(1, 2), 2, 1), heights = c(4.0, 1.9))
  par(mar = c(4.0, 4.2, 2.4, 1.0), mgp = c(2.3, 0.8, 0), cex.main = 0.95)
  yv <- 1e3 * curves[, c(main5, comp2), drop = FALSE]
  matplot(Ds, yv, type = "n", log = "y", xlab = "heaping grid D",
          ylab = expression(ISE %*% 10^3 ~ (bimodal)), main = "Error against coarsening")
  for (k in comp2) lines(Ds, 1e3 * curves[, k], col = colc[k], lty = ltyc[k], lwd = 1.4)
  for (k in main5) { lines(Ds, 1e3 * curves[, k], col = col5[k], lty = lty5[k], lwd = 2)
                     points(Ds, 1e3 * curves[, k], col = col5[k], pch = 16, cex = 0.55) }
  labs <- c(naive = "naive KDE", deconv = "deconvolution", sem = "SEM (Kernelheaping)",
            imput = "imput (Heitjan-Rubin)", combined = "combined de-heaping",
            deheap = "de-heaping component", super = "superposition component")
  keep <- c(main5, comp2); keep <- keep[sapply(keep, function(k) !all(is.na(curves[, k])))]
  legend_cell(unname(labs[keep]), c(col5, colc)[keep], c(lty5, ltyc)[keep],
              lwds = ifelse(keep %in% comp2, 1.4, 2), ncol = 2, cex = 0.58)
  dev.off()
  saveRDS(list(provenance = list(script = "make_figures_v10_7.R", part = 2,
                 generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                 r_version = R.version.string, platform = R.version$platform,
                 n = n, nseed = NSEED, M = M, dx = dx, Ds = Ds, target = "bimodal",
                 kernelheaping = has_km,
                 seed_rule = "20260627 + 1000*replicate + round(100*D)"),
               mean_ise = curves, sd_ise = sds),
          file.path(OUT, "fig2_ise_v10_7.rds"))
  cat(sprintf("  -> %s\n  -> %s\n", file.path(FIG, "fig_ise_grid.pdf"),
              file.path(OUT, "fig2_ise_v10_7.rds")))
}

## ============================================================================
## PART 3. Figure 4 and the NHANES table at full precision
## ============================================================================
if ("3" %in% PARTS) {
  if (!nzchar(NHDIR) || !file.exists(file.path(NHDIR, "DEMO_J.xpt"))) {
    cat("\n  (part 3 skipped: set NHANES_DIR to the directory holding DEMO_J.xpt and BMX_J.xpt)\n")
  } else {
    hr("3  Figure 4 and Table 2, NHANES controlled coarsening at full precision")
    ## Two of the five columns draw from the random number generator: the SEM baseline
    ## through Kernelheaping and the imputation baseline through its M = 8 draws. Part 3
    ## previously ran on whatever state the session happened to hold, so Table 2 was
    ## reproducible only to the precision at which those draws stop mattering. Anchoring
    ## on the seed of record makes it reproducible by construction instead.
    set.seed(NHSEED)
    suppressWarnings(suppressMessages(library(foreign)))
    KG2LB <- 2.2046226218
    demo <- read.xport(file.path(NHDIR, "DEMO_J.xpt"))[, c("SEQN", "RIDAGEYR")]
    bmx  <- read.xport(file.path(NHDIR, "BMX_J.xpt"))[, c("SEQN", "BMXWT")]
    d <- merge(demo, bmx, by = "SEQN"); d <- d[d$RIDAGEYR >= 20 & !is.na(d$BMXWT), ]
    wt <- d$BMXWT * KG2LB; n <- length(wt)
    M <- 4096; grid <- seq(min(wt) - 40, max(wt) + 40, length.out = M + 1)[-(M + 1)]
    dx <- grid[2] - grid[1]; w <- .fftw(M, dx)
    href <- 1.06 * stats::sd(wt) * n^(-1/5)
    fref <- .renorm(.reconstruct(stats::fft(bin_prob(wt, grid)$p), exp(-0.5 * (href * w)^2), dx), dx)
    ise <- function(fh) sum((fh - fref)^2) * dx
    Ds <- c(5, 10, 20, 30, 40)
    cat(sprintf("  n = %d, M = %d, dx = %.4f lb, reference bandwidth %.4f lb\n", n, M, dx, href))
    cat(sprintf("  %3s %10s %10s %10s %10s %10s %8s\n", "D", "naive", "deconv", "SEM", "imput", "combined", "Dhat"))
    rows <- list()
    for (D in Ds) {
      y <- D * round(wt / D)
      r <- list(D = D,
                naive = ise(naive_kde(y, grid)), deconv = ise(deconv_kde(y, D, grid)),
                imput = ise(heitjan_mi(y, D, grid, M = 8)), deheap = ise(deheap_kde(y, D, grid)),
                super = ise(superpose_kde(y, D, grid)), combined = ise(as.numeric(adkde(y, D, grid))))
      fs <- if (has_km) sem_kde(y, D, grid) else NULL
      r$sem <- if (!is.null(fs)) ise(fs) else NA_real_
      r$Dhat <- heap_grid(y, grid, near = D)
      rows[[as.character(D)]] <- r
      cat(sprintf("  %3d %10.4f %10.4f %10s %10.4f %10.4f %8.2f\n", D, 1e3*r$naive, 1e3*r$deconv,
                  if (is.na(r$sem)) "   -   " else sprintf("%.4f", 1e3*r$sem),
                  1e3*r$imput, 1e3*r$combined, r$Dhat))
    }
    saveRDS(list(provenance = list(script = "make_figures_v10_7.R", part = 3,
                   generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                   r_version = R.version.string, platform = R.version$platform,
                   n = n, M = M, dx = dx, href = href, Ds = Ds, kernelheaping = has_km,
                   seed = NHSEED,
                   dhat_mode = "verifying, imposed grid supplied as the hint"),
                 rows = rows), file.path(OUT, "nhanes_v10_7.rds"))

    cur <- sapply(c("naive", "deconv", "sem", "imput", "combined"),
                  function(k) sapply(as.character(Ds), function(dd) rows[[dd]][[k]]))
    pdf(file.path(FIG, "fig_nhanes_weight.pdf"), width = 7.2, height = 3.9)
    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE), heights = c(4.0, 1.35))
    par(mar = c(4.0, 4.2, 2.4, 1.0), mgp = c(2.3, 0.8, 0), cex.main = 0.95)
    cols <- c(naive = "grey45", deconv = "#2166ac", sem = "#1b7837", imput = "#762a83", combined = "#b2182b")
    lts  <- c(naive = 3, deconv = 2, sem = 4, imput = 5, combined = 1)
    ## a logarithmic axis on the table's own quantity, so the ticks are Table 2's numbers
    ## the y-range comes from the data, not from an imposed floor, so every column of
    ## Table 2 is visible on the panel; a fixed floor clipped the two finest deconv and
    ## imput values out of view in the v10.7 draft
    yv <- 1e3 * cur[is.finite(cur) & cur > 0]
    plot(range(Ds), range(yv), type = "n", log = "y",
         xlab = "imposed grid D (lb)",
         ylab = expression(ISE %*% 10^3), main = "(A) robustness to coarsening")
    for (k in names(cols)) if (!all(is.na(cur[, k])))
      lines(Ds, 1e3 * cur[, k], col = cols[k], lty = lts[k], lwd = 2)
    D <- 30; y <- D * round(wt / D)
    plot(grid, fref, type = "l", lwd = 2.2, xlim = c(90, 320), xlab = "weight (lb)",
         ylab = "density", main = "(B) weight rounded to 30 lb")
    lines(grid, naive_kde(y, grid), lty = 3, lwd = 1.6, col = "grey45")
    lines(grid, as.numeric(adkde(y, D, grid)), lwd = 1.8, col = "#b2182b")
    ## drop any method that produced nothing, so a run without Kernelheaping does not show
    ## a legend entry with no curve behind it
    labA <- c(naive = "naive KDE", deconv = "deconvolution", sem = "SEM (Kernelheaping)",
              imput = "imput (Heitjan-Rubin)", combined = "combined de-heaping")
    keepA <- names(labA)[!apply(is.na(cur[, names(labA), drop = FALSE]), 2, all)]
    legend_cell(unname(labA[keepA]), cols[keepA], lts[keepA], ncol = 2, cex = 0.58)
    legend_cell(c("no-heaping reference", "naive KDE", "combined de-heaping"),
                c("black", "grey45", "#b2182b"), c(1, 3, 1), ncol = 1, cex = 0.6)
    dev.off()
    cat("  READ THIS AS: panel (A) now plots the same quantity Table 2 prints, on a log axis,\n")
    cat("  so a tick of 0.01 in the figure is 0.010 in the table. The stored rds carries every\n")
    cat("  value at full precision, so Table 2 can be emitted rather than typed.\n")
    cat(sprintf("  -> %s\n  -> %s\n", file.path(FIG, "fig_nhanes_weight.pdf"),
                file.path(OUT, "nhanes_v10_7.rds")))
  }
}

cat(sprintf("\nelapsed %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
