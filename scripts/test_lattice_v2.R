#!/usr/bin/env Rscript
## test_lattice_v2.R
## Decision 1 deliverable. Validates the repaired mixed-grain reader, compares it against
## the shipped one, reruns B5 on the repaired reader, and reruns the B6 attenuation test
## on the repaired reader.
##
## The B6 verdict you saw was computed with the shipped reader, which we now know is not
## the algorithm Methods describes. So the attenuation question is still open and is
## re-asked here against the repaired reader. Do not withdraw the attenuation sentence on
## the strength of the earlier run.
##
## Four parts:
##   1  unit tests on the base-R nonnegative least squares
##   2  known-truth simulation, shipped reader against repaired reader
##   3  both readers on the real NHANES counts, against the numbers in the paper
##   4  bootstrap standard errors from the repaired reader, replacing B5
##
## Usage (bash, not zsh):
##     NHANES_DIR=/Users/mitch/src/KDE-AD-HEAPING/data/NHANES bash run_test_lattice_v2.sh
##
## Environment: PKG_DIR, LATTICE_V2 (path to heap_lattice_v2.R), NHANES_DIR,
##              NSIM (default 200), NBOOT (default 1000), OUTDIR.

t0 <- Sys.time()
.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) dirname(normalizePath(f)) else "." }
HERE <- .sd()
pkg <- Sys.getenv("PKG_DIR", unset = "/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping")
if (!dir.exists(file.path(pkg, "R"))) stop("no R/ under PKG_DIR=", pkg)
for (f in list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
v2 <- Sys.getenv("LATTICE_V2", unset = file.path(HERE, "heap_lattice_v2.R"))
if (!file.exists(v2)) stop("heap_lattice_v2.R not found at ", v2, " (set LATTICE_V2)")
source(v2)

NH    <- Sys.getenv("NHANES_DIR", unset = "")
NSIM  <- as.integer(Sys.getenv("NSIM", "200"))
NBOOT <- as.integer(Sys.getenv("NBOOT", "1000"))
OUT   <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
bar <- function(s) cat("\n", strrep("=", 76), "\n ", s, "\n", strrep("=", 76), "\n", sep = "")
res <- list()

GR <- c(1, 5, 10, 20)
TRUEW <- c(0.41, 0.18, 0.12, 0.16); TRUEU <- 0.13
PAPER <- c(0.41, 0.18, 0.12, 0.16)   # the weights printed in the manuscript
PAPER_U <- 0.13

## ---------------------------------------------------------------- 1, unit tests
bar("1  nonnegative least squares, unit tests")
set.seed(1)
X <- matrix(runif(60), 15, 4); bt <- c(0.4, 0, 0.25, 0.1)
e1 <- max(abs(.nnls(X, as.vector(X %*% bt)) - bt))
X2 <- matrix(rnorm(60), 15, 4); y2 <- as.vector(X2 %*% c(1, -1, 2, 0) + rnorm(15, 0, .01))
w2 <- .nnls(X2, y2)
cat(sprintf("  exact recovery of a nonnegative solution : max error %.2e   %s\n",
            e1, if (e1 < 1e-8) "pass" else "FAIL"))
cat(sprintf("  nonnegativity under an indefinite target : min weight %.3e   %s\n",
            min(w2), if (min(w2) >= -1e-12) "pass" else "FAIL"))
res$unit <- c(exact_err = e1, min_w = min(w2))

## --------------------------------------------- 2, known-truth simulation
bar("2  known-truth simulation, shipped reader against repaired reader")
wf <- c(TRUEW, TRUEU) / sum(c(TRUEW, TRUEU)); nn <- 1019
mkbase <- function(pool, m) {
  if (is.null(pool)) return(pmax(1, round(rgamma(m, shape = 2.2, scale = 7))))
  ag <- ifelse(pool %% 20 == 0, 20, ifelse(pool %% 10 == 0, 10, ifelse(pool %% 5 == 0, 5, 1)))
  dj <- pmax(0.5, pool + runif(length(pool), -ag/2, ag/2))
  mu <- mean(dj); v <- stats::var(dj)
  pmax(0.5, rgamma(m, shape = mu^2/v, scale = v/mu))
}
cig <- NULL
if (nzchar(NH) && file.exists(file.path(NH, "SMQ_J.xpt"))) {
  suppressWarnings(suppressMessages(library(foreign)))
  smq <- read.xport(file.path(NH, "SMQ_J.xpt"))
  cn <- intersect(c("SMD650","SMD641","SMQ020"), names(smq))
  if (length(cn)) { v <- suppressWarnings(as.numeric(smq[[cn[1]]]))
    cig <- v[is.finite(v) & v > 0 & v < 200] }
}
cat(sprintf("  base: %s\n\n", if (is.null(cig)) "gamma stand-in, NHANES_DIR not set"
            else sprintf("gamma fitted to %d dequantized real counts", length(cig))))
OM <- NM <- matrix(NA_real_, NSIM, 4); OU <- NU <- numeric(NSIM); SH <- matrix(NA_real_, NSIM, 3)
for (s in seq_len(NSIM)) {
  set.seed(20260627 + s)
  base <- mkbase(cig, nn)
  comp <- sample(seq_along(wf), nn, TRUE, wf)
  y <- pmax(1, round(base))
  for (j in seq_along(GR)) { g <- GR[j]; sel <- comp == j
    if (g > 1) y[sel] <- pmax(g, g * round(base[sel]/g)) }
  o <- tryCatch(heap_lattice(y, grains = GR), error = function(e) NULL)
  if (!is.null(o)) { OM[s, ] <- o$weights; OU[s] <- o$unrounded }
  nw <- tryCatch(heap_lattice2(y, grains = GR), error = function(e) NULL)
  if (!is.null(nw) && nw$fired) { NM[s, ] <- nw$weights; NU[s] <- nw$residual }
  SH[s, ] <- c(mean(y %% 5 == 0), mean(y %% 10 == 0), mean(y %% 20 == 0))
}
ms <- colMeans(SH, na.rm = TRUE)
gate <- all(abs(ms - c(0.570, 0.430, 0.220)) < 0.06)
cat(sprintf("  simulated shares 5/10/20 : %.3f %.3f %.3f   %s\n\n", ms[1], ms[2], ms[3],
            if (gate) "matches the real data, usable" else "does NOT match, treat with caution"))
cat(sprintf("  %-10s %8s | %9s %7s | %9s %7s\n", "quantity", "true", "shipped", "ratio", "repaired", "ratio"))
for (j in seq_along(GR))
  cat(sprintf("  grain %-4d %8.3f | %9.3f %7.2f | %9.3f %7.2f\n", GR[j], TRUEW[j],
      mean(OM[,j], na.rm=TRUE), mean(OM[,j], na.rm=TRUE)/TRUEW[j],
      mean(NM[,j], na.rm=TRUE), mean(NM[,j], na.rm=TRUE)/TRUEW[j]))
cat(sprintf("  %-10s %8.3f | %9.3f %7.2f | %9.3f %7.2f\n", "residual", TRUEU,
    mean(OU, na.rm=TRUE), mean(OU, na.rm=TRUE)/TRUEU,
    mean(NU, na.rm=TRUE), mean(NU, na.rm=TRUE)/TRUEU))
cat(sprintf("\n  shipped residual identical to its grain-one weight : %s\n",
            isTRUE(all.equal(OU, OM[,1]))))
cat(sprintf("  repaired residual identical to its grain-one weight: %s  (correlation %.3f)\n",
            isTRUE(all.equal(NU, NM[,1])), suppressWarnings(cor(NU, NM[,1], use="complete.obs"))))
rat <- sapply(seq_along(GR), function(j) mean(NM[,j], na.rm=TRUE)/TRUEW[j])
cat(sprintf("\n  ATTENUATION TEST on the repaired reader. Ratios: %s\n",
            paste(sprintf("%.2f", rat), collapse = " ")))
cat(sprintf("  %s\n", if (all(rat < 1))
  "All below one and in the same direction: the attenuation account in the paper is supported."
  else "They straddle one: the attenuation account is NOT supported and the sentence should go."))
res$sim <- list(shares = ms, gate = gate, old_w = OM, old_u = OU, new_w = NM, new_u = NU,
                ratios = rat, nsim = NSIM)

## ------------------------------------------- 3, both readers on the real counts
if (!is.null(cig)) {
  bar("3  both readers on the real NHANES counts, against the published numbers")
  o <- heap_lattice(cig, grains = GR)
  nw <- heap_lattice2(cig, grains = GR)
  cat(sprintf("  n = %d\n\n", length(cig)))
  cat(sprintf("  %-10s %10s %10s %10s\n", "quantity", "paper", "shipped", "repaired"))
  for (j in seq_along(GR))
    cat(sprintf("  grain %-4d %10.3f %10.3f %10.3f\n", GR[j], PAPER[j], o$weights[j], nw$weights[j]))
  cat(sprintf("  %-10s %10.3f %10.3f %10.3f\n", "residual", PAPER_U, o$unrounded, nw$residual))
  cat(sprintf("\n  detection statistic lambda from the repaired reader: %.2f (fired: %s)\n",
              nw$lambda, nw$fired))
  cat("\n  The paper's row is Python output. If the repaired R column now lands near it,\n")
  cat("  the one-language claim in Methods becomes true and the published numbers stand.\n")
  cat("  If it does not, the difference has to be reconciled before v10.3 is cut.\n")
  res$real <- list(n = length(cig), paper = c(PAPER, PAPER_U),
                   shipped = c(o$weights, o$unrounded),
                   repaired = c(nw$weights, nw$residual), lambda = nw$lambda)

  ## ------------------------------- 4, bootstrap on the repaired reader
  bar("4  bootstrap standard errors from the repaired reader, replacing B5")
  BW <- matrix(NA_real_, NBOOT, 4); BU <- numeric(NBOOT)
  set.seed(20260627)
  for (bq in seq_len(NBOOT)) {
    r <- tryCatch(heap_lattice2(sample(cig, length(cig), TRUE), grains = GR),
                  error = function(e) NULL)
    if (!is.null(r) && r$fired) { BW[bq, ] <- r$weights; BU[bq] <- r$residual }
  }
  cat(sprintf("  %-10s %10s %10s %22s\n", "quantity", "estimate", "boot SE", "95% percentile CI"))
  for (j in seq_along(GR)) {
    q <- stats::quantile(BW[,j], c(.025,.975), na.rm = TRUE)
    cat(sprintf("  grain %-4d %10.3f %10.3f %10.3f to %8.3f\n", GR[j], nw$weights[j],
                stats::sd(BW[,j], na.rm=TRUE), q[1], q[2]))
  }
  q <- stats::quantile(BU, c(.025,.975), na.rm = TRUE)
  cat(sprintf("  %-10s %10.3f %10.3f %10.3f to %8.3f\n", "residual", nw$residual,
              stats::sd(BU, na.rm=TRUE), q[1], q[2]))
  d <- BW[,3] - BW[,4]
  cat(sprintf("\n  ten-grain minus twenty-grain: %.3f, boot SE %.3f, %.1f SE from zero\n",
              mean(d, na.rm=TRUE), stats::sd(d, na.rm=TRUE),
              abs(mean(d, na.rm=TRUE))/stats::sd(d, na.rm=TRUE)))
  cat(sprintf("  correlation(residual, grain-one weight) across resamples: %.3f\n",
              suppressWarnings(cor(BU, BW[,1], use="complete.obs"))))
  cat("  That correlation was exactly 1.000 with the shipped reader. Anything well below\n")
  cat("  one confirms the residual is now a fitted quantity.\n")
  res$boot <- list(w = BW, u = BU, nboot = NBOOT)
} else {
  bar("3 and 4 SKIPPED")
  cat("  Set NHANES_DIR to the directory holding SMQ_J.xpt to run them.\n")
}

saveRDS(list(provenance = list(script = "test_lattice_v2.R",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = R.version.string, platform = R.version$platform,
  pkg_dir = normalizePath(pkg), lattice_v2 = normalizePath(v2),
  nsim = NSIM, nboot = NBOOT, grains = GR,
  seed_rule = "20260627 + index"), results = res),
  file.path(OUT, "test_lattice_v2.rds"))
cat(sprintf("\n-> %s\nelapsed %.1f s\n", file.path(OUT, "test_lattice_v2.rds"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
