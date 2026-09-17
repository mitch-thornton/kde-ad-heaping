#!/usr/bin/env Rscript
## verify_readers_v10_3.R
## Step 2 of the v10.3 plan. Two questions, both answered by running the shipped package
## unmodified on your own machine.
##
##   Q1. Does heap_grid in blind mode (near = NULL) recover the imposed grid, or does it
##       return a fixed wrong value? The paper reports recoveries of 0.500, 1.000 and
##       2.000 with zero variance across seeds. This checks whether those come from the
##       blind branch or from the verifying branch (near = D).
##
##   Q2. Does heap_detect with its shipped defaults fire on the sixteen Table 1 cells?
##       The paper reports the detector recovering imposed grids; this checks the defaults.
##
## Nothing here modifies the package. It only calls it.
##
## Usage (you use zsh, so go through the wrapper):
##     bash run_verify_readers_v10_3.sh
## or directly:
##     PKG_DIR=/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping \
##     MW15_R=/Users/mitch/src/KDE-AD-HEAPING/mw15.R \
##     Rscript verify_readers_v10_3.R

t0 <- Sys.time()
pkg <- Sys.getenv("PKG_DIR", unset = "/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping")
if (!dir.exists(file.path(pkg, "R"))) stop("no R/ under PKG_DIR=", pkg)
for (f in list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

NSEED <- as.integer(Sys.getenv("NSEED", "10"))
OUT   <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## the four Table 1 targets, same definitions the benchmark uses
TARGETS <- list(
  Gaussian = list(w = 1, mu = 0, sd = 1),
  Bimodal  = list(w = c(.5,.5), mu = c(-1.2,1.2), sd = c(.5,.5)),
  Kurtotic = list(w = c(2/3,1/3), mu = c(0,0), sd = c(1,.1)),
  Skewed   = list(w = rep(.2,5), mu = c(0,.5,1.0833,1.4167,1.6875),
                  sd = c(1,2/3,4/9,8/27,16/81)))
rmix <- function(n, t) { k <- sample(seq_along(t$w), n, TRUE, t$w); rnorm(n, t$mu[k], t$sd[k]) }

M <- 2048; grid <- seq(-10, 10, length.out = M + 1)[-(M + 1)]; dx <- grid[2] - grid[1]
n <- 4000

cat("verify_readers_v10_3\n")
cat(sprintf("  package : %s\n  R       : %s\n  seeds   : %d\n  grid    : M=%d dx=%.6f on [-10,10)\n\n",
            normalizePath(pkg), R.version.string, NSEED, M, dx))

## ---------------------------------------------------------------- Q1, heap_grid
cat("== Q1. heap_grid, blind mode against verifying mode ==\n")
cat("The paper's sentence is: on the bimodal target this recovers grids of 0.5, 1.0 and 2.0\n")
cat("as 0.500, 1.000 and 2.000 with zero variance across seeds.\n\n")
cat(sprintf("%-10s %6s | %10s %10s | %10s %10s | %s\n",
            "target","true D","hinted mean","hinted sd","blind mean","blind sd","blind matches?"))
q1 <- list()
for (nm in names(TARGETS)) {
  for (D in c(0.5, 1.0, 2.0)) {
    h <- numeric(NSEED); b <- numeric(NSEED)
    for (s in seq_len(NSEED)) {
      set.seed(20260627 + 1000*s + round(D*100))
      y <- D * round(rmix(n, TARGETS[[nm]]) / D)
      h[s] <- tryCatch(heap_grid(y, grid, near = D), error = function(e) NA_real_)
      b[s] <- tryCatch(heap_grid(y, grid),            error = function(e) NA_real_)
    }
    ok <- isTRUE(abs(mean(b, na.rm = TRUE) - D) < 0.05 * D)
    q1[[paste(nm, D)]] <- list(target = nm, D = D, hinted = h, blind = b, blind_ok = ok)
    cat(sprintf("%-10s %6.2f | %10.4f %10.4f | %10.4f %10.4f | %s\n", nm, D,
                mean(h, na.rm = TRUE), stats::sd(h), mean(b, na.rm = TRUE), stats::sd(b),
                if (ok) "YES" else "no"))
  }
}
nblind <- sum(sapply(q1, function(r) r$blind_ok))
nhint  <- sum(sapply(q1, function(r) abs(mean(r$hinted, na.rm = TRUE) - r$D) < 0.05 * r$D))
cat(sprintf("\n  blind mode recovers the imposed grid in %d of %d cases\n", nblind, length(q1)))
cat(sprintf("  verifying mode recovers it in %d of %d cases\n", nhint, length(q1)))
cat("\n  READ THIS AS: if the verifying column reproduces the paper's 0.500 / 1.000 / 2.000\n")
cat("  and the blind column does not, then the reported recoveries came from the branch\n")
cat("  that is given the expected grid, and the mixed-grain sentence needs correcting.\n\n")

## ------------------------------------------------------------- Q2, heap_detect
cat("== Q2. heap_detect with shipped defaults, on the sixteen Table 1 cells ==\n\n")
cat(sprintf("%-10s %6s %12s %12s %12s\n", "target", "D", "fire rate", "mean D_hat", "rel. error"))
q2 <- list()
for (nm in names(TARGETS)) {
  for (D in c(0.25, 0.5, 1.0, 1.5)) {
    fired <- logical(NSEED); dh <- numeric(NSEED)
    for (s in seq_len(NSEED)) {
      set.seed(20260627 + 1000*s + round(D*100))
      y <- D * round(rmix(n, TARGETS[[nm]]) / D)
      r <- tryCatch(heap_detect(y), error = function(e) list(detected = FALSE, D_hat = NA_real_))
      fired[s] <- isTRUE(r$detected); dh[s] <- r$D_hat
    }
    fr <- mean(fired); md <- mean(dh, na.rm = TRUE)
    q2[[paste(nm, D)]] <- list(target = nm, D = D, fire_rate = fr, dhat = dh)
    cat(sprintf("%-10s %6.2f %12.2f %12s %12s\n", nm, D, fr,
                if (is.nan(md)) "-" else sprintf("%.3f", md),
                if (is.nan(md)) "-" else sprintf("%+.1f%%", 100*(md - D)/D)))
  }
}
cells_fired <- sum(sapply(q2, function(r) r$fire_rate > 0.5))
cat(sprintf("\n  the detector fires in a majority of replicates in %d of 16 cells\n", cells_fired))
cat("\n  READ THIS AS: the paper reports the detector recovering imposed grids of 0.4, 0.8\n")
cat("  and 1.6. If the defaults fire in few of these sixteen cells, then the reported\n")
cat("  detections used different settings, and Methods has to state which.\n\n")

## -------------------------------------------------------------------- verdict
cat("== verdict ==\n")
v1 <- if (nblind <= 1 && nhint >= 9)
  "CONFIRMED: the reported grid recoveries come from the verifying branch, not blind detection." else
  if (nblind >= 9) "NOT CONFIRMED: blind mode works here. The audit finding does not reproduce." else
  "MIXED: neither branch behaves as the audit described. Send me the table above."
v2 <- if (cells_fired <= 4)
  "CONFIRMED: heap_detect with shipped defaults does not fire on most Table 1 cells." else
  if (cells_fired >= 12) "NOT CONFIRMED: the detector fires broadly here. The audit finding does not reproduce." else
  "MIXED: partial firing. Send me the table above."
cat("  Q1 ", v1, "\n")
cat("  Q2 ", v2, "\n\n")

saveRDS(list(
  provenance = list(script = "verify_readers_v10_3.R",
    generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    r_version = R.version.string, platform = R.version$platform,
    pkg_dir = normalizePath(pkg), nseed = NSEED, n = n, M = M, dx = dx,
    seed_rule = "20260627 + 1000*replicate + round(100*D)"),
  q1_heap_grid = q1, q2_heap_detect = q2,
  verdict = list(q1 = v1, q2 = v2, blind_ok = nblind, hinted_ok = nhint, cells_fired = cells_fired)),
  file.path(OUT, "verify_readers_v10_3.rds"))

cat(sprintf("-> %s\nelapsed %.1f s\n", file.path(OUT, "verify_readers_v10_3.rds"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
