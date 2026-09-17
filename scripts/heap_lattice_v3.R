## heap_lattice_v3.R
## Third mixed-grain reader. Standalone base R, no package dependency.
##
## WHAT CHANGED, AND WHY
##
## The known-truth test of 15 September showed that neither the released reader nor the
## repaired reader (heap_lattice_v2.R) recovers the grain weights at the cigarette sample.
## Working through the characteristic function of each grain class shows why, and it is
## not the design matrix. For an integer sample whose reports rounded to grain g, the
## characteristic function is exactly one at every replica center 2*pi*m/g. At a center
## of the coarsest grain that is NOT a center of any finer grain, the coarser class
## contributes its full weight and every finer class contributes the Fourier transform of
## its residue-class distribution modulo the coarse grain. That contribution vanishes only
## when the finer class is spread evenly across the residue classes, which is the case
## when the base density is smooth at the scale of the coarse grain. When the base is
## sharp at that scale, as cigarette counts with a mean near twelve are at a grain of
## twenty, the finer classes alias into the coarse grain's exclusive centers and the
## inversion is biased. That is the identifiability limit of the paper applied to the
## grain mixture: structure sharper than twice the coarsest grain is not separable.
##
## So the reader is written the way the theory says, with nothing else in it:
##   * the empirical characteristic function is evaluated directly on the integer sample,
##     with no binning, at the exclusive replica centers of each grain;
##   * the sampling floor 1/n is subtracted from the squared amplitude;
##   * the Moebius differences from coarsest to finest give the weights;
##   * the unit share is one minus the sum of the rounding weights (for integer data the
##     unit comb has amplitude exactly one, so the unit class and any unrounded residual
##     are one observational class and are reported as one number).
## This reader is unbiased within simulation error in its regime, and the test script maps the regime. It does not
## make the cigarette decomposition work, because the cigarette data are outside the
## regime, and the test script shows that too.
##
## The regime statistic is kappa = (standard deviation of the base) / (coarsest grain).
## The known-truth sweep in test_lattice_v3.R gives recovery against kappa.

.ecf_amp2 <- function(y, w) {
  ## squared modulus of the empirical characteristic function at each w, floor-corrected
  n <- length(y)
  a2 <- vapply(w, function(wi) Mod(mean(exp(1i * wi * y)))^2, numeric(1))
  pmax((a2 - 1 / n) / (1 - 1 / n), 0)
}

.exclusive_centers <- function(g, grains) {
  ## multiples m of 2*pi/g, 1 <= m < g, that are not replica centers of any FINER grain;
  ## every center of g is automatically a center of each coarser multiple of g, and those
  ## coarser grains are the ones whose weights add to the amplitude read there
  others <- grains[grains < g]
  ms <- seq_len(g - 1)
  keep <- vapply(ms, function(m) all((m * others) %% g != 0), logical(1))
  2 * pi * ms[keep] / g
}

#' Mixed-grain reader by Moebius inversion at exclusive replica centers
#' @param y integer-valued heaped sample
#' @param grains candidate grains; must include 1 and be a set closed enough that every
#'   grain has at least one exclusive center (any set of distinct integers works)
#' @param nboot bootstrap resamples for standard errors, 0 for none
#' @return list with grains, weights (the unit grain's weight is the unit share),
#'   amplitudes A_g at the exclusive centers, boot_se, and boot_w (the resample matrix)
heap_lattice3 <- function(y, grains = c(1, 5, 10, 20), nboot = 0, seed = 20260627) {
  y <- as.numeric(y); grains <- sort(unique(as.integer(grains)), decreasing = TRUE)
  stopifnot(1 %in% grains)
  read_once <- function(yy) {
    A <- vapply(grains, function(g) {
      if (g == 1) return(1)
      sqrt(mean(.ecf_amp2(yy, .exclusive_centers(g, grains))))
    }, numeric(1))
    w <- numeric(length(grains)); names(w) <- grains
    for (i in seq_along(grains)) {
      g <- grains[i]
      if (g == 1) { w[i] <- max(1 - sum(w[grains > 1]), 0); next }
      coarser <- grains[grains > g & grains %% g == 0]
      w[i] <- max(A[i] - sum(w[as.character(coarser)]), 0)
    }
    list(A = setNames(A, grains), w = w)
  }
  r <- read_once(y)
  out <- list(grains = grains, weights = r$w, amplitudes = r$A, n = length(y))
  if (nboot > 0) {
    set.seed(seed)
    bw <- t(replicate(nboot, read_once(sample(y, replace = TRUE))$w))
    out$boot_se <- apply(bw, 2, stats::sd); out$boot_w <- bw
  }
  out
}

#' Model-based refinement, NOT adopted for the paper.
#' Estimates the base by dequantizing each report within the rounding cell of a grain
#' drawn from its posterior under the current weights and base, then refits the weights
#' by a nonnegative least-squares regression of the empirical characteristic function on
#' the characteristic functions of the rounded base, over a dense frequency grid. This is
#' the stochastic-imputation lineage the paper contrasts itself with, and is included only
#' to show that the information is present in the data when a base model is supplied.
heap_lattice_refine <- function(y, grains = c(1, 5, 10, 20), iters = 6, nfreq = 400, seed = 20260627) {
  set.seed(seed)
  y <- as.numeric(y); grains <- sort(unique(as.integer(grains))); n <- length(y)
  wgrid <- seq(0.05, 2 * pi - 0.05, length.out = nfreq)
  ecf_c <- function(x) vapply(wgrid, function(wi) mean(exp(1i * wi * x)), complex(1))
  round_class <- function(x, g) pmax(g, g * round(x / g))
  apparent <- function(v) vapply(v, function(yi) max(grains[yi %% grains == 0]), numeric(1))
  nnls <- function(X, b) {                      # Lawson-Hanson, small problem
    p <- ncol(X); P <- logical(p); x <- rep(0, p); it <- 0L
    repeat {
      wv <- as.vector(crossprod(X, b - X %*% x))
      if (!any(!P) || max(wv[!P]) <= 1e-10 || it > 200) break
      j <- which(!P)[which.max(wv[!P])]; P[j] <- TRUE
      repeat {
        it <- it + 1L; s <- rep(0, p); s[P] <- qr.solve(X[, P, drop = FALSE], b)
        if (min(s[P]) > 1e-10) { x <- s; break }
        neg <- P & s <= 1e-10; alpha <- min(x[neg] / (x[neg] - s[neg])); x <- x + alpha * (s - x)
        P[P & abs(x) < 1e-10] <- FALSE; if (!any(P)) break
      }
    }
    pmax(x, 0)
  }
  fit_w <- function(base) {
    X <- sapply(grains, function(g) ecf_c(round_class(base, g))); b <- ecf_c(y)
    nnls(rbind(Re(X), Im(X)), c(Re(b), Im(b)))
  }
  g0 <- apparent(y); base <- pmax(0.5, y + runif(n, -g0 / 2, g0 / 2))
  traj <- matrix(NA_real_, iters, length(grains), dimnames = list(NULL, grains))
  for (it in seq_len(iters)) {
    w <- fit_w(base); traj[it, ] <- w
    xs <- sort(base)
    newb <- numeric(n)
    for (i in seq_len(n)) {
      yi <- y[i]; cand <- grains[yi %% grains == 0]
      cells <- lapply(cand, function(g) c(if (yi == g) 0 else yi - g / 2, yi + g / 2))
      mass <- vapply(cells, function(cl) sum(xs >= cl[1] & xs < cl[2]), numeric(1))
      pr <- w[match(cand, grains)] * mass + 1e-12; pr <- pr / sum(pr)
      cl <- cells[[sample.int(length(cand), 1, prob = pr)]]
      seg <- xs[xs >= cl[1] & xs < cl[2]]
      newb[i] <- if (length(seg) > 5) sample(seg, 1) + runif(1, -0.25, 0.25) else runif(1, cl[1], cl[2])
    }
    base <- pmax(0.5, newb)
  }
  list(grains = grains, weights = setNames(traj[iters, ], grains), trajectory = traj)
}
