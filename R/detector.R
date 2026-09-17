## detector.R -- higher-order comb detection in the spectral basis (self-contained base R).
##
## The power spectrum s = |phi_hat|^2 of a lattice variable is periodic with period K = M/P
## (P = D/dx cells per grid), hence invariant under a shift by K. The symmetry defect of the
## shift-by-K quotient cyclic group, evaluated on the rank-1 spectral covariance s s^T, has the
## closed form
##       eps_K^2 = 1 - (1/g) * sum_{j=0}^{g-1} a(jK)^2 / a(0)^2,   g = M/K,
## where a(l) = sum_i s_i s_{i+l} is the circular autocorrelation of the (de-meaned) spectrum,
## computed by the Wiener-Khinchin identity a = Re(ifft(|fft(s)|^2)). This is the fourth-order
## (covariance-of-the-power-spectrum) group-matching statistic used in the paper. The detected
## grid is the smallest period whose defect falls into a low band (the largest passing cyclic
## subgroup); on exact data no non-trivial period passes.

# Symmetry defect of the shift-by-K quotient group on the spectral covariance (closed form).
.spectral_defect <- function(s, K) {
  M <- length(s); g <- M %/% K
  a <- Re(fft(Mod(fft(s))^2, inverse = TRUE)) / M     # circular autocorrelation
  lags <- (((0:(g - 1)) * K) %% M) + 1
  sqrt(max(1 - mean((a[lags] / a[1])^2), 0))
}

#' Detect the rounding grid as a group-matched atom in the spectral basis
#'
#' Higher-order comb detection: the symmetry defect of the shift-by-\code{K} quotient group on
#' the covariance of the power spectrum (a fourth-order statistic) is minimized at the true
#' replica period; the smallest period with low defect (the largest passing cyclic subgroup) is
#' the detected grid.
#' @param y heaped sample.
#' @param span numeric length-2 range for the internal analysis grid (defaults from the data).
#' @param Mgrid internal grid size (default 512).
#' @param cand_K candidate replica periods in spectral bins.
#' @return A list with \code{D_hat} (recovered grid, or NA), \code{K_hat}, \code{detected}
#'   (logical), and \code{defects} (symmetry defect per candidate period).
#' @examples
#' set.seed(1)
#' x <- ifelse(runif(4000) < 0.5, rnorm(4000, -1.2, 0.5), rnorm(4000, 1.2, 0.5))
#' heap_detect(0.8 * round(x / 0.8), span = c(-12.8, 12.8))$D_hat
#' @export
heap_detect <- function(y, span = NULL, Mgrid = 512, cand_K = c(8,16,24,32,48,64,96)) {
  if (is.null(span)) { r <- range(y); pad <- 0.15 * diff(r) + 1e-9; span <- c(r[1] - pad, r[2] - pad + diff(r) + 2*pad) }
  grid <- seq(span[1], span[2], length.out = Mgrid + 1)[-(Mgrid + 1)]; dx <- grid[2] - grid[1]
  s <- Mod(fft(bin_prob(y, grid)$p))^2; s <- s - mean(s)
  defs <- sapply(cand_K, function(K) .spectral_defect(s, K))
  dmin <- min(defs); dmax <- max(defs); thr <- dmin + 0.35 * (dmax - dmin)
  passing <- cand_K[defs <= thr]
  detected <- length(passing) > 0 && (dmax - dmin) / (dmax + 1e-12) > 0.15
  Khat <- if (detected) min(passing) else NA
  list(D_hat = if (detected) Mgrid * dx / Khat else NA_real_, K_hat = Khat,
       detected = detected, defects = `names<-`(defs, cand_K))
}

#' Mixed-grain reader by Moebius inversion at exclusive replica centers.
#'
#' Reads the share of an integer-valued sample that was rounded to each of several grains.
#' A class of reports rounded to grain \code{g} is a lattice variable with spacing \code{g},
#' so its characteristic function is exactly one at every multiple of \code{2*pi/g}. At a
#' replica center of \code{g} that is not a center of any finer grain, that class contributes
#' its full weight while each finer class contributes the Fourier transform of its
#' residue-class distribution modulo \code{g}. The Moebius differences of the amplitudes at
#' those exclusive centers are therefore the grain weights. The empirical characteristic
#' function is evaluated directly on the sample with no binning, and the sampling floor
#' \code{1/n} is removed from the squared amplitude.
#'
#' The finer classes' contributions vanish only when the base density is smooth at the scale
#' of the coarsest grain. When it is sharp at that scale they alias into the coarse grain's
#' exclusive centers and the reading is biased, most on the coarsest grain. The regime is set
#' by \code{kappa}, the standard deviation of the dequantized base over the coarsest grain.
#' On a known-truth sweep the reader recovers every weight to within about four percent once
#' \code{kappa} exceeds roughly three, and recovers the unit share, and hence the heaped
#' fraction, at every \code{kappa} tested. Below that the grain split should be read as
#' descriptive rather than as a measurement.
#'
#' @param y integer-valued heaped sample.
#' @param grains candidate grains; must include 1 (e.g. \code{c(1, 5, 10, 20)}).
#' @param nboot bootstrap resamples for standard errors; 0 for none.
#' @param seed seed for the bootstrap.
#' @return A list with \code{grains}, \code{weights} (the weight on grain 1 is the unit
#'   share, so one minus it is the heaped fraction), the exclusive-center \code{amplitudes},
#'   \code{n}, and when \code{nboot > 0} the \code{boot_se} and the resample matrix
#'   \code{boot_w}.
#' @examples
#' set.seed(1)
#' base <- round(rnorm(4000, 100, 30))
#' g <- sample(c(1, 5, 10, 20), 4000, replace = TRUE, prob = c(.4, .3, .2, .1))
#' y <- g * round(base / g)
#' heap_lattice(y, grains = c(1, 5, 10, 20))
#' @export
heap_lattice <- function(y, grains = c(1, 5, 10, 20), nboot = 0, seed = 20260627) {
  y <- as.numeric(y)
  grains <- sort(unique(as.integer(grains)), decreasing = TRUE)
  stopifnot(1 %in% grains)
  read_once <- function(yy) {
    A <- vapply(grains, function(g) {
      if (g == 1) return(1)
      sqrt(mean(.ecf_amp2(yy, .exclusive_centers(g, grains))))
    }, numeric(1))
    names(A) <- grains
    w <- numeric(length(grains)); names(w) <- grains
    for (i in seq_along(grains)) {
      g <- grains[i]
      if (g == 1) { w[i] <- max(1 - sum(w[grains > 1]), 0); next }
      coarser <- grains[grains > g & grains %% g == 0]
      w[i] <- max(A[i] - sum(w[as.character(coarser)]), 0)
    }
    list(A = A, w = w)
  }
  r <- read_once(y)
  out <- list(grains = grains, weights = r$w, amplitudes = r$A, n = length(y))
  if (nboot > 0) {
    set.seed(seed)
    bw <- t(replicate(nboot, read_once(sample(y, replace = TRUE))$w))
    out$boot_se <- apply(bw, 2, stats::sd)
    out$boot_w <- bw
  }
  out
}

## squared modulus of the empirical characteristic function at each frequency, with the
## sampling floor 1/n removed
.ecf_amp2 <- function(y, w) {
  n <- length(y)
  a2 <- vapply(w, function(wi) Mod(mean(exp(1i * wi * y)))^2, numeric(1))
  pmax((a2 - 1 / n) / (1 - 1 / n), 0)
}

## multiples m of 2*pi/g, 1 <= m < g, that are not replica centers of any finer grain
.exclusive_centers <- function(g, grains) {
  others <- grains[grains < g]
  ms <- seq_len(g - 1)
  keep <- vapply(ms, function(m) all((m * others) %% g != 0), logical(1))
  2 * pi * ms[keep] / g
}
