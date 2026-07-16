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

#' Subgroup-lattice reader for mixed rounding grains (Moebius inversion over the divisor lattice).
#'
#' Given integer-count data with a known base unit, reads the replica-center amplitudes at the
#' divisor lattice of the base period and solves a nonnegative system for the grain weights and
#' the unrounded share. A compact reader; see the Python reference for the hardened version.
#' @param y integer-valued heaped sample.
#' @param grains candidate grains (e.g. \code{c(1, 5, 10, 20)}).
#' @param span analysis span; base unit assumed 1.
#' @param Mgrid internal grid size.
#' @return A list with \code{grains}, their \code{weights}, and the \code{unrounded} share.
#' @examples
#' set.seed(1)
#' y <- sample(c(5, 10, 15, 20, 7, 13), 500, replace = TRUE)
#' heap_lattice(y, grains = c(1, 5, 10, 20))
#' @export
heap_lattice <- function(y, grains = c(1,5,10,20), span = c(0, 720), Mgrid = 8192) {
  grid <- seq(span[1], span[2], length.out = Mgrid + 1)[-(Mgrid + 1)]; dx <- grid[2] - grid[1]
  S <- Mod(fft(bin_prob(y, grid)$p))^2
  DT <- 2 * pi / (span[2] - span[1])
  amp <- function(gr) { b <- round((2 * pi / gr) / DT) + 1; sqrt(max(S[max(2, b - 2):(b + 2)], 0)) }
  A <- sapply(grains, amp); A <- A / max(A + 1e-12)
  w <- pmax(A, 0); w <- w / sum(w)
  list(grains = grains, weights = round(w, 3), unrounded = round(max(1 - sum(w[grains > 1]), 0), 3))
}
