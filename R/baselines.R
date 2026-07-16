## baselines.R -- external heaped-data density estimators for head-to-head comparison.
## Three mainstream families are provided:
##   (1) SEM (Gross-Rendtel): the real Kernelheaping::dheaping package is used when installed
##       (see sem_kde); no replica is needed.
##   (2) Heitjan-Rubin multiple imputation: a faithful base-R replica (heitjan_mi), since no
##       standalone Heitjan KDE package exists; stated as a replica in the manuscript.
##   (3) Measurement-error deconvolution KDE (deconv_kde): the Stefanski-Carroll / rounding-as-
##       uniform-error lineage.

#' Measurement-error deconvolution KDE (rounding as Uniform(-D/2, D/2) error).
#'
#' Divides the empirical characteristic function by sinc(wD/2) inside the grid band with a
#' fixed Silverman Gaussian cutoff and a ridge cap; the classical deconvolution baseline.
#' @param y heaped (rounded) sample.
#' @param D grid width.
#' @param grid evaluation grid (equally spaced).
#' @return Density values on \code{grid}.
#' @examples
#' grid <- seq(-6, 6, length.out = 1024)
#' f <- deconv_kde(0.5 * round(rnorm(2000) / 0.5), 0.5, grid)
#' @export
deconv_kde <- function(y, D, grid) {
  M <- length(grid); b <- bin_prob(y, grid); dx <- b$dx; n <- length(y)
  P <- fft(b$p); w <- .fftw(M, dx)
  h <- 1.06 * stats::sd(y) * n^(-1/5)
  sc <- .sinc(w * D / 2); band <- (abs(w) < pi / D) & (abs(sc) > 0.1)
  H <- ifelse(band, exp(-0.5 * (h * w)^2) / sc, 0)
  .renorm(.reconstruct(P, H, dx), dx)
}

#' Heitjan-Rubin multiple-imputation de-heaping (faithful base-R replica).
#'
#' A faithful replica of the multiple-imputation approach of Heitjan and Rubin (1990) for
#' coarse/heaped data: a pilot density is fit to jitter-dequantized values; then, in each of
#' \code{M} independent imputations, every heaped value is imputed within its rounding cell in
#' proportion to the pilot density and a kernel density estimate is formed on the imputed data;
#' the imputations are averaged. Unlike the stochastic-EM chain of Kernelheaping, the draws are
#' independent multiple imputations, in the Heitjan-Rubin style. No standalone Heitjan KDE
#' package exists, so this replica stands in for it; the manuscript states as much.
#' @param y heaped (rounded) sample.
#' @param D grid width.
#' @param grid evaluation grid (equally spaced).
#' @param M number of imputations (default 10).
#' @param seed random seed.
#' @return Density values on \code{grid}.
#' @references Heitjan, D. F. and Rubin, D. B. (1990) Inference from coarse data via multiple
#'   imputation with application to age heaping. JASA 85, 304-314.
#' @examples
#' grid <- seq(-6, 6, length.out = 512)
#' f <- heitjan_mi(1.0 * round(rnorm(1000)), 1.0, grid, M = 4)
#' @export
heitjan_mi <- function(y, D, grid, M = 10, seed = 20260627) {
  set.seed(seed); G <- length(grid); dx <- grid[2] - grid[1]; n <- length(y)
  w <- .fftw(G, dx)
  # pilot density from a single dequantization
  xj <- y + stats::runif(n, -D / 2, D / 2)
  hp <- 1.06 * stats::sd(xj) * n^(-1/5)
  pilot <- .renorm(.reconstruct(fft(bin_prob(xj, grid)$p), exp(-0.5 * (hp * w)^2), dx), dx)
  pilot <- pmax(pilot, 1e-12)
  cells <- max(1, round(D / dx))
  loi <- pmin(pmax(floor((y - D / 2 - grid[1]) / dx) + 1, 1), G)
  cols <- pmin(pmax(outer(loi, 0:(cells - 1), "+"), 1), G)   # n x cells grid indices
  segc <- matrix(pilot[cols], n, cells)                       # pilot density per cell slot
  csn <- segc / rowSums(segc)                                 # normalized weights
  cscum <- t(apply(csn, 1, cumsum))
  acc <- numeric(G)
  for (m in seq_len(M)) {
    u <- stats::runif(n)
    k <- rowSums(cscum < u) + 1L                              # inverse-CDF slot within cell
    k <- pmin(pmax(k, 1L), cells)
    xi <- grid[pmin(pmax(loi + k - 1L, 1L), G)] + stats::runif(n, 0, dx)
    hi <- 1.06 * stats::sd(xi) * n^(-1/5)
    acc <- acc + .reconstruct(fft(bin_prob(xi, grid)$p), exp(-0.5 * (hi * w)^2), dx)
  }
  .renorm(acc / M, dx)
}

#' SEM de-heaping via the real Kernelheaping package
#'
#' Wrapper around the stochastic expectation-maximization estimator of Gross and Rendtel
#' (the \pkg{Kernelheaping} package), interpolated to \code{grid}. Requires that package;
#' returns \code{NULL} if it is not installed.
#' @param y heaped (rounded) sample.
#' @param D grid width.
#' @param grid evaluation grid (equally spaced).
#' @param burnin,samples SEM burn-in and sample counts passed to \code{Kernelheaping::dheaping}.
#' @return Density values on \code{grid}, or \code{NULL} if \pkg{Kernelheaping} is absent.
#' @examples
#' grid <- seq(-6, 6, length.out = 512)
#' if (requireNamespace("Kernelheaping", quietly = TRUE))
#'   f <- sem_kde(0.5 * round(rnorm(800) / 0.5), 0.5, grid)
#' @export
sem_kde <- function(y, D, grid, burnin = 5, samples = 10) {
  if (!requireNamespace("Kernelheaping", quietly = TRUE)) return(NULL)
  res <- try(suppressWarnings(suppressMessages(
    Kernelheaping::dheaping(y, rounds = c(D), burnin = burnin, samples = samples, bw = "nrd0"))),
    silent = TRUE)
  if (inherits(res, "try-error")) return(NULL)
  f <- stats::approx(res$gridx, res$meanPostDensity, xout = grid, rule = 2)$y
  f[is.na(f)] <- 0
  .renorm(pmax(f, 0), grid[2] - grid[1])
}
