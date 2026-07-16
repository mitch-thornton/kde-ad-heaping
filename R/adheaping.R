## adheaping -- characteristic-function de-heaping density estimation (base R).
## Companion R implementation of the algebraic-diversity de-heaping methods (Thornton,
## arXiv:2606.15450). Pure base R (uses stats::fft); no external dependencies. The Python
## scripts in the bundle are the canonical reference; this port reproduces them.
## Seed of record: 20260627.

## ---- FFT / spectral plumbing ------------------------------------------------

# Angular FFT frequencies (numpy fftfreq convention), times 2*pi.
.fftw <- function(M, dx) {
  k <- 0:(M - 1)
  f <- ifelse(k < (M + 1) %/% 2, k, k - M) / (M * dx)
  2 * pi * f
}

.sinc <- function(u) { out <- sin(u) / u; out[u == 0] <- 1; out }

#' Bin a sample to a probability vector on an equal grid
#'
#' @param x numeric sample.
#' @param grid numeric vector of equally spaced cell centers.
#' @return A list with elements \code{p} (cell probabilities) and \code{dx} (cell width).
#' @examples
#' grid <- seq(-5, 5, length.out = 256)
#' b <- bin_prob(rnorm(500), grid)
#' sum(b$p)
#' @export
bin_prob <- function(x, grid) {
  M <- length(grid); dx <- grid[2] - grid[1]
  idx <- pmin(pmax(floor((x - grid[1]) / dx + 0.5) + 1, 1), M)
  p <- tabulate(idx, nbins = M) / length(x)
  list(p = p, dx = dx)
}

.reconstruct <- function(P, H, dx) {
  f <- Re(fft(P * H, inverse = TRUE)) / length(P) / dx
  pmax(f, 0)
}
.renorm <- function(f, dx) { Z <- sum((f[-1] + f[-length(f)]) / 2) * dx; if (Z > 0) f / Z else f }

.residue_floor <- function(S, w, frac = 0.5) {
  wa <- abs(w); thr <- quantile(wa, 1 - frac)
  stats::median(S[wa >= thr]) / log(2)
}

## ---- core estimators --------------------------------------------------------

#' Naive Gaussian kernel density estimate on a grid
#'
#' Silverman-bandwidth Gaussian kernel density estimate, shown for reference. On heaped
#' data it inherits the rounding comb.
#' @param x numeric sample.
#' @param grid evaluation grid (equally spaced).
#' @return Density values on \code{grid}.
#' @examples
#' grid <- seq(-6, 6, length.out = 512)
#' f <- naive_kde(rnorm(1000), grid)
#' @export
naive_kde <- function(x, grid) {
  M <- length(grid); b <- bin_prob(x, grid); dx <- b$dx
  w <- .fftw(M, dx); n <- length(x)
  h <- 1.06 * stats::sd(x) * n^(-1/5)
  .renorm(.reconstruct(fft(b$p), exp(-0.5 * (h * w)^2), dx), dx)
}

#' Tuning-free de-heaping density estimator (box deconvolution + residue-floor Wiener).
#'
#' Rounding to grid \code{D} is convolution of the density with a width-D box followed
#' by lattice sampling; in the identifiable band |w| < pi/D the empirical characteristic
#' function is divided by sinc(wD/2) (de-Sheppard) and Wiener-tapered against the
#' sampling floor 1/n, amplified by the deconvolution. No bandwidth is chosen by hand.
#' @param y heaped (rounded) sample.
#' @param D grid width.
#' @param grid evaluation grid (equally spaced).
#' @return Density values on \code{grid}.
#' @examples
#' grid <- seq(-6, 6, length.out = 1024)
#' x <- ifelse(runif(2000) < 0.5, rnorm(2000, -1.2, 0.5), rnorm(2000, 1.2, 0.5))
#' f <- deheap_kde(0.5 * round(x / 0.5), 0.5, grid)
#' @export
deheap_kde <- function(y, D, grid) {
  M <- length(grid); b <- bin_prob(y, grid); dx <- b$dx; n <- length(y)
  P <- fft(b$p); S <- Mod(P)^2
  w <- .fftw(M, dx); sc <- .sinc(w * D / 2); wc <- pi / D
  band <- (abs(w) < wc) & (abs(sc) > 1e-2)
  sc2 <- pmax(sc^2, 1e-4); floor <- (1 / n) / sc2
  Sdh <- ifelse(band, S / sc2, 0)
  strip <- pmax(Sdh - floor, 0)
  g <- ifelse(band, strip / (strip + floor), 0)
  H <- ifelse(band, g / sc, 0)
  .renorm(.reconstruct(P, H, dx), dx)
}

## simple base-R Gaussian-mixture EM with BIC selection (broad-component base) -------
.gmm_em <- function(x, k, iters = 120) {
  n <- length(x); mu <- stats::quantile(x, probs = (seq_len(k)) / (k + 1))
  v <- rep(stats::var(x), k); pi_ <- rep(1 / k, k)
  for (it in seq_len(iters)) {
    dens <- sapply(seq_len(k), function(j) pi_[j] * stats::dnorm(x, mu[j], sqrt(v[j])))
    if (is.null(dim(dens))) dens <- matrix(dens, ncol = k)
    tot <- rowSums(dens) + 1e-300; r <- dens / tot
    Nk <- colSums(r) + 1e-8
    pi_ <- Nk / n; mu <- colSums(r * x) / Nk
    v <- pmax(colSums(r * outer(x, mu, "-")^2) / Nk, 1e-4)
  }
  ll <- sum(log(rowSums(sapply(seq_len(k), function(j) pi_[j] * stats::dnorm(x, mu[j], sqrt(v[j])))) + 1e-300))
  list(pi = pi_, mu = mu, v = v, bic = -2 * ll + (3 * k - 1) * log(n))
}
.broad_base <- function(x, grid, n, cwidth = 1.5, kmax = 6) {
  best <- NULL
  for (k in 1:kmax) { m <- .gmm_em(x, k); if (is.null(best) || m$bic < best$bic) best <- m }
  hs <- 0.9 * min(stats::sd(x), diff(stats::quantile(x, c(.25, .75))) / 1.349) * n^(-1/5)
  broad <- sqrt(best$v) >= cwidth * hs
  if (!any(broad)) broad <- seq_along(best$v) == which.max(best$v)
  rowSums(sapply(which(broad), function(j) best$pi[j] * stats::dnorm(grid, best$mu[j], sqrt(best$v[j]))))
}
.deheap_residual_mass <- function(rm, n, w, D) {
  P <- fft(rm); sc <- .sinc(w * D / 2); band <- (abs(w) < pi / D) & (abs(sc) > 1e-2)
  sc2 <- pmax(sc^2, 1e-4); floor <- (1 / n) / sc2
  Sdh <- ifelse(band, Mod(P)^2 / sc2, 0); strip <- pmax(Sdh - floor, 0)
  g <- ifelse(band, strip / (strip + floor), 0); H <- ifelse(band, g / sc, 0)
  Re(fft(P * H, inverse = TRUE)) / length(P)
}

#' Superposition de-heaping estimator
#'
#' Smooth Bayesian-information-criterion Gaussian-mixture base (broad components) plus a
#' band-limited de-heaped residual on the leftover mass; the base carries the non-identifiable
#' band beyond \code{pi/D}.
#' @param y heaped (rounded) sample.
#' @param D grid width.
#' @param grid evaluation grid (equally spaced).
#' @param seed random seed for the jitter dequantization.
#' @return Density values on \code{grid}.
#' @examples
#' grid <- seq(-6, 6, length.out = 1024)
#' x <- ifelse(runif(2000) < 0.5, rnorm(2000, -1.2, 0.5), rnorm(2000, 1.2, 0.5))
#' f <- superpose_kde(1.0 * round(x / 1.0), 1.0, grid)
#' @export
superpose_kde <- function(y, D, grid, seed = 20260627) {
  set.seed(seed); M <- length(grid); b <- bin_prob(y, grid); dx <- b$dx; n <- length(y)
  x <- y + stats::runif(n, -D / 2, D / 2)
  base <- .broad_base(x, grid, n)
  w <- .fftw(M, dx)
  sharp <- .deheap_residual_mass(b$p - base * dx, n, w, D) / dx
  .renorm(pmax(base + sharp, 0), dx)
}

.band_capacity <- function(y, D, grid) {
  M <- length(grid); b <- bin_prob(y, grid); dx <- b$dx; n <- length(y)
  S <- Mod(fft(b$p))^2; w <- .fftw(M, dx); sc <- .sinc(w * D / 2); wc <- pi / D
  band <- (abs(w) < wc) & (abs(sc) > 1e-2)
  Sdh <- ifelse(band, S / pmax(sc^2, 1e-4), 0); floor <- (1 / n) / pmax(sc^2, 1e-4)
  coh <- pmax(Sdh - floor, 0) * band; tot <- sum(coh)
  if (tot <= 0) return(0)
  sum(coh[band & abs(w) > 0.75 * wc]) / tot
}
.superpose_iter <- function(y, D, grid, passes = 3, seed = 20260627) {
  set.seed(seed); M <- length(grid); b <- bin_prob(y, grid); dx <- b$dx; n <- length(y)
  w <- .fftw(M, dx); x <- y + stats::runif(n, -D / 2, D / 2); f <- NULL
  for (it in seq_len(passes)) {
    base <- .broad_base(x, grid, n)
    sharp <- .deheap_residual_mass(b$p - base * dx, n, w, D) / dx
    f <- .renorm(pmax(base + sharp, 0), dx)
    cdf <- cumsum(f) * dx; cdf <- cdf / cdf[M]
    x <- stats::approx(cdf, grid, xout = stats::runif(n), rule = 2)$y
  }
  f
}

#' Combined AD heaping-KDE: band-capacity gate selects de-heap / superposition / iterated.
#'
#' The single tuning-free estimator. A band-capacity statistic (de-heaped in-band power in
#' the outer quarter of the identifiable band) selects the de-heaping estimator when the
#' density fits inside the band, the superposition when it exceeds the band, and the
#' iterated superposition at the coarsest grids.
#' @param y heaped (rounded) sample.
#' @param D grid width.
#' @param grid evaluation grid (equally spaced).
#' @param rho_lo,rho_hi band-capacity thresholds selecting the component.
#' @param seed random seed passed to the superposition components.
#' @return Density values on \code{grid}, with attributes \code{pick} (the selected component)
#'   and \code{rho} (the band-capacity statistic).
#' @examples
#' grid <- seq(-6, 6, length.out = 1024)
#' x <- ifelse(runif(2000) < 0.5, rnorm(2000, -1.2, 0.5), rnorm(2000, 1.2, 0.5))
#' f <- adkde(0.5 * round(x / 0.5), 0.5, grid)
#' attr(f, "pick")
#' @export
adkde <- function(y, D, grid, rho_lo = 0.06, rho_hi = 0.14, seed = 20260627) {
  rho <- .band_capacity(y, D, grid)
  if (rho < rho_lo) return(structure(deheap_kde(y, D, grid), pick = "deheap", rho = rho))
  if (rho < rho_hi) return(structure(superpose_kde(y, D, grid, seed), pick = "superpose", rho = rho))
  structure(.superpose_iter(y, D, grid, 3, seed), pick = "super-iter", rho = rho)
}

## ---- comb readers -----------------------------------------------------------

#' Recover the rounding grid from the leading spectral replica
#'
#' The leading tooth of the power spectrum sits at \code{2*pi/D}, so the grid is
#' \code{Dhat = 2*pi / w_peak}.
#' @param y heaped (rounded) sample.
#' @param grid evaluation grid (equally spaced).
#' @param near optional grid hint; if supplied the tooth is located in a frequency window
#'   around \code{2*pi/near} (a verification that the comb sits at the expected grid); if NULL
#'   the first strong tooth beyond the signal lobe is taken (fully blind).
#' @return The recovered grid width \code{Dhat} (or NA if no tooth is found).
#' @examples
#' grid <- seq(-6, 6, length.out = 2048)
#' x <- ifelse(runif(4000) < 0.5, rnorm(4000, -1.2, 0.5), rnorm(4000, 1.2, 0.5))
#' heap_grid(0.5 * round(x / 0.5), grid, near = 0.5)
#' @export
heap_grid <- function(y, grid, near = NULL) {
  M <- length(grid); b <- bin_prob(y, grid); dx <- b$dx
  S <- Mod(fft(b$p))^2; w <- .fftw(M, dx)
  half <- 2:(M %/% 2); wp <- w[half]; Sp <- S[half]
  if (!is.null(near)) {
    tgt <- 2 * pi / near; win <- wp > 0.5 * tgt & wp < 1.6 * tgt
    if (!any(win)) return(NA_real_)
    return(2 * pi / wp[win][which.max(Sp[win])])
  }
  sm <- stats::filter(Sp, rep(1, 9) / 9, sides = 2); sm[is.na(sm)] <- Sp[is.na(sm)]
  Smax <- max(Sp); ke <- which(sm < 0.02 * Smax); ke <- if (length(ke)) ke[1] else 1
  cand <- which(seq_along(Sp) > ke & Sp > 0.1 * Smax)
  if (!length(cand)) return(NA_real_)
  2 * pi / wp[cand[1]]
}

#' Estimate the heaped fraction from the comb tooth height
#'
#' Under partial heaping the leading tooth has height about \code{p^2} over the residue floor,
#' so \code{phat = sqrt(S_peak - floor)}.
#' @param y heaped (rounded) sample.
#' @param D grid width.
#' @param grid evaluation grid (equally spaced).
#' @return The estimated heaped fraction in \code{[0, 1]}.
#' @examples
#' grid <- seq(-6, 6, length.out = 2048)
#' x <- rnorm(4000); heaped <- runif(4000) < 0.4
#' y <- x; y[heaped] <- round(x[heaped])
#' heap_fraction(y, 1.0, grid)
#' @export
heap_fraction <- function(y, D, grid) {
  M <- length(grid); b <- bin_prob(y, grid); dx <- b$dx
  S <- Mod(fft(b$p))^2; w <- .fftw(M, dx); bfloor <- .residue_floor(S, w)
  target <- 2 * pi / D; win <- abs(w) > 0.85 * target & abs(w) < 1.15 * target
  min(sqrt(max(max(S[win]) - bfloor, 0)), 1)
}
