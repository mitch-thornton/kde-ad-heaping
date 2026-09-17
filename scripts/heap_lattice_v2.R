## heap_lattice_v2.R
## Repaired mixed-grain reader for the adheaping package, matching the algorithm the
## paper's Methods section describes and that python/exp_cig_lattice_probe_v1.py
## implements. Drop this into R/detector.R in place of the current heap_lattice, or
## source it after the package to override.
##
## WHAT WAS WRONG WITH THE SHIPPED VERSION
##
## The shipped reader did this:
##     A <- sapply(grains, amp); A <- A / max(A + 1e-12)
##     w <- pmax(A, 0); w <- w / sum(w)
##     unrounded <- max(1 - sum(w[grains > 1]), 0)
## Two defects follow. First, the weights are normalized to sum to one, so
## 1 - sum(w[grains > 1]) is identically w[grains == 1] and the residual share carries
## no information. Second, there is no inversion at all. The weights are normalized
## replica amplitudes, which double counts: a report rounded to twenty also lands on the
## ten lattice and the five lattice, so the raw amplitude at a coarse replica center
## contains the weight of every grain whose period divides it. Undoing that overlap is
## the whole point of the divisor-lattice inversion.
##
## WHAT THIS VERSION DOES
##
## Following the Python reference, amplitudes are read at replica centers, a design
## matrix records which grains contribute to each center, and the weights solve a
## nonnegative least-squares system. That is the Moebius inversion generalized to the
## case where some levels are unobservable. Two details matter and are kept:
##   * amplitudes are read only BEYOND the signal extent, because at these scales the
##     coarsest comb's fundamental sits inside the signal band and is confounded there
##   * the weights are NOT renormalized unless they exceed one, so the residual share
##     1 - sum(w) is a fitted quantity rather than an alias of the unit weight
##
## The base-R nonnegative least squares below is Lawson and Hanson's active-set
## algorithm. For a five-column problem it is exact and takes microseconds, so the
## package gains no new dependency.

## ---- Lawson-Hanson nonnegative least squares, base R -----------------------
.nnls <- function(X, y, tol = 1e-10, maxit = 100L) {
  n <- ncol(X); P <- logical(n); w0 <- rep(0, n)
  x <- rep(0, n); it <- 0L
  repeat {
    r <- y - X %*% x
    w0 <- as.vector(crossprod(X, r))
    if (!any(!P) || max(w0[!P], -Inf) <= tol || it >= maxit) break
    j <- which(!P)[which.max(w0[!P])]
    P[j] <- TRUE
    repeat {
      it <- it + 1L
      s <- rep(0, n)
      Xp <- X[, P, drop = FALSE]
      sp <- tryCatch(qr.solve(Xp, y), error = function(e) rep(0, sum(P)))
      s[P] <- sp
      if (min(s[P]) > tol || it >= maxit) { x <- s; break }
      neg <- P & (s <= tol)
      alpha <- min(x[neg] / (x[neg] - s[neg]))
      x <- x + alpha * (s - x)
      P[P & abs(x) < tol] <- FALSE
      if (!any(P)) break
    }
    if (it >= maxit) break
  }
  pmax(x, 0)
}

#' Subgroup-lattice reader for mixed rounding grains (nonnegative divisor-lattice inversion)
#'
#' Reads replica-center amplitudes at the divisor lattice of the base period and solves a
#' nonnegative least-squares system for the grain weights, which is the Moebius inversion
#' generalized to the case where the coarsest levels are confounded with the signal band
#' and so are excluded from the observation set. The residual share is one minus the sum
#' of the fitted weights and is a fitted quantity, not a restatement of the unit weight.
#' @param y integer-valued heaped sample.
#' @param grains candidate grains, stated a priori from the instrument's natural units.
#' @param span analysis span; wide enough that replica centers beyond the signal extent
#'   are resolvable.
#' @param Mgrid internal grid size.
#' @param wmin weights below this are treated as absent.
#' @return A list with \code{grains}, their \code{weights}, the \code{residual} share, the
#'   detection statistic \code{lambda}, and \code{fired}.
#' @export
heap_lattice2 <- function(y, grains = c(1, 2, 5, 10, 20), span = c(-330, 390),
                          Mgrid = 8192, wmin = 0.08) {
  grid <- seq(span[1], span[2], length.out = Mgrid + 1)[-(Mgrid + 1)]
  b <- bin_prob(y, grid)
  S <- Mod(stats::fft(b$p))^2
  nh <- Mgrid %/% 2
  S <- S[seq_len(nh)]
  logS <- log(S + 1e-300)
  ## rolling median, width 12, edges held
  rollmed <- stats::filter(logS, rep(1, 12)/12, sides = 2)
  rollmed[is.na(rollmed)] <- logS[is.na(rollmed)]
  rollmed <- as.numeric(rollmed)

  amp_at <- function(bn, halfwin = 2) {
    lo <- max(2, bn - halfwin); hi <- min(nh, bn + halfwin)
    sqrt(max(max(S[lo:hi]), 0))
  }

  ## provisional floor, then the signal extent, then the floor re-measured beyond it
  hi1 <- min(700, nh - 1)
  vprov <- as.numeric(stats::quantile(rollmed[40:hi1], 0.20, na.rm = TRUE))
  ext <- 8L
  for (k in 8:(nh - 6)) {
    if (all(rollmed[k:(k + 4)] < vprov + 2.0)) { ext <- k; break }
  }
  lo2 <- min(ext + 10, nh - 11)
  vfloor <- as.numeric(stats::quantile(rollmed[lo2:(nh - 10)], 0.20, na.rm = TRUE))
  noise_amp <- exp(vfloor / 2)

  ## the unit comb must actually be present
  spanw <- span[2] - span[1]
  P1 <- spanw / 1.0
  if (round(P1) + 5 >= nh)
    return(list(grains = grains, weights = rep(NA_real_, length(grains)),
                residual = NA_real_, lambda = 0, fired = FALSE))
  a_unit <- amp_at(as.integer(round(P1)))
  lambda <- 2 * log(max(a_unit, 1e-12)) - vfloor
  if (a_unit < max(wmin, 4 * noise_amp))
    return(list(grains = grains, weights = rep(0, length(grains)),
                residual = 1, lambda = lambda, fired = FALSE))

  ## observation bins: replica multiples beyond the signal extent
  Pg <- setNames(P1 / grains, as.character(grains))
  obs <- list()
  for (gi in seq_along(grains)) {
    per <- Pg[gi]; s <- 1
    while (s * per < nh - 6) {
      bn <- as.integer(round(s * per))
      if (bn > ext + 6) obs[[as.character(bn)]] <- union(obs[[as.character(bn)]], gi)
      s <- s + 1
    }
  }
  bins <- sort(as.integer(names(obs)))
  if (length(bins) < length(grains))
    return(list(grains = grains, weights = rep(NA_real_, length(grains)),
                residual = NA_real_, lambda = lambda, fired = FALSE))

  X <- matrix(0, length(bins), length(grains))
  yv <- numeric(length(bins))
  for (i in seq_along(bins)) {
    bn <- bins[i]; yv[i] <- amp_at(bn)
    for (gi in seq_along(grains)) {
      q <- bn / Pg[gi]
      if (abs(q - round(q)) < 0.02) X[i, gi] <- 1
    }
  }
  w <- .nnls(X, yv)
  w[w <= wmin] <- 0
  if (sum(w) > 1) w <- w / sum(w)
  list(grains = grains, weights = round(w, 3),
       residual = round(max(0, 1 - sum(w)), 3),
       lambda = lambda, fired = TRUE)
}
