#' adheaping: characteristic-function de-heaping density estimation
#'
#' Tuning-free kernel density estimation for heaped and rounded data. See the vignette
#' \code{vignette("adheaping")}. The methods derive from the spectral-decomposition kernel
#' density estimation of Thornton (arXiv:2606.15450). The main entry point is \code{\link{adkde}}.
#'
#' @keywords internal
#' @importFrom stats fft sd dnorm runif quantile median var approx filter
"_PACKAGE"
