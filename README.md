# adheaping — density estimation for heaped data

<!-- badges: start -->
[![R-CMD-check](https://github.com/mitch-thornton/kde-ad-heaping/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mitch-thornton/kde-ad-heaping/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- CRAN badge (enable after acceptance):
[![CRAN status](https://www.r-pkg.org/badges/version/adheaping)](https://CRAN.R-project.org/package=adheaping) -->
<!-- badges: end -->

Tuning-free kernel density estimation for **heaped and rounded data**. The methods extend the
spectral-decomposition kernel density estimation of Thornton (arXiv:2606.15450) to heaped and
rounded data. Rounding to a grid of width `D` is convolution of the
density with a width-`D` box followed by lattice sampling, so the density is recovered by
deconvolving the known box inside the grid-Nyquist band `|w| < pi/D`; beyond that band it is
not identifiable.

## Install

```r
# from source (this repository is an R package at its root):
# install.packages("remotes")
remotes::install_github("mitch-thornton/kde-ad-heaping")
```

## Quick start

```r
library(adheaping)
set.seed(20260627)
x <- ifelse(runif(4000) < 0.5, rnorm(4000, -1.2, 0.5), rnorm(4000, 1.2, 0.5))
D <- 0.5; y <- D * round(x / D)                # heap to a grid of 0.5
grid <- seq(-6, 6, length.out = 2048)

f  <- adkde(y, D, grid)                         # tuning-free combined de-heaping estimator
attr(f, "pick")                                 # which component the band-capacity gate chose
heap_grid(y, grid, near = D)                    # recover the grid from the comb
heap_detect(y, span = c(-12.8, 12.8))$D_hat     # blind spectral detection
```

See `vignette("adheaping")` for a worked example.

## What is here

- R package (root): `R/`, `man/`, `tests/`, `vignettes/`, `DESCRIPTION`, `NAMESPACE`.
- `python/` — the canonical Python reference scripts (the R library reproduces them to ~3e-16).
- `benchmark.R`, `nhanes_R.R`, `berlin_R.R` — reproduce the three benchmark tables (the last two
  need the public NHANES / Berlin data; see `data/DATA.md`).
- `cran-comments.md`, `NEWS.md`, `.github/workflows/` — CRAN-submission and CI support.

## Methods compared

naive Silverman KDE; the AD de-heaping and superposition components; the combined estimator;
measurement-error deconvolution (`deconv_kde`); the Heitjan-Rubin multiple-imputation replica
(`heitjan_mi`); and the real `Kernelheaping` SEM (`sem_kde`) when that package is installed.

## Citation

Please cite the companion methods article:

> M. A. Thornton (2026). Kernel density estimation by spectral decomposition: adaptive
> tapering, superposition, and a data-driven noise floor. arXiv:2606.15450.

```bibtex
@misc{thornton2026adkde,
  author        = {Thornton, Mitchell A.},
  title         = {Kernel density estimation by spectral decomposition: adaptive tapering, superposition, and a data-driven noise floor},
  year          = {2026},
  eprint        = {2606.15450},
  archivePrefix = {arXiv}
}
```

## License and patents

Code under the MIT License (`LICENSE`). See `PATENTS.md` for the patent notice.

Contact: Mitchell A. Thornton — mitch@smu.edu
