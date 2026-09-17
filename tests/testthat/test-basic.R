test_that("de-heaping estimator integrates to one", {
  grid <- seq(-10, 10, length.out = 2048)
  set.seed(20260627)
  x <- ifelse(runif(3000) < 0.5, rnorm(3000, -1.2, 0.5), rnorm(3000, 1.2, 0.5))
  y <- 0.5 * round(x / 0.5)
  f <- deheap_kde(y, 0.5, grid)
  expect_true(all(f >= 0))
  trap <- sum((f[-1] + f[-length(f)]) / 2) * (grid[2] - grid[1])
  expect_equal(trap, 1, tolerance = 1e-6)
})

test_that("combined estimator returns a finite density and a pick attribute", {
  grid <- seq(-10, 10, length.out = 2048)
  set.seed(20260627)
  x <- ifelse(runif(3000) < 0.5, rnorm(3000, -1.2, 0.5), rnorm(3000, 1.2, 0.5))
  f <- adkde(0.5 * round(x / 0.5), 0.5, grid)
  expect_true(is.finite(sum(as.numeric(f))))
  expect_true(attr(f, "pick") %in% c("deheap", "superpose", "super-iter"))
})

test_that("grid recovery and spectral detector work", {
  grid <- seq(-10, 10, length.out = 2048)
  set.seed(20260627)
  x <- ifelse(runif(4000) < 0.5, rnorm(4000, -1.2, 0.5), rnorm(4000, 1.2, 0.5))
  expect_equal(heap_grid(0.5 * round(x / 0.5), grid, near = 0.5), 0.5, tolerance = 0.05)
  d <- heap_detect(0.8 * round(x / 0.8), span = c(-12.8, 12.8))
  expect_true(d$detected)
  expect_equal(d$D_hat, 0.8, tolerance = 0.05)
})

## Helper for the mixed-grain tests: the construction of the known-truth sweep, a gamma
## base with sd three quarters of its mean, so that kappa = sd(base) / 20 is set by the
## mean alone. Absolute tolerances are used throughout, since these are shares.
.mg_sample <- function(mu, n, seed, truth = c(0.4, 0.3, 0.2, 0.1)) {
  set.seed(seed)
  sd <- 0.75 * mu
  base <- pmax(0.5, rgamma(n, shape = mu^2 / sd^2, scale = sd^2 / mu))
  g <- sample(c(1, 5, 10, 20), n, replace = TRUE, prob = truth)
  pmax(g, g * round(base / g))
}

test_that("the mixed-grain reader returns a well formed reading", {
  r <- heap_lattice(.mg_sample(192, 4000, 20260627), grains = c(1, 5, 10, 20))
  expect_named(r, c("grains", "weights", "amplitudes", "n"))
  expect_equal(r$n, 4000L)
  expect_true(all(r$weights >= 0))
  expect_lt(sum(r$weights), 1 + 1e-8)
  expect_equal(unname(r$amplitudes["1"]), 1)
})

test_that("the mixed-grain reader recovers every weight well inside its regime", {
  ## mu = 192 gives kappa = 7.2, comfortably above the threshold of about three that the
  ## documentation states. Averaged over draws, as the published sweep is.
  truth <- c(`1` = 0.4, `5` = 0.3, `10` = 0.2, `20` = 0.1)
  W <- vapply(1:25, function(s) {
    r <- heap_lattice(.mg_sample(192, 4000, 20260627 + s), grains = c(1, 5, 10, 20))
    r$weights[names(truth)]
  }, numeric(4))
  m <- rowMeans(W)
  for (k in names(truth)) expect_lt(abs(m[[k]] - truth[[k]]), 0.03)
})

test_that("the unit share survives below the regime where the grain split does not", {
  ## mu = 12 is the cigarette scale, kappa = 0.45, where the documentation says the grain
  ## split is descriptive only. The unit share, and hence the heaped fraction, still holds.
  W <- vapply(1:25, function(s) {
    r <- heap_lattice(.mg_sample(12, 4000, 20260627 + s), grains = c(1, 5, 10, 20))
    r$weights[c("1", "5", "10", "20")]
  }, numeric(4))
  m <- rowMeans(W)
  expect_lt(abs(m[[1]] - 0.4), 0.05)
  expect_true(all(m >= 0))
})

test_that("the mixed-grain reader reports bootstrap standard errors when asked", {
  r <- heap_lattice(.mg_sample(192, 2000, 20260627), grains = c(1, 5, 10, 20), nboot = 50)
  expect_true(all(is.finite(r$boot_se)))
  expect_equal(dim(r$boot_w), c(50L, 4L))
  expect_true(all(r$boot_se > 0))
})
