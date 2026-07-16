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
