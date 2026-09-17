## mw15.R
## The fifteen normal-mixture test densities of Marron and Wand (1992),
## "Exact Mean Integrated Squared Error", Annals of Statistics 20(2), 712-736, Table 1.
##
## Each entry is a list with a name and three equal-length vectors, w (weights, summing
## to one), mu (component means) and sd (component standard deviations). The density is
## sum_j w[j] * dnorm(x, mu[j], sd[j]).
##
## Sourced by benchmark_mw15_v10.R. Self-checks live in verify_mw15.R.

MW15 <- list()

MW15[[1]] <- list(index = 1, name = "Gaussian",
  w = 1, mu = 0, sd = 1)

MW15[[2]] <- list(index = 2, name = "Skewed unimodal",
  w  = c(1/5, 1/5, 3/5),
  mu = c(0, 1/2, 13/12),
  sd = c(1, 2/3, 5/9))

## #3 strongly skewed: sum_{l=0}^{7} (1/8) N( 3*((2/3)^l - 1), (2/3)^(2l) )
local({
  l <- 0:7
  MW15[[3]] <<- list(index = 3, name = "Strongly skewed",
    w = rep(1/8, 8), mu = 3 * ((2/3)^l - 1), sd = (2/3)^l)
})

MW15[[4]] <- list(index = 4, name = "Kurtotic unimodal",
  w  = c(2/3, 1/3), mu = c(0, 0), sd = c(1, 1/10))

MW15[[5]] <- list(index = 5, name = "Outlier",
  w  = c(1/10, 9/10), mu = c(0, 0), sd = c(1, 1/10))

MW15[[6]] <- list(index = 6, name = "Bimodal",
  w  = c(1/2, 1/2), mu = c(-1, 1), sd = c(2/3, 2/3))

MW15[[7]] <- list(index = 7, name = "Separated bimodal",
  w  = c(1/2, 1/2), mu = c(-3/2, 3/2), sd = c(1/2, 1/2))

MW15[[8]] <- list(index = 8, name = "Skewed bimodal",
  w  = c(3/4, 1/4), mu = c(0, 3/2), sd = c(1, 1/3))

MW15[[9]] <- list(index = 9, name = "Trimodal",
  w  = c(9/20, 9/20, 1/10), mu = c(-6/5, 6/5, 0), sd = c(3/5, 3/5, 1/4))

## #10 claw: (1/2) N(0,1) + sum_{l=0}^{4} (1/10) N( l/2 - 1, (1/10)^2 )
local({
  l <- 0:4
  MW15[[10]] <<- list(index = 10, name = "Claw",
    w = c(1/2, rep(1/10, 5)), mu = c(0, l/2 - 1), sd = c(1, rep(1/10, 5)))
})

## #11 double claw:
##   (49/100) N(-1,(2/3)^2) + (49/100) N(1,(2/3)^2)
##   + sum_{l=0}^{6} (1/350) N( (l-3)/2, (1/100)^2 )
local({
  l <- 0:6
  MW15[[11]] <<- list(index = 11, name = "Double claw",
    w  = c(49/100, 49/100, rep(1/350, 7)),
    mu = c(-1, 1, (l - 3)/2),
    sd = c(2/3, 2/3, rep(1/100, 7)))
})

## #12 asymmetric claw:
##   (1/2) N(0,1) + sum_{l=-2}^{2} (2^(1-l)/31) N( l + 1/2, (2^(-l)/10)^2 )
local({
  l <- -2:2
  MW15[[12]] <<- list(index = 12, name = "Asymmetric claw",
    w  = c(1/2, 2^(1 - l)/31),
    mu = c(0, l + 1/2),
    sd = c(1, (2^(-l))/10))
})

## #13 asymmetric double claw:
##   sum_{m=0}^{1} (46/100) N(2m-1,(2/3)^2)
##   + sum_{l=1}^{3} (1/300) N(-l/2,(1/100)^2)
##   + sum_{l=1}^{3} (7/300) N( l/2,(7/100)^2)
local({
  m <- 0:1; l <- 1:3
  MW15[[13]] <<- list(index = 13, name = "Asymmetric double claw",
    w  = c(rep(46/100, 2), rep(1/300, 3), rep(7/300, 3)),
    mu = c(2*m - 1, -l/2, l/2),
    sd = c(rep(2/3, 2), rep(1/100, 3), rep(7/100, 3)))
})

## #14 smooth comb:
##   sum_{l=0}^{5} (2^(5-l)/63) N( (65 - 96*(1/2)^l)/21, (32/63)^2 / 2^(2l) )
local({
  l <- 0:5
  MW15[[14]] <<- list(index = 14, name = "Smooth comb",
    w  = 2^(5 - l)/63,
    mu = (65 - 96 * (1/2)^l)/21,
    sd = (32/63)/2^l)
})

## #15 discrete comb:
##   sum_{l=0}^{2} (2/7) N( (12l-15)/7, (2/7)^2 )
##   + sum_{l=8}^{10} (1/21) N( 2l/7, (1/21)^2 )
local({
  a <- 0:2; b <- 8:10
  MW15[[15]] <<- list(index = 15, name = "Discrete comb",
    w  = c(rep(2/7, 3), rep(1/21, 3)),
    mu = c((12*a - 15)/7, 2*b/7),
    sd = c(rep(2/7, 3), rep(1/21, 3)))
})

names(MW15) <- sapply(MW15, function(d) d$name)

## density, sampler, and characteristic function for a normal mixture
mw_d <- function(x, d) {
  out <- numeric(length(x))
  for (j in seq_along(d$w)) out <- out + d$w[j] * dnorm(x, d$mu[j], d$sd[j])
  out
}
mw_r <- function(n, d) {
  k <- sample(seq_along(d$w), n, TRUE, d$w); rnorm(n, d$mu[k], d$sd[k])
}
mw_phi <- function(t, d) {
  out <- complex(length.out = length(t))
  for (j in seq_along(d$w)) out <- out + d$w[j] * exp(1i * d$mu[j] * t - 0.5 * (d$sd[j] * t)^2)
  out
}

## the smallest component standard deviation, which sets the grid resolution a density needs
mw_min_sd <- function(d) min(d$sd)
