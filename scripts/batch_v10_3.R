#!/usr/bin/env Rscript
## batch_v10_3.R
## Step 3 of the v10.3 plan. Six independent work items in one driver, each writing its
## own results file with a provenance block so no number reaches the paper by hand.
##
##   B5  bootstrap standard errors for the NHANES cigarette grain weights      (needs data)
##   B6  known-grain simulation through the same lattice reader, to test the
##       attenuation explanation for the 14 to 26 percent forward-check shortfall
##   B7  sensitivity of the two band-capacity thresholds on the Table 1 grid
##   B8  deconvolution baseline with a data-driven bandwidth, and with an oracle
##       bandwidth as an upper bound on how much that baseline could improve
##   B9  ablation of integrated squared error against the cutoff fraction of pi/D
##   B11 NHANES sensitivity at several reference bandwidths                    (needs data)
##
## B7 is computed efficiently. The band-capacity statistic and the three component
## estimates do not depend on the thresholds, so they are computed once per cell and
## replicate and every threshold pair is then evaluated from the cached components. A
## full two-dimensional sweep therefore costs no more than a single setting.
##
## Usage (you use zsh, so go through the wrapper):
##     bash run_batch_v10_3.sh                 # everything the available data allows
##     BLOCKS=B7,B9 bash run_batch_v10_3.sh    # a subset
##
## Environment:
##   PKG_DIR     kde-ad-heaping checkout   (default /Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping)
##   NHANES_DIR  directory holding DEMO_J.xpt, BMX_J.xpt, SMQ_J.xpt (B5 and B11 only)
##   BLOCKS      comma list, default B5,B6,B7,B8,B9,B11
##   NSEED       replicates for B7, B8, B9   (default 20)
##   NBOOT       bootstrap resamples for B5  (default 1000)
##   NSIM        replicates for B6           (default 200)
##   OUTDIR      results directory           (default $PKG_DIR/results)

t0 <- Sys.time()
pkg <- Sys.getenv("PKG_DIR", unset = "/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping")
if (!dir.exists(file.path(pkg, "R"))) stop("no R/ under PKG_DIR=", pkg)
for (f in list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

BLOCKS <- strsplit(Sys.getenv("BLOCKS", "B5,B6,B7,B8,B9,B11"), ",")[[1]]
NSEED  <- as.integer(Sys.getenv("NSEED", "20"))
NBOOT  <- as.integer(Sys.getenv("NBOOT", "1000"))
NSIM   <- as.integer(Sys.getenv("NSIM", "200"))
OUT    <- Sys.getenv("OUTDIR", unset = file.path(pkg, "results"))
NH     <- Sys.getenv("NHANES_DIR", unset = "")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
has <- function(b) b %in% BLOCKS
res <- list()

prov <- list(script = "batch_v10_3.R",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = R.version.string, platform = R.version$platform,
  pkg_dir = normalizePath(pkg), blocks = paste(BLOCKS, collapse = ","),
  nseed = NSEED, nboot = NBOOT, nsim = NSIM,
  seed_rule = "20260627 + 1000*replicate + round(100*D)")

TARGETS <- list(
  Gaussian = list(w = 1, mu = 0, sd = 1),
  Bimodal  = list(w = c(.5,.5), mu = c(-1.2,1.2), sd = c(.5,.5)),
  Kurtotic = list(w = c(2/3,1/3), mu = c(0,0), sd = c(1,.1)),
  Skewed   = list(w = rep(.2,5), mu = c(0,.5,1.0833,1.4167,1.6875),
                  sd = c(1,2/3,4/9,8/27,16/81)))
dmix <- function(x,t) rowSums(sapply(seq_along(t$w), function(j) t$w[j]*dnorm(x,t$mu[j],t$sd[j])))
rmix <- function(n,t) { k <- sample(seq_along(t$w), n, TRUE, t$w); rnorm(n, t$mu[k], t$sd[k]) }
M <- 2048; grid <- seq(-10, 10, length.out = M + 1)[-(M + 1)]; dx <- grid[2] - grid[1]
n <- 4000; Ds <- c(0.25, 0.5, 1.0, 1.5)
ise <- function(fh, ft) sum((fh - ft)^2) * dx
banner <- function(s) cat("\n", strrep("=", 78), "\n ", s, "\n", strrep("=", 78), "\n", sep = "")

## ============================================================ B6, attenuation
if (has("B6")) {
  banner("B6  known-grain simulation through the lattice reader")
  cat("Simulates integer counts with KNOWN grain weights at the NHANES cigarette sample\n")
  cat("size, pushes them through the same heap_lattice reader, and compares recovered\n")
  cat("weights and forward-predicted round-number shares against the truth. If the\n")
  cat("shortfall reproduces, the attenuation explanation becomes evidence.\n\n")
  GRAINS <- c(1, 5, 10, 20)
  TRUEW  <- c(0.41, 0.18, 0.12, 0.16); TRUEU <- 0.13
  wfull  <- c(TRUEW, TRUEU) / sum(c(TRUEW, TRUEU))
  nn <- 1019
  ## The base count distribution matters. If the NHANES smoking file is available, the
  ## simulation draws its unrounded base from the real reported counts, which makes this a
  ## test of the reader rather than a test of an invented count shape. Without it a gamma
  ## stand-in is used and the result should be read as indicative only.
  base_pool <- NULL
  if (nzchar(NH) && file.exists(file.path(NH, "SMQ_J.xpt"))) {
    suppressWarnings(suppressMessages(library(foreign)))
    .smq <- read.xport(file.path(NH, "SMQ_J.xpt"))
    .cn <- intersect(c("SMD650","SMD641","SMQ020"), names(.smq))
    if (length(.cn)) {
      v <- suppressWarnings(as.numeric(.smq[[.cn[1]]]))
      v <- v[is.finite(v) & v > 0 & v < 200]
      if (length(v) > 200) base_pool <- v
    }
  }
  cat(if (is.null(base_pool))
        "base counts: gamma stand-in, NHANES_DIR not set. Indicative only.\n\n" else
        sprintf("base counts: resampled from %d real NHANES reports.\n\n", length(base_pool)))
  rec <- matrix(NA_real_, NSIM, length(GRAINS)); unr <- numeric(NSIM)
  sh5 <- numeric(NSIM); sh10 <- numeric(NSIM); sh20 <- numeric(NSIM)
  for (s in seq_len(NSIM)) {
    set.seed(20260627 + s)
    base <- if (is.null(base_pool)) pmax(1, round(rgamma(nn, shape = 2.2, scale = 7)))
            else sample(base_pool, nn, TRUE)
    comp <- sample(seq_along(wfull), nn, TRUE, wfull)
    y <- base
    for (j in seq_along(GRAINS)) {
      g <- GRAINS[j]; sel <- comp == j
      if (g > 1) y[sel] <- pmax(g, g * round(base[sel] / g))
    }
    r <- tryCatch(heap_lattice(y, grains = GRAINS), error = function(e) NULL)
    if (!is.null(r)) { rec[s, ] <- r$weights; unr[s] <- r$unrounded }
    sh5[s]  <- mean(y %%  5 == 0); sh10[s] <- mean(y %% 10 == 0); sh20[s] <- mean(y %% 20 == 0)
  }
  cat(sprintf("%-16s %10s %10s %10s\n", "quantity", "true", "recovered", "ratio"))
  for (j in seq_along(GRAINS))
    cat(sprintf("grain %-10d %10.3f %10.3f %10.3f\n", GRAINS[j], TRUEW[j],
                mean(rec[,j], na.rm=TRUE), mean(rec[,j], na.rm=TRUE)/TRUEW[j]))
  cat(sprintf("%-16s %10.3f %10.3f %10.3f\n", "residual share", TRUEU, mean(unr), mean(unr)/TRUEU))
  cat("\nForward check on the simulated data, where the truth IS known:\n")
  cat(sprintf("  observed share multiple of  5 : %.3f\n", mean(sh5)))
  cat(sprintf("  observed share multiple of 10 : %.3f\n", mean(sh10)))
  cat(sprintf("  observed share multiple of 20 : %.3f\n", mean(sh20)))
  cat("\nThe paper reports implied 0.47 / 0.32 / 0.19 against observed 0.57 / 0.43 / 0.22,\n")
  cat("a uniform shortfall of 14 to 26 percent. Compare the ratios above.\n")
  res$B6 <- list(grains = GRAINS, true_w = TRUEW, true_unrounded = TRUEU,
                 base_source = if (is.null(base_pool)) "gamma stand-in" else "real NHANES counts",
                 recovered = rec, unrounded = unr,
                 obs_shares = c(m5 = mean(sh5), m10 = mean(sh10), m20 = mean(sh20)), nsim = NSIM)
}

## ================================================= B7 and B9, cached components
need_cache <- has("B7") || has("B9")
if (need_cache) {
  banner("caching component estimates for B7 and B9")
  cat(sprintf("16 cells x %d replicates. The band-capacity statistic and the three\n", NSEED))
  cat("components do not depend on the thresholds, so every threshold pair in B7 is then free.\n\n")
  cache <- list()
  for (nm in names(TARGETS)) for (D in Ds) {
    key <- paste(nm, D); tt <- TARGETS[[nm]]; ft <- dmix(grid, tt)
    e_dh <- e_su <- e_si <- rho <- numeric(NSEED); ys <- vector("list", NSEED)
    for (s in seq_len(NSEED)) {
      set.seed(20260627 + 1000*s + round(D*100))
      y <- D * round(rmix(n, tt) / D); ys[[s]] <- y
      e_dh[s] <- ise(deheap_kde(y, D, grid), ft)
      e_su[s] <- ise(superpose_kde(y, D, grid), ft)
      e_si[s] <- ise(.superpose_iter(y, D, grid, 3), ft)
      rho[s]  <- .band_capacity(y, D, grid)
    }
    cache[[key]] <- list(target = nm, D = D, ft = ft, ys = ys,
                         deheap = e_dh, super = e_su, super_iter = e_si, rho = rho)
    cat(sprintf("  %-10s D=%4.2f  rho mean %.4f  deheap %8.3f  super %8.3f  iter %8.3f\n",
                nm, D, mean(rho), mean(e_dh)*1e3, mean(e_su)*1e3, mean(e_si)*1e3))
  }
}

## =============================================== B7, band-capacity threshold sweep
if (has("B7")) {
  banner("B7  sensitivity of the two band-capacity thresholds")
  pick <- function(rho, lo, hi) if (rho < lo) "deheap" else if (rho < hi) "super" else "iter"
  total <- function(lo, hi) {
    tot <- 0
    for (k in names(cache)) { c1 <- cache[[k]]
      for (s in seq_len(NSEED)) {
        p <- pick(c1$rho[s], lo, hi)
        tot <- tot + switch(p, deheap = c1$deheap[s], super = c1$super[s], iter = c1$super_iter[s])
      } }
    tot / (length(cache) * NSEED) * 1e3
  }
  LO <- c(0.02, 0.04, 0.06, 0.08, 0.10, 0.12)
  HI <- c(0.10, 0.12, 0.14, 0.16, 0.20, 0.25)
  base <- total(0.06, 0.14)
  cat(sprintf("Mean ISE x1e3 over all 16 cells at the shipped thresholds (0.06, 0.14): %.4f\n\n", base))
  cat("Full grid, rows are rho_lo and columns rho_hi. Entries are mean ISE x1e3.\n\n")
  cat(sprintf("%8s", "lo\\hi")); for (h in HI) cat(sprintf("%10.2f", h)); cat("\n")
  G <- matrix(NA_real_, length(LO), length(HI), dimnames = list(LO, HI))
  for (i in seq_along(LO)) { cat(sprintf("%8.2f", LO[i]))
    for (j in seq_along(HI)) { v <- if (LO[i] < HI[j]) total(LO[i], HI[j]) else NA
      G[i,j] <- v; cat(if (is.na(v)) sprintf("%10s","-") else sprintf("%10.4f", v)) }
    cat("\n") }
  rng <- range(G, na.rm = TRUE)
  cat(sprintf("\nAcross the whole grid the mean ISE ranges %.4f to %.4f, a spread of %.1f%% of the\n",
              rng[1], rng[2], 100*(rng[2]-rng[1])/rng[1]))
  cat(sprintf("shipped-threshold value. Best cell is %.4f at rho_lo=%s, rho_hi=%s.\n",
              min(G, na.rm=TRUE), rownames(G)[which(G == min(G, na.rm=TRUE), arr.ind=TRUE)[1,1]],
              colnames(G)[which(G == min(G, na.rm=TRUE), arr.ind=TRUE)[1,2]]))
  cat("\nREAD THIS AS: a small spread means the results are insensitive to the thresholds,\n")
  cat("which is what Reviewer 2 says would be a strong point in the method's favour.\n")
  cat("It also settles the title question and the provenance concern in major point 5.\n")
  res$B7 <- list(lo = LO, hi = HI, grid_meanISE_x1e3 = G, shipped = base)
}

## =========================================== B9, ablation against the cutoff fraction
if (has("B9")) {
  banner("B9  ablation of ISE against the cutoff fraction of pi/D")
  deheap_cut <- function(y, D, grid, cut) {
    Mg <- length(grid); b <- bin_prob(y, grid); dxx <- b$dx; nn <- length(y)
    P <- fft(b$p); S <- Mod(P)^2; w <- .fftw(Mg, dxx)
    sc <- .sinc(w * D / 2); wc <- cut * pi / D
    band <- (abs(w) < wc) & (abs(sc) > 1e-2)
    sc2 <- pmax(sc^2, 1e-4); fl <- (1/nn)/sc2
    Sdh <- ifelse(band, S/sc2, 0); strip <- pmax(Sdh - fl, 0)
    g <- ifelse(band, strip/(strip + fl), 0); H <- ifelse(band, g/sc, 0)
    .renorm(.reconstruct(P, H, dxx), dxx)
  }
  CUTS <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  A <- matrix(NA_real_, length(cache), length(CUTS),
              dimnames = list(names(cache), sprintf("%.1f", CUTS)))
  cat(sprintf("%-16s", "cell")); for (cc in CUTS) cat(sprintf("%10.1f", cc)); cat("\n")
  for (k in names(cache)) { c1 <- cache[[k]]
    v <- sapply(CUTS, function(cc)
      mean(sapply(seq_len(NSEED), function(s) ise(deheap_cut(c1$ys[[s]], c1$D, grid, cc), c1$ft))) * 1e3)
    A[k, ] <- v
    cat(sprintf("%-16s", k)); cat(sprintf("%10.3f", v)); cat("\n") }
  cm <- colMeans(A)
  cat(sprintf("\n%-16s", "mean over cells")); cat(sprintf("%10.3f", cm)); cat("\n")
  cat(sprintf("\nBest cutoff fraction on average: %.1f (mean ISE %.3f) against the shipped 1.0 (%.3f).\n",
              CUTS[which.min(cm)], min(cm), cm[length(cm)]))
  res$B9 <- list(cuts = CUTS, ise_x1e3 = A, mean_by_cut = cm)
}

## ======================================= B8, deconv with a data-driven bandwidth
if (has("B8")) {
  banner("B8  deconvolution baseline with a data-driven bandwidth")
  cat("The shipped deconv baseline uses a fixed Silverman Gaussian cutoff. Reviewer 2\n")
  cat("point 10 says that understates it. Two alternatives are run here. The first is\n")
  cat("least-squares cross-validation, which is data driven and reportable. The second is\n")
  cat("an ORACLE bandwidth chosen by minimising the true ISE, which no practitioner could\n")
  cat("use and which therefore bounds how much the baseline could possibly improve.\n\n")
  deconv_h <- function(y, D, grid, h) {
    Mg <- length(grid); b <- bin_prob(y, grid); dxx <- b$dx
    P <- fft(b$p); w <- .fftw(Mg, dxx)
    sc <- .sinc(w * D / 2); band <- (abs(w) < pi/D) & (abs(sc) > 0.1)
    H <- ifelse(band, exp(-0.5*(h*w)^2)/sc, 0)
    list(f = .renorm(.reconstruct(P, H, dxx), dxx), H = H)
  }
  lscv <- function(y, D, grid, hs) {
    nn <- length(y); sapply(hs, function(h) {
      o <- deconv_h(y, D, grid, h); f <- o$f
      K0 <- sum(Re(o$H)) / (length(grid) * dx)          # effective kernel at zero lag
      idx <- pmax(1, pmin(length(grid), round((y - grid[1])/dx) + 1))
      floo <- (nn * f[idx] - K0) / (nn - 1)
      sum(f^2)*dx - 2*mean(floo)
    })
  }
  HS <- exp(seq(log(0.02), log(1.5), length.out = 24))
  cat(sprintf("%-10s %5s %10s %10s %10s %10s %10s\n",
              "target","D","fixed","LSCV","oracle","h fixed","h oracle"))
  B8 <- list()
  for (nm in names(TARGETS)) for (D in Ds) {
    tt <- TARGETS[[nm]]; ft <- dmix(grid, tt)
    ef <- el <- eo <- hf <- ho <- numeric(NSEED)
    for (s in seq_len(NSEED)) {
      set.seed(20260627 + 1000*s + round(D*100))
      y <- D * round(rmix(n, tt) / D)
      hf[s] <- 1.06 * stats::sd(y) * n^(-1/5)
      ef[s] <- ise(deconv_kde(y, D, grid), ft)
      cv <- tryCatch(lscv(y, D, grid, HS), error = function(e) rep(NA_real_, length(HS)))
      hl <- if (all(is.na(cv))) hf[s] else HS[which.min(cv)]
      el[s] <- ise(deconv_h(y, D, grid, hl)$f, ft)
      isev <- sapply(HS, function(h) ise(deconv_h(y, D, grid, h)$f, ft))
      ho[s] <- HS[which.min(isev)]; eo[s] <- min(isev)
    }
    B8[[paste(nm,D)]] <- list(target=nm, D=D, fixed=ef, lscv=el, oracle=eo, h_fixed=hf, h_oracle=ho)
    cat(sprintf("%-10s %5.2f %10.3f %10.3f %10.3f %10.4f %10.4f\n", nm, D,
                mean(ef)*1e3, mean(el)*1e3, mean(eo)*1e3, mean(hf), mean(ho)))
  }
  cat("\nREAD THIS AS: if the oracle column is still worse than the combined estimator in a\n")
  cat("cell, then no bandwidth choice rescues the baseline there and the concern is closed\n")
  cat("for that cell regardless of which selector is used.\n")
  res$B8 <- B8
}

## ========================================== B5 and B11, the NHANES-dependent items
nh_ok <- nzchar(NH) && all(file.exists(file.path(NH, c("DEMO_J.xpt","BMX_J.xpt","SMQ_J.xpt"))))
if ((has("B5") || has("B11")) && !nh_ok) {
  banner("B5 and B11 SKIPPED")
  cat("Set NHANES_DIR to the directory holding DEMO_J.xpt, BMX_J.xpt and SMQ_J.xpt.\n")
  cat("Those three files are the ones data/DATA.md names. Re-run with, for example:\n")
  cat("  NHANES_DIR=/Users/mitch/data/nhanes BLOCKS=B5,B11 bash run_batch_v10_3.sh\n")
}
if (nh_ok) {
  suppressWarnings(suppressMessages(library(foreign)))
  if (has("B5")) {
    banner("B5  bootstrap standard errors for the cigarette grain weights")
    smq <- read.xport(file.path(NH, "SMQ_J.xpt"))
    cn <- intersect(c("SMD650","SMD641","SMQ020"), names(smq))
    cig <- suppressWarnings(as.numeric(smq[[cn[1]]]))
    cig <- cig[is.finite(cig) & cig > 0 & cig < 200]
    cat(sprintf("using column %s, n = %d respondents with a positive daily count\n\n", cn[1], length(cig)))
    GR <- c(1,5,10,20)
    pt <- heap_lattice(cig, grains = GR)
    BW <- matrix(NA_real_, NBOOT, length(GR)); BU <- numeric(NBOOT)
    set.seed(20260627)
    for (bq in seq_len(NBOOT)) {
      r <- tryCatch(heap_lattice(sample(cig, length(cig), TRUE), grains = GR), error = function(e) NULL)
      if (!is.null(r)) { BW[bq, ] <- r$weights; BU[bq] <- r$unrounded }
    }
    cat(sprintf("%-14s %10s %10s %22s\n", "quantity", "estimate", "boot SE", "95% percentile CI"))
    for (j in seq_along(GR)) {
      q <- quantile(BW[,j], c(.025,.975), na.rm=TRUE)
      cat(sprintf("grain %-8d %10.3f %10.3f %10.3f to %8.3f\n",
                  GR[j], pt$weights[j], sd(BW[,j], na.rm=TRUE), q[1], q[2]))
    }
    q <- quantile(BU, c(.025,.975), na.rm=TRUE)
    cat(sprintf("%-14s %10.3f %10.3f %10.3f to %8.3f\n", "residual", pt$unrounded,
                sd(BU, na.rm=TRUE), q[1], q[2]))
    d10 <- BW[,3] - BW[,4]
    cat(sprintf("\nten-grain minus twenty-grain: %.3f, boot SE %.3f, %.1f SE from zero\n",
                mean(d10, na.rm=TRUE), sd(d10, na.rm=TRUE),
                abs(mean(d10, na.rm=TRUE))/sd(d10, na.rm=TRUE)))
    cat("That last line answers whether the ten and twenty grains are distinguishable.\n")
    res$B5 <- list(n = length(cig), point = pt, boot_w = BW, boot_u = BU, nboot = NBOOT, column = cn[1])
  }
  if (has("B11")) {
    banner("B11  NHANES sensitivity to the reference bandwidth")
    demo <- read.xport(file.path(NH,"DEMO_J.xpt"))[, c("SEQN","RIDAGEYR")]
    bmx  <- read.xport(file.path(NH,"BMX_J.xpt"))[, c("SEQN","BMXWT")]
    dd <- merge(demo, bmx, by="SEQN"); dd <- dd[dd$RIDAGEYR >= 20 & !is.na(dd$BMXWT), ]
    wt <- dd$BMXWT * 2.2046226218; nw <- length(wt)
    Mg <- 4096; gw <- seq(min(wt)-40, max(wt)+40, length.out = Mg+1)[-(Mg+1)]
    dxw <- gw[2]-gw[1]; ww <- .fftw(Mg, dxw)
    h0 <- 1.06*sd(wt)*nw^(-1/5)
    MULT <- c(0.5, 0.75, 1.0, 1.5, 2.0)
    cat(sprintf("n = %d, Silverman reference bandwidth h0 = %.4f lb\n\n", nw, h0))
    cat(sprintf("%8s %8s", "D (lb)", "mult")); cat(sprintf("%10s","naive"))
    cat(sprintf("%10s","deconv")); cat(sprintf("%10s","imput")); cat(sprintf("%10s\n","combined"))
    B11 <- list()
    for (m in MULT) {
      href <- m*h0
      fref <- .renorm(.reconstruct(fft(bin_prob(wt, gw)$p), exp(-0.5*(href*ww)^2), dxw), dxw)
      isr <- function(f) sum((f-fref)^2)*dxw
      for (D in c(5,10,20,30,40)) {
        y <- D*round(wt/D)
        v <- c(naive = isr(naive_kde(y, gw)), deconv = isr(deconv_kde(y, D, gw)),
               imput = isr(heitjan_mi(y, D, gw, M = 8)),
               combined = isr(as.numeric(adkde(y, D, gw))))
        B11[[paste(m,D)]] <- v
        cat(sprintf("%8.0f %8.2f", D, m)); cat(sprintf("%10.4f", v*1e3)); cat("\n")
      }
    }
    cat("\nREAD THIS AS: if the ordering of the methods is stable across the multiplier,\n")
    cat("then the Table 2 conclusions do not depend on the reference smoothing choice.\n")
    res$B11 <- list(mult = MULT, h0 = h0, n = nw, values = B11)
  }
}

prov$elapsed_sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
saveRDS(list(provenance = prov, results = res), file.path(OUT, "batch_v10_3.rds"))
cat(sprintf("\n-> %s\nblocks run: %s\nelapsed %.1f s\n", file.path(OUT, "batch_v10_3.rds"),
            paste(names(res), collapse = ", "), prov$elapsed_sec))
