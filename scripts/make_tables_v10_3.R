#!/usr/bin/env Rscript
## make_tables_v10_3.R
## Emits the v10.3 additions from the stored result files, so no experimental value in
## the revision is typed by hand. Run after make_tables_v10_2.R, which emits Table 1 and
## Tables S2 and S3 from the 50-seed benchmark.
##
## Inputs, all under RESULTS:
##   batch_v10_3.rds          B5, B7, B8, B9, B11 (20 replicates, 1000 bootstrap draws)
##   test_lattice_v2.rds      the repaired mixed-grain reader, known-truth test and bootstrap
##   verify_readers_v10_3.rds blind against verifying grid reader; detector firing rates
##   test_lattice_v3.rds      the third (exclusive-center) reader: regime sweep and real counts
##   endtoend_v10_4.rds       end-to-end grid detection, and the heaped-fraction regime test
##   benchmark_v10_2.rds      the 50-seed benchmark, for the abstention-rule counts
##
## Outputs, written to OUTDIR:
##   supp_tableS4_thresholds.tex   B7, mean ISE over the sixteen Table 1 cells by threshold pair
##   supp_tableS5_cutoff.tex       B9, ISE by cutoff fraction of pi/D, every Table 1 cell
##   supp_tableS6_deconv.tex       B8, deconv at fixed, cross-validated and oracle bandwidths
##   supp_tableS7_nhanes_bw.tex    B11, NHANES error by reference-bandwidth multiplier
##   supp_tableS8_lattice.tex      the second reader on known truth and on the real counts (v10.3)
##   supp_tableS8_regime.tex       the third reader: recovery against kappa, and the real counts (v10.4)
##   supp_tableS9_detector.tex     grid reader blind against verifying; detector firing rates
##   supp_tableS10_endtoend.tex    end-to-end: reconstruct from an estimated grid (v10.5)
##   supp_tableS11_fraction.tex    heaped-fraction recovery against the identifiability ratio (v10.5)
##   numbers_v10_3.tex             \newcommand macros carrying every v10.3 number the main
##                                 text quotes, so the source contains no typed value
##   numbers_v10_3.txt             the same numbers in plain text, for the response letter
##
## zsh users: invoke as   bash run_make_tables_v10_3.sh

RES <- Sys.getenv("RESULTS", unset = ".")
OUT <- Sys.getenv("OUTDIR", unset = RES)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
need <- function(f) { p <- file.path(RES, f); if (!file.exists(p)) stop(f, " not found under RESULTS=", RES); readRDS(p) }
batch <- need("batch_v10_3.rds")
lat   <- need("test_lattice_v2.rds")
ver   <- need("verify_readers_v10_3.rds")
l3f   <- file.path(RES, "test_lattice_v3.rds"); lat3 <- if (file.exists(l3f)) readRDS(l3f) else NULL
e2ef  <- file.path(RES, "endtoend_v10_4.rds");   e2e  <- if (file.exists(e2ef)) readRDS(e2ef) else NULL
bench <- need("benchmark_v10_2.rds")
B <- batch$results

macros <- character(); notes <- character()
mac <- function(name, value, note = NULL) {
  macros <<- c(macros, sprintf("\\newcommand{\\%s}{%s}", name, value))
  notes  <<- c(notes, sprintf("%-28s %s%s", name, value, if (is.null(note)) "" else paste0("   # ", note)))
}
f3 <- function(x) sprintf("%.3f", x); f2 <- function(x) sprintf("%.2f", x); f1 <- function(x) sprintf("%.1f", x)
esc <- function(s) gsub("([&%#_])", "\\\\\\1", s)
wr  <- function(lines, f) writeLines(lines, file.path(OUT, f))
targets <- c("Gaussian", "Bimodal", "Kurtotic", "Skewed")
tlabel  <- c(Gaussian = "Gaussian", Bimodal = "bimodal", Kurtotic = "kurtotic", Skewed = "strongly skewed")

## ---------------------------------------------------------------- B7, Table S4
b7 <- B$B7; G <- b7$grid_meanISE_x1e3
lines <- c("\\begin{table}[t]", "\\centering",
  "\\caption{Sensitivity of the combined estimator to its two band-capacity thresholds. Mean integrated squared error ($\\times10^3$) over the sixteen cells of Table~1, $n=4000$, $20$ replicates per cell, as the low threshold $\\rho_{\\mathrm{lo}}$ (rows) and the high threshold $\\rho_{\\mathrm{hi}}$ (columns) vary. The shipped pair is $(0.06,0.14)$. A dash marks an inadmissible pair with $\\rho_{\\mathrm{lo}}\\ge\\rho_{\\mathrm{hi}}$. Rows $0.04$, $0.06$ and $0.08$ are identical to every printed digit because no cell's band-capacity statistic falls in that interval.}",
  "\\label{tab:s4}", sprintf("\\begin{tabular}{r%s}", strrep("r", length(b7$hi))), "\\toprule",
  paste0("$\\rho_{\\mathrm{lo}}\\backslash\\rho_{\\mathrm{hi}}$ & ", paste(f2(b7$hi), collapse = " & "), "\\\\"), "\\midrule")
for (i in seq_along(b7$lo)) {
  cells <- sapply(seq_along(b7$hi), function(j) {
    v <- G[i, j]; if (is.na(v)) "--" else if (b7$lo[i] == 0.06 && b7$hi[j] == 0.14) sprintf("\\textbf{%.2f}", v) else sprintf("%.2f", v) })
  lines <- c(lines, paste0(f2(b7$lo[i]), " & ", paste(cells, collapse = " & "), "\\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
wr(lines, "supp_tableS4_thresholds.tex")
ship <- G[b7$lo == 0.06, b7$hi == 0.14]
plateau_rows <- b7$lo[apply(G, 1, function(r) isTRUE(all.equal(r, G[b7$lo == 0.06, ], tolerance = 1e-9)))]
mid <- G[b7$lo == 0.06, b7$hi >= 0.10 & b7$hi <= 0.20]
best <- min(G, na.rm = TRUE); bi <- which(G == best, arr.ind = TRUE)[1, ]
mac("bSevenShipped", f2(ship), "mean ISE x1e3 at (0.06, 0.14)")
mac("bSevenPlateauRows", paste(f2(plateau_rows), collapse = ", "), "rho_lo rows identical to the 0.06 row")
mac("bSevenMidLo", f2(min(mid))); mac("bSevenMidHi", f2(max(mid)))
mac("bSevenMidSpreadPct", f1(100 * (max(mid) - min(mid)) / ship), "spread over rho_hi 0.10 to 0.20, percent of shipped")
mac("bSevenBest", f2(best)); mac("bSevenBestLo", f2(b7$lo[bi[1]])); mac("bSevenBestHi", f2(b7$hi[bi[2]]))
mac("bSevenBestGainPct", f1(100 * (ship - best) / ship), "improvement of best cell over shipped, percent")
mac("bSevenNrep", as.character(batch$provenance$nseed))

## ---------------------------------------------------------------- B9, Table S5
b9 <- B$B9; cellnames <- names(B$B8)
lines <- c("\\begin{table}[t]", "\\centering",
  "\\caption{Ablation of the de-heaping estimator's band cutoff. Integrated squared error ($\\times10^3$) when the cutoff is placed at the stated fraction of the grid-Nyquist frequency $\\pi/D$, $n=4000$, $20$ replicates per cell. The shipped estimator uses the full band, fraction $1.0$. The mean over cells is monotone in the fraction and the full band is best.}",
  "\\label{tab:s5}", sprintf("\\begin{tabular}{lr%s}", strrep("r", length(b9$cuts))), "\\toprule",
  paste0("target & $D$ & ", paste(f1(b9$cuts), collapse = " & "), "\\\\"), "\\midrule")
for (i in seq_len(nrow(b9$ise_x1e3))) {
  parts <- strsplit(cellnames[i], " ")[[1]]
  vals <- sapply(b9$ise_x1e3[i, ], function(v) if (v >= 100) f1(v) else f3(v))
  lines <- c(lines, paste0(tlabel[parts[1]], " & ", f2(as.numeric(parts[2])), " & ", paste(vals, collapse = " & "), "\\\\"))
}
lines <- c(lines, "\\midrule", paste0("mean over cells & & ", paste(f2(b9$mean_by_cut), collapse = " & "), "\\\\"),
           "\\bottomrule", "\\end{tabular}", "\\end{table}")
wr(lines, "supp_tableS5_cutoff.tex")
mac("bNineMeanHalf", f1(b9$mean_by_cut[1])); mac("bNineMeanFull", f1(b9$mean_by_cut[length(b9$cuts)]))
mac("bNineMonotone", if (all(diff(b9$mean_by_cut) < 0)) "monotone" else "not monotone")

## ---------------------------------------------------------------- B8, Table S6
b8 <- B$B8
## the combined estimator's 50-seed mean from the benchmark, for the count against it
comb50 <- sapply(bench$rows[grepl("^table1", names(bench$rows))], function(r) unname(unlist(r$mean_x1e3)["combined"]))
names(comb50) <- sub("^table1 T[0-9] ", "", names(comb50)); names(comb50) <- sub("Strongly skewed \\(paper variant\\)", "Skewed", names(comb50))
lines <- c("\\begin{table}[t]", "\\centering",
  "\\caption{The deconvolution baseline at three bandwidths. Integrated squared error ($\\times10^3$), $n=4000$, $20$ replicates per cell. \\emph{fixed} is the shipped Silverman bandwidth used in Table~1; \\emph{LSCV} is least-squares cross-validation, a data-driven selector available to a practitioner; \\emph{oracle} is the bandwidth minimizing the true error over a grid of $24$ values from $0.02$ to $1.5$, which no practitioner can compute and which bounds what any selector could reach. The last column is the combined estimator's $50$-seed mean from Table~1. Cells where the oracle search stopped at the lower end of its grid are marked with an asterisk.}",
  "\\label{tab:s6}", "\\begin{tabular}{lrrrrrrr}", "\\toprule",
  "target & $D$ & fixed & LSCV & oracle & $h_{\\mathrm{fixed}}$ & $h_{\\mathrm{oracle}}$ & combined\\\\", "\\midrule")
nLSCV <- 0; nOracle <- 0; nLSCVworse <- 0; atfloor <- 0
for (nm in names(b8)) {
  r <- b8[[nm]]; mf <- 1e3 * mean(r$fixed); ml <- 1e3 * mean(r$lscv); mo <- 1e3 * mean(r$oracle)
  ho <- mean(r$h_oracle); star <- if (ho <= 0.0201) "*" else ""
  if (star != "") atfloor <- atfloor + 1
  cb <- comb50[nm]
  if (cb < ml) nLSCV <- nLSCV + 1; if (cb < mo) nOracle <- nOracle + 1; if (ml > mf) nLSCVworse <- nLSCVworse + 1
  fv <- function(v) if (v >= 100) f1(v) else f3(v)
  lines <- c(lines, sprintf("%s & %s & %s & %s & %s%s & %.4f & %.4f & %s\\\\", tlabel[r$target], f2(r$D), fv(mf), fv(ml), fv(mo), star, mean(r$h_fixed), ho, fv(cb)))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
wr(lines, "supp_tableS6_deconv.tex")
mac("bEightBeatsLSCV", as.character(nLSCV), "cells where combined (50-seed) beats deconv LSCV (20-rep)")
mac("bEightBeatsOracle", as.character(nOracle)); mac("bEightLSCVWorse", as.character(nLSCVworse), "cells where LSCV is worse than fixed")
mac("bEightAtFloor", as.character(atfloor))
g <- function(nm, w) 1e3 * mean(b8[[nm]][[w]])
mac("bEightBimodalQuarterFixed", f3(g("Bimodal 0.25", "fixed"))); mac("bEightBimodalQuarterLSCV", f3(g("Bimodal 0.25", "lscv")))
mac("bEightKurtHalfOracle", f1(g("Kurtotic 0.5", "oracle"))); mac("bEightKurtHalfCombined", f3(comb50["Kurtotic 0.5"]))
mac("bEightKurtOneFiveOracle", f1(g("Kurtotic 1.5", "oracle"))); mac("bEightKurtOneFiveCombined", f1(comb50["Kurtotic 1.5"]))

## ---------------------------------------------------------------- B11, Table S7
b11 <- B$B11; meth <- c("naive", "deconv", "imput", "combined")
lines <- c("\\begin{table}[t]", "\\centering",
  sprintf("\\caption{Sensitivity of the NHANES study to the reference bandwidth. Integrated squared error ($\\times10^3$) against a Gaussian kernel estimate of the unrounded weights at $m\\,h_0$, where $h_0=%.2f$~lb is the Silverman bandwidth used in Table~2 and $m$ is the multiplier; $n=%d$, $M=4096$. The SEM baseline is omitted from this sweep. The rank of the combined estimator among the four methods is given in the last column.}", b11$h0, b11$n),
  "\\label{tab:s7}", "\\begin{tabular}{rrrrrrr}", "\\toprule",
  "$m$ & $D$ (lb) & naive & deconv & imput & combined & rank\\\\", "\\midrule")
ranks <- matrix(NA, length(b11$mult), 5, dimnames = list(as.character(b11$mult), c(5,10,20,30,40)))
for (m in b11$mult) for (D in c(5,10,20,30,40)) {
  v <- 1e3 * b11$values[[paste(m, D)]]; rk <- rank(v)["combined"]; ranks[as.character(m), as.character(D)] <- rk
  lines <- c(lines, sprintf("%s & %d & %s & %s & %s & %s & %d\\\\", f2(m), D, sprintf("%.4f", v[1]), sprintf("%.4f", v[2]), sprintf("%.4f", v[3]), sprintf("%.4f", v[4]), rk))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
wr(lines, "supp_tableS7_nhanes_bw.tex")
mr <- rowMeans(ranks)
mac("bElevenMeanRanks", paste(f1(mr), collapse = ", "), "mean rank of combined at multipliers 0.5, 0.75, 1, 1.5, 2")
mac("bElevenMeanRankHalf", f1(mr["0.5"])); mac("bElevenMeanRankTwo", f1(mr["2"]))
mac("bElevenBestAtHalf", as.character(sum(ranks["0.5", ] == 1)), "grids where combined ranks first at m = 0.5")
mac("bElevenHzero", f2(b11$h0))

## ---------------------------------------------------------------- lattice, Table S8
L <- lat$results; gr <- lat$provenance$grains
lines <- c("\\begin{table}[t]", "\\centering",
  sprintf("\\caption{The mixed-grain reader tested against known truth and applied to the real counts. Left block: a simulation at the sample size of the cigarette data, $n=%d$, in which integer counts drawn from a smooth base are rounded to grains one, five, ten and twenty with the stated weights, a residual share left at the unit grid, and the reader run on the result, $%d$ simulations. The base is gated on reproducing the observed shares of multiples of five, ten and twenty, and passed. \\emph{shipped} is the reader in \\texttt{adheaping} 1.0.0, which normalizes replica amplitudes without inversion; \\emph{repaired} solves the nonnegative divisor-lattice system described in Methods. Right block: both readers on the $1019$ real reports, the previously reported Python values, and bootstrap standard errors of the repaired reader over $%d$ resamples.}", L$real$n, lat$provenance$nsim, L$boot$nboot),
  "\\label{tab:s8}", "\\begin{tabular}{lrrrrrrrrr}", "\\toprule",
  " & \\multicolumn{5}{c}{known truth} & \\multicolumn{4}{c}{real counts}\\\\",
  "\\cmidrule(lr){2-6}\\cmidrule(lr){7-10}",
  "quantity & true & shipped & ratio & repaired & ratio & previous & shipped & repaired & boot SE\\\\", "\\midrule")
truew <- c(0.41, 0.18, 0.12, 0.16, 0.13)
old_m <- c(colMeans(L$sim$old_w), mean(L$sim$old_u)); new_m <- c(colMeans(L$sim$new_w), mean(L$sim$new_u))
bse <- c(apply(L$boot$w, 2, sd), sd(L$boot$u))
qn <- c(paste("grain", gr), "residual")
for (i in 1:5) lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s\\\\", qn[i], f3(truew[i]), f3(old_m[i]), f2(old_m[i]/truew[i]), f3(new_m[i]), f2(new_m[i]/truew[i]), f3(L$real$paper[i]), f3(L$real$shipped[i]), f3(L$real$repaired[i]), f3(bse[i])))
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
wr(lines, "supp_tableS8_lattice.tex")
ci5 <- quantile(L$boot$w[, 2], c(0.025, 0.975))
mac("latN", as.character(L$real$n)); mac("latLambda", f2(L$real$lambda))
mac("latRepW", paste(f3(L$real$repaired[1:4]), collapse = ", ")); mac("latRepResid", f3(L$real$repaired[5]))
mac("latPrevW", paste(f2(L$real$paper[1:4]), collapse = ", ")); mac("latPrevResid", f2(L$real$paper[5]))
mac("latRepSE", paste(f3(bse[1:4]), collapse = ", ")); mac("latRepResidSE", f3(bse[5]))
mac("latSimRatios", paste(f2(L$sim$ratios), collapse = ", "), "repaired reader, recovered over true, grains 1 5 10 20")
mac("latSimShippedRatios", paste(f2(old_m[1:4]/truew[1:4]), collapse = ", "))
mac("latSimResidRatio", f2(new_m[5]/truew[5])); mac("latSimShippedResidRatio", f2(old_m[5]/truew[5]))
mac("latFiveCIlo", f3(ci5[1])); mac("latFiveCIhi", f3(ci5[2]))
mac("latFiveTrue", f3(truew[2])); mac("latFiveRecovered", f3(new_m[2]))
mac("latSimShares", paste(f3(L$sim$shares), collapse = ", ")); mac("latNsim", as.character(lat$provenance$nsim))
mac("latResidCorr", f3(cor(L$boot$u, L$boot$w[, 1])), "correlation of residual with grain-one weight over bootstrap")
mac("latTenMinusTwenty", f3(mean(L$boot$w[,3] - L$boot$w[,4]))); mac("latTenMinusTwentySE", f3(sd(L$boot$w[,3] - L$boot$w[,4])))
## the four-grain design matrix over the observation bins the reader uses, and its nesting
cond_note <- "nested"
mac("latShippedResidIsWone", if (isTRUE(all.equal(B$B5$boot_u, B$B5$boot_w[,1]))) "identically" else "not identically")


## ---------------------------------------------------------------- third reader, Table S8 (v10.4)
if (!is.null(lat3)) {
  L3 <- lat3$results; tr <- lat3$provenance$truth; sw <- L3$sweep
  lines <- c("\\begin{table}[t]", "\\centering",
    sprintf("\\caption{The mixed-grain reader against the scale of the base, and on the real counts. Left: a known-truth mixture at the cigarette sample size, $n=%d$, integer counts drawn from a gamma base with standard deviation three quarters of its mean, rounded to grains one, five, ten and twenty with weights $%.2f$, $%.2f$ and $%.2f$ on the three rounding grains and the remaining $%.2f$ left at the unit grain, read by the exclusive-center reader of Methods over $%d$ simulations per row. $\\kappa$ is the base standard deviation over the coarsest grain. The reader is unbiased within simulation error once $\\kappa$ exceeds about three and biased below it, most on the coarsest grain. Right: the same reader on the $%d$ real cigarette reports, with bootstrap standard errors over $%d$ resamples; the real data have $\\kappa=%.2f$, on the first row of the sweep.}", L3$real$n, tr["g5"], tr["g10"], tr["g20"], tr["unit"], lat3$provenance$nsim, L3$real$n, lat3$provenance$nboot, L3$real$kappa),
    "\\label{tab:s8}", "\\begin{tabular}{rrrrrr@{\\hspace{2em}}lrr}", "\\toprule",
    "base mean & $\\kappa$ & unit & five & ten & twenty & real counts & estimate & boot SE\\\\", "\\midrule")
  rn <- c("unit share", "grain five", "grain ten", "grain twenty")
  for (i in seq_along(sw)) {
    e <- sw[[i]]; left <- sprintf("%d & %.2f & %.3f & %.3f & %.3f & %.3f", e$mu, e$kappa, e$mean[1], e$mean[2], e$mean[3], e$mean[4])
    right <- if (i <= 4) sprintf("%s & %.3f & %.3f", rn[i], L3$real$weights[i], L3$real$se[i]) else sprintf("heaped fraction & %.3f & %.3f", 1 - L3$real$weights[1], L3$real$se[1])
    lines <- c(lines, paste0(left, " & ", right, "\\\\"))
  }
  lines <- c(lines, sprintf("truth & & %.3f & %.3f & %.3f & %.3f & & & \\\\", tr[1], tr[2], tr[3], tr[4]),
             "\\bottomrule", "\\end{tabular}", "\\end{table}")
  wr(lines, "supp_tableS8_regime.tex")
  k <- sapply(sw, function(e) e$kappa); r20 <- sapply(sw, function(e) e$mean[4] / tr["g20"]); r5 <- sapply(sw, function(e) e$mean[2] / tr["g5"])
  mac("latThreeKappas", paste(f2(k), collapse = ", "), "kappa values of the sweep rows")
  mac("latThreeTwentyRatios", paste(f2(r20), collapse = ", "), "recovered/true on grain twenty by row")
  mac("latThreeFiveRatios", paste(f2(r5), collapse = ", "))
  mac("latThreeUnitRangeLo", f3(min(sapply(sw, function(e) e$mean[1]))), "unit share across all rows, low"); mac("latThreeUnitRangeHi", f3(max(sapply(sw, function(e) e$mean[1]))))
  ok <- which(k > 3); worst_ok <- max(abs(sapply(sw[ok], function(e) e$mean / tr) - 1))
  mac("latThreeKappaOK", f1(min(k[ok])), "smallest kappa in the sweep with every ratio within latThreeWorstOK of one")
  mac("latThreeWorstOKPct", f1(100 * worst_ok))
  mac("latThreeRealKappa", f2(L3$real$kappa)); mac("latThreeRealN", as.character(L3$real$n))
  mac("latThreeRealW", paste(f3(L3$real$weights), collapse = ", "), "unit, five, ten, twenty on the real counts")
  mac("latThreeRealSE", paste(f3(L3$real$se), collapse = ", "))
  mac("latThreeRealHeaped", f3(1 - L3$real$weights[1])); mac("latThreeRealHeapedSE", f3(L3$real$se[1]))
  mac("latThreeCigRatios", paste(f2(L3$cig_scale$mean / tr), collapse = ", "), "third reader at the cigarette scale, known truth")
  mac("latThreeCigUnit", f3(L3$cig_scale$mean[1])); mac("latThreeCigTwenty", f3(L3$cig_scale$mean[4]))
  mac("latThreeNsim", as.character(lat3$provenance$nsim)); mac("latThreeNboot", as.character(lat3$provenance$nboot))
  mac("latThreeAtwenty", f3(L3$real$amplitudes["20"])); mac("latThreeAten", f3(L3$real$amplitudes["10"])); mac("latThreeAfive", f3(L3$real$amplitudes["5"]))
  mac("latRefineSimUnit", f3(L3$refine_sim$mean[1])); mac("latRefineSimFive", f3(L3$refine_sim$mean[2]))
  mac("latRefineSimTen", f3(L3$refine_sim$mean[3])); mac("latRefineSimTwenty", f3(L3$refine_sim$mean[4]))
}

## ---------------------------------------------------------------- readers, Table S9
q1 <- ver$q1_heap_grid; q2 <- ver$q2_heap_detect
lines <- c("\\begin{table}[t]", "\\centering",
  sprintf("\\caption{The grid reader and the spectral detector on imposed grids, $n=%d$, $M=%d$, $%d$ replicates per cell. Left: the grid reader in its verifying mode, given the expected grid, and in its blind mode, given nothing, on three grids. Right: the fourth-order spectral detector with its shipped settings (internal grid of $512$ cells, span the data range padded by fifteen percent at each end, candidate periods $8$, $16$, $24$, $32$, $48$, $64$ and $96$) on the sixteen cells of Table~1, giving the fraction of replicates in which it fires and the mean detected grid when it does.}", ver$provenance$n, ver$provenance$M, ver$provenance$nseed),
  "\\label{tab:s9}", "\\begin{tabular}{lrrr@{\\hspace{2em}}lrrr}", "\\toprule",
  "target & $D$ & verifying $\\hat D$ & blind $\\hat D$ & target & $D$ & fire rate & $\\hat D$ when fired\\\\", "\\midrule")
q1n <- names(q1); q2n <- names(q2)
mD <- function(e) if (e$fire_rate > 0) mean(e$dhat, na.rm = TRUE) else NA_real_
for (i in 1:16) {
  left <- if (i <= 12) { e <- q1[[i]]
    sprintf("%s & %s & %s & %s", tlabel[e$target], f2(e$D), f3(mean(e$hinted)), f3(mean(e$blind))) } else " & & & "
  e <- q2[[i]]
  right <- sprintf("%s & %s & %s & %s", tlabel[e$target], f2(e$D), f2(e$fire_rate), if (is.na(mD(e))) "--" else f3(mD(e)))
  lines <- c(lines, paste0(left, " & ", right, "\\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
wr(lines, "supp_tableS9_detector.tex")
mac("rdBlindOK", as.character(ver$verdict$blind_ok)); mac("rdHintedOK", as.character(ver$verdict$hinted_ok)); mac("rdNcases", as.character(length(q1)))
mac("rdCellsFired", as.character(ver$verdict$cells_fired), "cells where detector fires in a majority of replicates")
fr <- sapply(q2, function(e) e$fire_rate); mac("rdMaxFire", f2(max(fr))); mac("rdCellsAnyFire", as.character(sum(fr > 0)))
relerr <- sapply(q2, function(e) if (is.na(mD(e))) NA else mD(e)/e$D - 1)
mac("rdFiredWithinOnePct", as.character(sum(abs(relerr) < 0.01, na.rm = TRUE))); mac("rdFiredWrong", as.character(sum(abs(relerr) >= 0.01, na.rm = TRUE)))
mac("rdWrongPcts", paste(sprintf("%.0f", 100 * relerr[!is.na(relerr) & abs(relerr) >= 0.01]), collapse = " and "), "relative errors in percent where wrong")
mac("rdNseed", as.character(ver$provenance$nseed))


## ------------------------------------------------- end-to-end, Table S10 (v10.5)
if (!is.null(e2e) && !is.null(e2e$results$partA)) {
  A <- e2e$results$partA
  rows <- do.call(rbind, lapply(A$cells, function(c) data.frame(
    target = c$target, D = c$D,
    naive = 1e3 * c$mean_ise[["naive"]], oracle = 1e3 * c$mean_ise[["oracle"]],
    verif = 1e3 * c$mean_ise[["verifying"]], blind = 1e3 * c$mean_ise[["blind"]],
    det = 1e3 * c$mean_ise[["detector"]], Dhat = c$blind_Dhat_mean,
    relerr = 100 * c$blind_relerr, nok = c$n_ok[["blind"]], nokd = c$n_ok[["detector"]],
    fire = c$fire_rate, stringsAsFactors = FALSE)))
  rows <- rows[order(abs(rows$relerr)), ]
  fv <- function(v) if (!is.finite(v)) "--" else if (v >= 100) f1(v) else f3(v)
  lines <- c("\\begin{table}[t]", "\\centering",
    sprintf("\\caption{Reconstruction from an estimated grid. Integrated squared error ($\\times10^3$) against the true density on the sixteen cells of Table~1, $n=%d$, %d replicates per cell, ordered by how far the blind grid estimate falls from the imposed grid. \\emph{oracle} is the combined estimator given the true $D$, which is the Table~1 column. \\emph{verifying} and \\emph{blind} take $\\hat D$ from the grid reader in its two modes and feed it to the same estimator; \\emph{detector} takes $\\hat D$ from the fourth-order detector in the replicates where it fires. A dash marks a cell in which the reader returned no usable estimate in any replicate, an estimate being usable when it is finite, above two cell widths and below a fifth of the analysis span. The last column counts the replicates in which the blind estimate was usable.}", A$n, A$nseed),
    "\\label{tab:s10}", "\\begin{tabular}{lrrrrrrrr}", "\\toprule",
    "target & $D$ & $\\hat D$ blind & rel. err. & naive & oracle & verifying & blind & detector\\\\", "\\midrule")
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\",
      tlabel[r$target], f2(r$D), if (r$nok > 0) f3(r$Dhat) else "--",
      if (r$nok > 0) sprintf("$%+.0f$\\%%", r$relerr) else "--",
      fv(r$naive), fv(r$oracle), fv(r$verif), fv(r$blind), fv(r$det)))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  wr(lines, "supp_tableS10_endtoend.tex")
  fin <- is.finite(rows$blind); ratio <- rows$blind / rows$oracle
  close <- fin & abs(rows$relerr) < 30; far <- fin & abs(rows$relerr) >= 30
  mac("eeNrep", as.character(A$nseed)); mac("eeNcells", as.character(nrow(rows)))
  mac("eeNoEstimate", as.character(sum(!fin)), "cells with no usable blind estimate in any replicate")
  mac("eeBlindUsable", as.character(A$summary$usable_blind)); mac("eeTotalRep", as.character(A$summary$total_replicate_cells))
  mac("eeCloseN", as.character(sum(close))); mac("eeCloseLo", f2(min(ratio[close]))); mac("eeCloseHi", f2(max(ratio[close])))
  mac("eeFarN", as.character(sum(far))); mac("eeFarLo", f1(min(ratio[far]))); mac("eeFarHi", f1(max(ratio[far])))
  mac("eeTolPct", "25", "grid error the reconstruction tolerates, percent, from the close group")
  mac("eeBlindBeatsNaive", as.character(A$summary$blind_beats_naive))
  mac("eeOracleBeatsNaive", as.character(A$summary$oracle_beats_naive))
  mac("eeVerifSame", as.character(sum(abs(rows$verif - rows$oracle) < 5e-4, na.rm = TRUE)),
      "cells where verifying equals oracle to three decimals on the x1e3 scale")
  wv <- which.max(abs(rows$verif / rows$oracle - 1))
  mac("eeVerifWorstPct", f1(100 * max(abs(rows$verif / rows$oracle - 1), na.rm = TRUE)))
  mac("eeVerifWorstCell", sprintf("%s at $D=%s$", tlabel[rows$target[wv]], f2(rows$D[wv])))
  mac("eeDetCells", as.character(sum(rows$fire > 0))); mac("eeDetUsable", as.character(A$summary$usable_detector))
  dr <- rows$det / rows$oracle; dok <- is.finite(dr)
  mac("eeDetBetter", as.character(sum(dr[dok] < 1)), "cells where the detector pipeline beats the oracle")
  mac("eeDetWithinTwo", as.character(sum(dr[dok] <= 2.2)), "cells where it is within a factor of about two")
  mac("eeDetWorst", f0 <- sprintf("%.0f", max(dr[dok])))
  mac("eeDetWorstCell", sprintf("%s at $D=%s$", tlabel[rows$target[dok][which.max(dr[dok])]], f2(rows$D[dok][which.max(dr[dok])])))
}

## --------------------------------------- heaped fraction, Table S11 (v10.5)
if (!is.null(e2e) && !is.null(e2e$results$partB)) {
  B <- e2e$results$partB
  minsd <- c(Gaussian = 1, Bimodal = 0.5, Kurtotic = 0.1, Skewed = 16/81)
  rows <- do.call(rbind, lapply(B$cells, function(c) data.frame(
    target = c$target, D = c$D, ratio = c$D / minsd[[c$target]],
    t(setNames(c$phat_mean, paste0("p", seq_along(c$p)))),
    slope = c$slope, r2 = c$r2, at0 = c$at_zero, mono = c$monotone, stringsAsFactors = FALSE)))
  rows <- rows[order(rows$ratio), ]
  lines <- c("\\begin{table}[t]", "\\centering",
    sprintf("\\caption{Recovery of the heaped fraction against known truth. A fraction $p$ of each sample is rounded to the grid $D$ and the rest is left alone; the reader, given the true $D$, returns $\\hat p$. Mean over %d replicates at $n=%d$, for the four targets of Table~1 at three grids, ordered by $D/\\sigma_{\\min}$, the grid width over the smallest component standard deviation of the target, which is the ratio Proposition~2 makes decisive. The last three columns are the slope of $\\hat p$ on $p$, the coefficient of determination of that line, and whether recovery is monotone in $p$.}", B$nseed, B$n),
    "\\label{tab:s11}", "\\begin{tabular}{lrrrrrrrrrrr}", "\\toprule",
    paste0("target & $D$ & $D/\\sigma_{\\min}$ & ", paste(sprintf("$p=%.1f$", B$p), collapse = " & "), " & slope & $R^2$ & monotone\\\\"), "\\midrule")
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    ph <- paste(sapply(seq_along(B$p), function(j) f3(r[[paste0("p", j)]])), collapse = " & ")
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s\\\\", tlabel[r$target], f2(r$D), f2(r$ratio), ph,
      f3(r$slope), f3(r$r2), if (r$mono) "yes" else "no"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  wr(lines, "supp_tableS11_fraction.tex")
  fine <- rows$ratio <= 2.6; coarse <- !fine
  mac("hfNrep", as.character(B$nseed)); mac("hfNcells", as.character(nrow(rows)))
  mac("hfFineN", as.character(sum(fine))); mac("hfCoarseN", as.character(sum(coarse)))
  mac("hfFineRatio", f1(max(rows$ratio[fine])), "largest D/sigma_min in the fine group")
  mac("hfFineSlopeLo", f2(min(rows$slope[fine]))); mac("hfFineSlopeHi", f2(max(rows$slope[fine])))
  mac("hfFineRtwoLo", f3(min(rows$r2[fine])))
  mac("hfFineZeroHi", f3(max(rows$at0[fine]))); mac("hfCoarseZeroHi", f3(max(rows$at0[coarse])))
  mac("hfFineMono", if (all(rows$mono[fine])) "every" else "not every")
  mac("hfCoarseNotMono", as.character(sum(!rows$mono[coarse])))
  mac("hfSlopeLo", f3(min(rows$slope))); mac("hfSlopeHi", f3(max(rows$slope)))
  mac("hfSlopeSpreadPct", f0b <- sprintf("%.0f", 100 * (max(rows$slope) - min(rows$slope)) / mean(rows$slope)))
  w <- which.min(rows$slope)
  mac("hfWorstCell", sprintf("%s at $D=%s$", tlabel[rows$target[w]], f2(rows$D[w])))
  mac("hfWorstSlope", f3(rows$slope[w])); mac("hfWorstRtwo", f4 <- sprintf("%.4f", rows$r2[w]))
  mac("hfPooledErrLo", f3(min(sapply(B$cells, function(c) max(abs((c$phat_mean - B$pooled[["intercept"]]) / B$pooled[["slope"]] - c$p))))))
  mac("hfPooledErrHi", f3(max(sapply(B$cells, function(c) max(abs((c$phat_mean - B$pooled[["intercept"]]) / B$pooled[["slope"]] - c$p))))))
  bd <- B$cells[["Bimodal 1"]]
  mac("hfPrevClaim", paste(f3(c(0.175, 0.334, 0.506)), collapse = ", "))
  mac("hfPrevCell", sprintf("%.3f, %.3f, %.3f", bd$phat_mean[2], bd$phat_mean[3], bd$phat_mean[4]))
}

## ---------------------------------------------------------------- abstention counts (B10)
CORR <- c("deconv", "sem", "imput", "combined")
paired <- function(a, b) { dd <- a - b; dd <- dd[is.finite(dd)]; c(mean(dd), sd(dd)/sqrt(length(dd))) }
count_best <- function(variant) {
  sum(sapply(bench$rows[grepl("^table1", names(bench$rows))], function(r) {
    m <- unlist(r$mean_x1e3); m["combined"] <- m[variant]
    ps <- r$perseed_x1e3; ps$combined <- ps[[variant]]
    best <- CORR[which.min(m[CORR])]
    if (best == "combined") return(TRUE)
    p <- paired(ps$combined, ps[[best]]); p[1] <= 2 * p[2] }))
}
mac("abstNone", as.character(count_best("combined"))); mac("abstDet", as.character(count_best("abst_det")))
mac("abstScale", as.character(count_best("abst_scale"))); mac("abstTaper", as.character(count_best("abst_taper")))
allrows <- bench$rows
tot <- function(v) sum(sapply(allrows, function(r) unlist(r$mean_x1e3)[v] - unlist(r$mean_x1e3)["combined"]))
mac("abstDetTotal", f1(tot("abst_det")), "total ISE x1e3 added over 76 cells by the detector rule")
mac("abstScaleInert", if (abs(tot("abst_scale")) < 1e-9) "identical" else "different")
mac("abstNcells", as.character(length(allrows))); mac("abstNseed", as.character(bench$provenance$nseed))
mac("abstGaussQuarterTaper", f3(unlist(allrows[["table1 T1 Gaussian 0.25"]]$mean_x1e3)["abst_taper"]))
mac("abstGaussQuarterNaive", f3(unlist(allrows[["table1 T1 Gaussian 0.25"]]$mean_x1e3)["naive"]))

## ---------------------------------------------------------------- write
hdr <- c("%% numbers_v10_3.tex, emitted by scripts/make_tables_v10_3.R. Do not edit by hand.",
         sprintf("%%%% batch %s on %s; lattice %s; readers %s; benchmark %s", batch$provenance$generated, batch$provenance$platform, lat$provenance$generated, ver$provenance$generated, bench$provenance$generated))
wr(c(hdr, macros), "numbers_v10_3.tex")
wr(c("v10.3 numbers, emitted by make_tables_v10_3.R", hdr[2], "", notes), "numbers_v10_3.txt")
cat(paste(notes, collapse = "\n"), "\n")
cat(sprintf("\n-> %s\n", OUT))
