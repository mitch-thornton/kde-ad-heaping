#!/usr/bin/env Rscript
## make_tableS12_v13.R
## Emits Supplementary Table S12 and the macros the main text and the supplement quote for
## the sample-size sensitivity analysis requested in round five. Nothing here is typed: every
## value is read out of results/sensitivity_n_v13.rds, which carries its own provenance block.
##
## Input, under RESULTS:
##   sensitivity_n_v13.rds   the sample-size sweep, written by scripts/sensitivity_n_v13.R
##   benchmark_v10_2.rds     the 50-seed benchmark, used only to re-check the anchor row here
##
## Outputs, written to OUTDIR:
##   supp_tableS12_n.tex     Table S12, the bimodal rows
##   numbers_v13.tex         \newcommand macros for both documents
##   numbers_v13.txt         the same numbers in plain text, for the response letter
##
## TARGET selects which density the table shows, by its Table 1 index. The run computes both
## T1 (Gaussian) and T2 (Bimodal); the editor asked for one representative scenario and the
## bimodal is the target the grid sweep of Fig. 2 uses, so T2 is the default. Setting
## TARGET=T1,T2 emits both, which has consequences discussed in the v13.0 notes and is not
## the shipped choice.
##
## zsh users: invoke as   bash run_make_tableS12_v13.sh

RES    <- Sys.getenv("RESULTS", unset = ".")
OUT    <- Sys.getenv("OUTDIR", unset = RES)
TARGET <- strsplit(Sys.getenv("TARGET", "T2"), ",")[[1]]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

need <- function(f) {
  p <- file.path(RES, f)
  if (!file.exists(p)) stop(f, " not found under RESULTS=", RES)
  readRDS(p)
}
S <- need("sensitivity_n_v13.rds")
P <- S$provenance

## ---- refuse to emit from a run that did not verify itself ------------------------------
## The run checks its own n = 4000 row against the published Table 1 cell. If that check did
## not run, or did not pass, the numbers must not reach a document.
an <- P$table1_anchor
if (!isTRUE(an$ran))
  stop("the stored run did not check itself against Table 1 (", an$note,
       "), so nothing is emitted from it")
if (!isTRUE(an$pass))
  stop("the stored run failed its own Table 1 anchor, so nothing is emitted from it")
if (isTRUE(P$smoke))
  stop("the stored run is a smoke run, so nothing is emitted from it")

## ---- and re-check it here, independently of what the run recorded ----------------------
## The run's own verdict is a claim in a file. This repeats the comparison from the two
## result files sitting side by side, so the table cannot be built on a stale assertion.
bench <- need("benchmark_v10_2.rds")
METH  <- c("naive", "deconv", "combined")
recheck <- list()
for (ti in unique(c(TARGET, "T2"))) {
  rk <- grep(sprintf("^%s ", ti), names(S$rows), value = TRUE)
  if (!length(rk)) stop("no rows for target ", ti, " in the stored run")
  nm <- S$rows[[rk[1]]]$name
  me <- S$rows[[sprintf("%s %s D%.2f n%d", ti, nm, P$D, P$n_draw)]]
  if (is.null(me)) stop("the stored run has no n = ", P$n_draw, " row for ", ti)
  bk <- sprintf("table1 %s %s %s", ti, nm, P$D)
  br <- bench$rows[[bk]]
  if (is.null(br)) stop("the stored benchmark has no row keyed '", bk, "'")
  d <- sapply(METH, function(k) abs(me$mean_x1e3[[k]] - br$mean_x1e3[[k]]))
  if (any(d > 5e-4))
    stop("the n = ", P$n_draw, " means for ", ti, " do not match the published Table 1 cell: ",
         paste(sprintf("%s %.6f", METH, d), collapse = ", "))
  recheck[[ti]] <- d
}

## ---- the rows for the selected target, in increasing n --------------------------------
pick <- function(ti) {
  ks <- grep(sprintf("^%s ", ti), names(S$rows), value = TRUE)
  rs <- S$rows[ks]
  rs[order(sapply(rs, function(r) r$n))]
}
sel <- unlist(lapply(TARGET, pick), recursive = FALSE)
if (!length(sel)) stop("TARGET=", paste(TARGET, collapse = ","), " selected no rows")

f3 <- function(x) sprintf("%.3f", x)
f4 <- function(x) sprintf("%.4f", x)

## ---- Table S12 -------------------------------------------------------------------------
multi <- length(TARGET) > 1L
hdr <- paste0(if (multi) "target & " else "",
              "$n$ & naive & deconv & combined & combined/naive & naive $-$ combined & SE\\\\")
body <- sapply(sel, function(r) sprintf("%s%d & %s & %s & %s & %s & %s & %s\\\\",
  if (multi) paste0(r$name, " & ") else "",
  r$n,
  f3(r$mean_x1e3[["naive"]]), f3(r$mean_x1e3[["deconv"]]), f3(r$mean_x1e3[["combined"]]),
  f3(r$ratio_combined_naive),
  f3(r$advantage_vs_naive_x1e3[["mean"]]), f3(r$advantage_vs_naive_x1e3[["se"]])))

ns     <- sapply(sel, function(r) r$n)
rat    <- sapply(sel, function(r) r$ratio_combined_naive)
wins   <- sapply(sel, function(r) r$wins_vs_naive)
picks  <- sapply(sel, function(r) paste(names(r$pick_counts), collapse = "+"))
nseed  <- sel[[1]]$nseed
allwin <- all(wins == nseed)
onepick <- length(unique(picks)) == 1L && identical(unique(picks), "deheap")
tgtname <- paste(unique(sapply(sel, function(r) r$name)), collapse = " and ")

cap <- sprintf(paste0(
  "Sample-size sensitivity at a fixed grid. Mean integrated squared error ",
  "($\\times10^{3}$) for the uncorrected kernel, the deconvolution baseline and the ",
  "combined estimator on the %s target at $D=%.2f$, over %d seeds, at five sample sizes. ",
  "The draws are nested: each seed draws the $n=%d$ sample of Table~1 under the same seed ",
  "rule and the smaller samples are the leading $n$ of that draw, so the comparison across ",
  "$n$ is paired and the last row reproduces the corresponding cell of Table~1 exactly. ",
  "The sixth column is the ratio of the combined estimator to leaving the data alone, and ",
  "the last two are the paired per-seed difference and its standard error. %s%s"),
  tolower(tgtname), P$D, nseed, P$n_draw,
  if (allwin) sprintf("The combined estimator is lower than the uncorrected estimate in all %d seeds at every sample size. ", nseed) else "",
  if (onepick) sprintf("The band-capacity statistic selected the de-heaping component in all %d seeds at every sample size, so the gate does not change its choice as the sample shrinks.", nseed) else "")

tex <- c("\\begin{table}[t]", "\\centering",
         sprintf("\\caption{%s}", cap),
         "\\label{tab:s12}",
         sprintf("\\begin{tabular}{%s}", paste0(if (multi) "l" else "", "rrrrrrr")),
         "\\toprule", hdr, "\\midrule", body, "\\bottomrule",
         "\\end{tabular}", "\\end{table}")

## ---- macros ----------------------------------------------------------------------------
macros <- character(0); notes <- character(0)
mac <- function(name, value, note = NULL) {
  macros <<- c(macros, sprintf("\\newcommand{\\%s}{%s}", name, value))
  notes  <<- c(notes, sprintf("%-22s %s%s", name, value, if (is.null(note)) "" else paste0("   ", note)))
}

lo <- which.min(ns); hi <- which.max(ns)
mac("snSeeds",     as.character(nseed))
mac("snD",         sprintf("%.1f", P$D))
mac("snNlevels",   as.character(length(unique(ns))))
mac("snNmin",      format(min(ns), big.mark = ",", trim = TRUE))
mac("snNmax",      format(max(ns), big.mark = ",", trim = TRUE))
mac("snTarget",    tolower(tgtname))
mac("snRatioAtMax", f3(rat[hi]), "combined / naive at the largest sample")
mac("snRatioAtMin", f3(rat[lo]), "combined / naive at the smallest sample")
mac("snRatioWiden", sprintf("%.1f", rat[lo] / rat[hi]),
    "factor by which the ratio widens from the largest to the smallest sample")
mac("snMonotone", if (all(diff(rat[order(ns)]) < 0)) "decreases monotonically with $n$" else
                  "does not decrease monotonically with $n$",
    "direction of the ratio against n")
mac("snWinsAll", as.character(nseed))
mac("snWinsEvery", if (allwin) "every" else "not every")
mac("snGainAtMin", f3(sapply(sel, function(r) r$advantage_vs_naive_x1e3[["mean"]])[lo]))
mac("snGainAtMinSE", f3(sapply(sel, function(r) r$advantage_vs_naive_x1e3[["se"]])[lo]))
mac("snGainSEmin", sprintf("%.0f", min(sapply(sel, function(r)
      r$advantage_vs_naive_x1e3[["mean"]] / r$advantage_vs_naive_x1e3[["se"]]))),
    "smallest ratio of paired gain to its standard error across the sweep")
dvr <- sapply(sel, function(r) r$ratio_deconv_naive)
mac("snDeconvRatioAtMax", f3(dvr[hi]))
mac("snDeconvRatioAtMin", f3(dvr[lo]))
mac("snPickAll", if (onepick) sprintf("de-heaping component in all %d seeds at every sample size", nseed)
                 else "more than one component across the sweep")
mac("snAnchorMaxDiff", format(max(unlist(recheck)), scientific = TRUE, digits = 2),
    "largest absolute difference between the n = 4000 means and the published Table 1 cell")

## ---- write -----------------------------------------------------------------------------
wr <- function(lines, f) { writeLines(lines, file.path(OUT, f)); cat(sprintf("wrote %s\n", f)) }
stamp <- c(sprintf("%%%% emitted by scripts/make_tableS12_v13.R from %s", "results/sensitivity_n_v13.rds"),
           sprintf("%%%% run %s on %s, %s; %s seeds; anchor passed", P$generated, P$platform,
                   P$r_version, P$nseed),
           "%% Do not edit by hand.")
wr(c(stamp, tex), "supp_tableS12_n.tex")
wr(c(sub("^%% emitted by", "%% numbers_v13.tex, emitted by", stamp[1]), stamp[-1], macros),
   "numbers_v13.tex")
wr(c("v13 sample-size numbers, emitted by make_tableS12_v13.R",
     sub("^%%%% ", "", stamp[2]), "", notes), "numbers_v13.txt")

cat("\n", paste(notes, collapse = "\n"), "\n", sep = "")
cat(sprintf("\ntarget(s) emitted: %s   rows: %d\n", paste(TARGET, collapse = ","), length(sel)))
cat(sprintf("anchor re-checked here, largest mean difference %s (x1e3)\n",
            format(max(unlist(recheck)), scientific = TRUE, digits = 2)))
cat(sprintf("-> %s\n", OUT))
