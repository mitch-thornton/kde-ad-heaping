#!/usr/bin/env Rscript
## make_tables_v10_2.R
## Emits every table in the v10 revision directly from benchmark_v10_2.rds, so no
## experimental value is typed by hand and a rerun puts the rerun's numbers into the PDF.
##
## Outputs, written to OUTDIR:
##   table1_v10.tex         main-text Table 1, four densities, 50 seeds
##   supp_tableS2_mw15.tex  supplementary, all fifteen Marron-Wand densities
##   supp_tableS3_paired.tex supplementary, paired differences against the best corrected
##                           method with their paired standard errors, so every tie
##                           decision in Table 1 and Table S2 can be checked
##   tables_v10_summary.txt  counts and provenance, for the response letter
##
## You use zsh. Invoke through the wrapper as
##     bash run_make_tables_v10_2.sh
## or directly as
##     RESULTS=/Users/mitch/src/KDE-AD-HEAPING/kde-ad-heaping/results Rscript make_tables_v10_2.R
##
## Tie rule: a corrected method is marked as tied with the best corrected method when the
## mean of the per-seed paired differences is at most two paired standard errors above
## zero. This is the rule the v9.4 caption stated. The paired standard errors are now
## emitted in Table S3 so the decision is checkable from the paper, which is what
## Reviewer 2 asked for in major point 7.

RES <- Sys.getenv("RESULTS", unset = ".")
OUT <- Sys.getenv("OUTDIR", unset = RES)
rdsfile <- file.path(RES, "benchmark_v10_2.rds")
if (!file.exists(rdsfile)) stop("benchmark_v10_2.rds not found under RESULTS=", RES)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

obj  <- readRDS(rdsfile)
prov <- obj$provenance
rows <- obj$rows

CORR <- c("deconv", "sem", "imput", "combined")
LBL  <- c(naive = "naive", deconv = "deconv", sem = "SEM", imput = "imput",
          deheap = "de-heap", super = "super", combined = "combined")

paired <- function(a, b) {
  dd <- a - b; dd <- dd[is.finite(dd)]
  if (length(dd) < 2) return(c(mean = NA_real_, se = NA_real_))
  c(mean = mean(dd), se = stats::sd(dd)/sqrt(length(dd)))
}
tieset <- function(r) {
  m <- unlist(r$mean_x1e3); av <- CORR[!is.na(m[CORR])]
  best <- av[which.min(m[av])]
  tie <- best
  for (cc in setdiff(av, best)) {
    p <- paired(r$perseed_x1e3[[cc]], r$perseed_x1e3[[best]])
    if (!is.na(p["se"]) && p["mean"] <= 2 * p["se"]) tie <- c(tie, cc)
  }
  list(best = best, tie = tie)
}
fmt <- function(v, sd, bold) {
  if (is.na(v)) return("--")
  s <- if (abs(v) >= 100) sprintf("%.1f", v) else sprintf("%.3f", v)
  d <- if (is.na(sd)) "" else sprintf("\\,(%s)", if (abs(sd) >= 100) sprintf("%.1f", sd) else sprintf("%.2f", sd))
  if (bold) sprintf("\\textbf{%s%s}", s, d) else paste0(s, d)
}
esc <- function(s) gsub("([&%#_])", "\\\\\\1", s)
## display names: the run labels the fourth Table 1 target explicitly as a paper variant;
## the table itself is narrower with the short form and the caption already says the
## bimodal and strongly skewed targets are custom rather than Marron-Wand members.
shortname <- function(s) sub(" \\(paper variant\\)$", "", s)

sel <- function(set) {
  k <- names(rows)[sapply(rows, function(r) r$set == set)]
  rows[k]
}
ordkey <- function(rs) {
  o <- order(sapply(rs, function(r) suppressWarnings(as.numeric(r$mw_member))),
             sapply(rs, function(r) r$D), na.last = TRUE)
  rs[o]
}

## ---------------------------------------------------------------- main Table 1
t1 <- sel("table1")
t1 <- t1[order(sapply(t1, function(r) r$index), sapply(t1, function(r) r$D))]
nwin <- 0
lines <- c()
prev <- ""
for (r in t1) {
  ts <- tieset(r); m <- unlist(r$mean_x1e3); s <- unlist(r$sd_x1e3)
  if ("combined" %in% ts$tie) nwin <- nwin + 1
  nm <- if (r$name == prev) "" else esc(shortname(r$name)); prev <- r$name
  if (nm != "" && length(lines)) lines <- c(lines, "\\midrule")
  cells <- c(fmt(m["naive"], s["naive"], FALSE),
             sapply(c("deconv","sem","imput"), function(c) fmt(m[c], s[c], c %in% ts$tie)),
             fmt(m["deheap"], s["deheap"], FALSE), fmt(m["super"], s["super"], FALSE),
             fmt(m["combined"], s["combined"], "combined" %in% ts$tie))
  lines <- c(lines, sprintf("%s & %.2f & %s\\\\", nm, r$D, paste(cells, collapse = " & ")))
}
cap1 <- sprintf(paste0(
  "Integrated squared error ($\\times10^3$) against the true density as the heaping grid $D$ coarsens. ",
  "Mean over %d seeds with the seed-to-seed standard deviation in parentheses, $n=%d$, computed in the ",
  "accompanying R library on a grid of $M=%d$ cells of width $\\Delta x=%.6f$. ",
  "The four targets are a selected set. Two are members 1 and 4 of the normal-mixture benchmark of ",
  "Marron and Wand; the bimodal and strongly skewed targets are custom and are defined in ",
  "Supplementary Table~S1. All fifteen members of that benchmark are evaluated in Supplementary Table~S2. ",
  "The three external methods are deconv, SEM and imput. The proposed method is the combined estimator, ",
  "whose two constituent components are set off by vertical rules. In each row the best of the four ",
  "heaping-corrected values is bold, together with any corrected value whose per-seed paired difference ",
  "from the best is within two paired standard errors; those paired differences and their standard errors ",
  "are reported in Supplementary Table~S3. The naive kernel is the uncorrected reference. ",
  "The combined estimator is the best or tied with the best in %d of the 16 cells."),
  prov$nseed, prov$n, prov$M, prov$dx, nwin)
writeLines(c("\\begin{table*}[t]", "\\centering", paste0("\\caption{", cap1, "}"),
  "\\label{tab:bench}", "\\footnotesize",
  "\\begin{tabular}{ll rrrr | rr | r}", "\\toprule",
  "density & $D$ & naive & deconv & SEM & imput & de-heap & super & combined\\\\",
  "\\midrule", lines, "\\bottomrule", "\\end{tabular}", "\\end{table*}"),
  file.path(OUT, "table1_v10.tex"))

## ------------------------------------------------- supplementary Table S2, MW15
mw <- ordkey(sel("mw15"))
nwin2 <- 0; lines2 <- c(); prev <- ""
for (r in mw) {
  ts <- tieset(r); m <- unlist(r$mean_x1e3); s <- unlist(r$sd_x1e3)
  if ("combined" %in% ts$tie) nwin2 <- nwin2 + 1
  nm <- if (r$name == prev) "" else sprintf("%s.~%s", r$mw_member, esc(shortname(r$name))); prev <- r$name
  if (nm != "" && length(lines2)) lines2 <- c(lines2, "\\midrule")
  cells <- c(fmt(m["naive"], s["naive"], FALSE),
             sapply(c("deconv","sem","imput"), function(c) fmt(m[c], s[c], c %in% ts$tie)),
             fmt(m["combined"], s["combined"], "combined" %in% ts$tie))
  lines2 <- c(lines2, sprintf("%s & %.2f & %s\\\\", nm, r$D, paste(cells, collapse = " & ")))
}
cap2 <- sprintf(paste0(
  "The comparison of Table~1 on all fifteen normal-mixture densities of ",
  "Marron and Wand. Integrated squared error ($\\times10^3$), mean over %d seeds with the seed-to-seed ",
  "standard deviation in parentheses, $n=%d$, same grid, same seeds and same estimators as Table~1. ",
  "Members 1 and 4 also appear in Table~1 and the values agree exactly. ",
  "Bolding follows Table~1. The combined estimator is the best or tied with the best in %d of the 60 cells."),
  prov$nseed, prov$n, nwin2)
writeLines(c("\\begin{table*}[t]", "\\centering", paste0("\\caption{", cap2, "}"),
  "\\label{tab:mw15}", "\\footnotesize",
  "\\begin{tabular}{ll rrrr r}", "\\toprule",
  "density & $D$ & naive & deconv & SEM & imput & combined\\\\",
  "\\midrule", lines2, "\\bottomrule", "\\end{tabular}", "\\end{table*}"),
  file.path(OUT, "supp_tableS2_mw15.tex"))

## ------------------------------- supplementary Table S3, paired differences
lines3 <- c(); prev <- ""
for (r in c(t1, mw)) {
  ts <- tieset(r)
  nm <- sprintf("%s %s", r$set, esc(shortname(r$name)))
  if (nm == prev) nm2 <- "" else nm2 <- nm; prev <- nm
  ds <- sapply(setdiff(CORR, ts$best), function(cc) {
    p <- paired(r$perseed_x1e3[[cc]], r$perseed_x1e3[[ts$best]])
    if (is.na(p["mean"])) "--" else sprintf("%.3f\\,(%.3f)", p["mean"], p["se"])
  })
  lines3 <- c(lines3, sprintf("%s & %.2f & %s & %s\\\\", nm2, r$D, LBL[ts$best],
                              paste(ds, collapse = " & ")))
}
cap3 <- paste0("Per-seed paired differences against the best corrected ",
  "method in each cell, mean with the paired standard error in parentheses, in the units of ",
  "Table~1. A value at or below twice its paired standard error is a tie and is bolded in ",
  "Table~1 and Supplementary Table~S2. Columns are the three corrected methods other than the ",
  "best one, in the fixed order deconv, SEM, imput, combined, with the best one omitted.")
writeLines(c("\\begin{table*}[t]", "\\centering", paste0("\\caption{", cap3, "}"),
  "\\label{tab:paired}", "\\scriptsize", "\\begin{tabular}{ll l rrr}", "\\toprule",
  "cell & $D$ & best & \\multicolumn{3}{c}{others, mean paired difference (paired SE)}\\\\",
  "\\midrule", lines3, "\\bottomrule", "\\end{tabular}", "\\end{table*}"),
  file.path(OUT, "supp_tableS3_paired.tex"))

## ------------------------------------------------------------------- summary
byD <- sapply(c(0.25, 0.5, 1.0, 1.5), function(DD) {
  cs <- Filter(function(r) r$set == "mw15" && r$D == DD, mw)
  sum(sapply(cs, function(r) "combined" %in% tieset(r)$tie))
})
sm <- c(
  "v10 table generation summary",
  sprintf("generated from %s", normalizePath(rdsfile)),
  sprintf("benchmark run %s on %s, %s", prov$generated, prov$platform, prov$r_version),
  sprintf("elapsed %.1f s, seeds %d, n %d, M %d, dx %.9f, grids %s",
          prov$elapsed_sec, prov$nseed, prov$n, prov$M, prov$dx, prov$Ds),
  sprintf("seed rule: %s", prov$seed_rule),
  sprintf("Kernelheaping available: %s", prov$kernelheaping),
  "",
  sprintf("Table 1: combined best or tied in %d of 16 cells", nwin),
  sprintf("Supplementary Table S2: combined best or tied in %d of 60 cells", nwin2),
  sprintf("  by grid: D=0.25 %d/15, D=0.50 %d/15, D=1.00 %d/15, D=1.50 %d/15",
          byD[1], byD[2], byD[3], byD[4]),
  sprintf("Total cells reported: %d", length(rows)))
writeLines(sm, file.path(OUT, "tables_v10_summary.txt"))
cat(paste(sm, collapse = "\n"), "\n\nwrote table1_v10.tex, supp_tableS2_mw15.tex, supp_tableS3_paired.tex, tables_v10_summary.txt to\n  ",
    normalizePath(OUT), "\n")
