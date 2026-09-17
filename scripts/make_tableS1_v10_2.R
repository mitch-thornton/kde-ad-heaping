#!/usr/bin/env Rscript
## make_tableS1_v10_2.R -- emits Supplementary Table S1, the definitions of every
## density used in the paper, from mw15.R and the two custom Table 1 targets, so the
## definitions in the supplement cannot drift from the ones the benchmark runs.
## Usage:  MW15_R=/path/to/mw15.R OUTDIR=/path/to/tables Rscript make_tableS1_v10_2.R
.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) dirname(normalizePath(f)) else "." }
source(Sys.getenv("MW15_R", unset = file.path(.sd(), "mw15.R")))
OUT <- Sys.getenv("OUTDIR", unset = ".")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fr <- function(v) paste(sapply(v, function(x) {
  r <- MASS::fractions(x); if (requireNamespace("MASS", quietly = TRUE)) as.character(r) else sprintf("%.4f", x)
}), collapse = ", ")
num <- function(v) paste(sprintf("%.4f", v), collapse = ", ")
row <- function(tag, nm, d) sprintf("%s & %s & %d & %s & %s & %s\\\\",
  tag, gsub("([&%#_])", "\\\\\\1", nm), length(d$w), num(d$w), num(d$mu), num(d$sd))
CUSTOM <- list(
  list(tag = "C1", nm = "Bimodal (custom, Table 1)",
       w = c(.5,.5), mu = c(-1.2,1.2), sd = c(.5,.5)),
  list(tag = "C2", nm = "Strongly skewed (custom, Table 1)",
       w = rep(.2,5), mu = c(0,.5,1.0833,1.4167,1.6875), sd = c(1,2/3,4/9,8/27,16/81)))
lines <- c(sapply(MW15, function(d) row(sprintf("MW%d", d$index), d$name, d)),
           "\\midrule",
           sapply(CUSTOM, function(d) row(d$tag, d$nm, d)))
cap <- paste0("Every density used in this work, as a normal mixture ",
  "$\\sum_j w_j\\,\\mathcal{N}(\\mu_j,\\sigma_j^2)$. Rows MW1 to MW15 are the fifteen members of ",
  "the benchmark of Marron and Wand and are evaluated in Supplementary Table~S2. Rows C1 and C2 ",
  "are the two custom targets used in Table~1 of the main text; they are not members of that ",
  "benchmark. The other two Table~1 targets are MW1 and MW4 exactly. These definitions are ",
  "emitted from the same source file the benchmark uses, so they cannot drift from it.")
writeLines(c("\\begin{table*}[t]", "\\centering", paste0("\\caption{", cap, "}"),
  "\\label{tab:defs}", "\\scriptsize", "\\begin{tabular}{ll r p{4.2cm} p{4.2cm} p{4.0cm}}", "\\toprule",
  "tag & density & $J$ & weights $w_j$ & means $\\mu_j$ & standard deviations $\\sigma_j$\\\\",
  "\\midrule", lines, "\\bottomrule", "\\end{tabular}", "\\end{table*}"),
  file.path(OUT, "supp_tableS1_defs.tex"))
cat("wrote supp_tableS1_defs.tex to", normalizePath(OUT), "\n")
