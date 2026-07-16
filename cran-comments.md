## Submission

New submission of the package `adheaping` (version 1.0.0).

## Test environments

* local: macOS Sequoia 15.3.2 (aarch64-apple-darwin), R 4.6.0
* win-builder R-release: R version 4.6.1 (2026-06-24 ucrt) -- Status: 1 NOTE
* win-builder R-devel: R Under development (unstable) (2026-07-15 r90261 ucrt) -- Status: 1 NOTE
* GitHub Actions (.github/workflows/R-CMD-check.yaml): ubuntu (devel/release/oldrel), macOS, windows -- all passing

## R CMD check results

0 errors | 0 warnings | 1 note

The single NOTE (win-builder, both R-release and R-devel) is the expected
new-submission NOTE:

* "New submission" -- this is the package's first submission to CRAN.
* "Possibly misspelled words in DESCRIPTION: De, de, Heitjan, deconvolution,
  deconvolving" -- these are correctly spelled: "De"/"de" are parts of the
  hyphenated method name "de-heaping", "Heitjan" is the proper name of a cited
  author (Heitjan-Rubin multiple imputation), and "deconvolution"/"deconvolving"
  are standard statistical terms.

## Reverse dependencies

None (new package).

## Notes for the reviewer

* Suggested packages (Kernelheaping, foreign) are used conditionally via
  requireNamespace(); the package's own tests and examples do not require them.
