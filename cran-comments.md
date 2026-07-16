## Submission

New submission of the package `adheaping` (version 1.0.0).

## Test environments

* local: macOS Sequoia 15.3.2 (aarch64-apple-darwin), R 4.6.0
* win-builder: R-devel and R-release        [run devtools::check_win_devel()/check_win_release(); paste results before submitting]
* R-hub v2 (GitHub Actions): linux, macos, windows   [run rhub::rhub_check(); confirm green]
* GitHub Actions (.github/workflows/R-CMD-check.yaml): ubuntu (devel/release/oldrel), macOS, windows

## R CMD check results

0 errors | 0 warnings | 1 note

* Note: "New submission." This is the package's first submission to CRAN.
* The package URL (https://github.com/mitch-thornton/kde-ad-heaping) resolves once the
  repository is public; it returned 404 only while the repository was private.

## Reverse dependencies

None (new package).

## Notes for the reviewer

* Suggested packages (Kernelheaping, foreign) are used conditionally via requireNamespace();
  the package's own tests and examples do not require them.
