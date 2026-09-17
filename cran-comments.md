## Submission

Update of `adheaping` from 1.0.0 to 1.1.0.

This release corrects one exported function. `heap_lattice()` is reimplemented to perform
the Moebius inversion its documentation always described. The 1.0.0 implementation
normalized replica amplitudes without inverting, and the `unrounded` share it returned was
identically its grain-one weight, so two quantities it reported separately were one number.
The new implementation was validated against known truth over a sweep of base scales, and
its regime of validity is now stated in the documentation.

The argument list of `heap_lattice()` changes from `(y, grains, span, Mgrid)` to
`(y, grains, nboot, seed)`. The package has no reverse dependencies, and the function is
the only one whose behavior changes. Every density estimator is untouched, so estimates
from 1.0.0 are unchanged.

Several documentation statements were also corrected against experiments run for the
accompanying paper, where a blind grid reading and an automatic comb detection turned out
to be weaker than the 1.0.0 documentation implied.

## Test environments

* local: macOS Sequoia 15.3.2 (aarch64-apple-darwin24.6.0), R 4.6.1, `R CMD check --as-cran`
* GitHub Actions, `.github/workflows/R-CMD-check.yaml`: macOS release, Windows release,
  Ubuntu R-devel, Ubuntu release, Ubuntu oldrel-1

## R CMD check results

0 errors | 0 warnings | 0 notes

The local run reports one further note, that HTML validation was skipped because the
installed HTML Tidy is older than R requires. That is a property of the local toolchain
rather than of the package.

## Reverse dependencies

None.

## Notes for the reviewer

* Suggested packages (Kernelheaping, foreign) are used conditionally via
  requireNamespace(); the package's own tests and examples do not require them.
* The new `heap_lattice()` gains three regression tests that assert its documented regime,
  including one that asserts the grain split is NOT accurate below that regime, so the
  stated limitation is itself tested.
