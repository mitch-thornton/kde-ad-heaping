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

* local: macOS Sequoia (aarch64-apple-darwin), R 4.6.1
* win-builder R-release
* win-builder R-devel
* GitHub Actions (.github/workflows/R-CMD-check.yaml): ubuntu (devel/release/oldrel), macOS, windows

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

None.

## Notes for the reviewer

* Suggested packages (Kernelheaping, foreign) are used conditionally via
  requireNamespace(); the package's own tests and examples do not require them.
* The new `heap_lattice()` gains three regression tests that assert its documented regime,
  including one that asserts the grain split is NOT accurate below that regime, so the
  stated limitation is itself tested.
