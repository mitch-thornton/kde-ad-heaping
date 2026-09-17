# adheaping 1.1.0

* `heap_lattice()` is reimplemented. It now performs the Moebius inversion the
  documentation always described, reading the empirical characteristic function directly
  at the exclusive replica centers of each grain with no binning and with the sampling
  floor removed. Its arguments change from `span` and `Mgrid` to `nboot` and `seed`, and
  it gains bootstrap standard errors.

  The 1.0.0 reader did not invert. It normalized replica amplitudes, and the `unrounded`
  share it returned was identically its grain-one weight, so the two quantities it
  reported separately were one number. On a known-truth simulation it was biased on every
  grain. Any weights obtained from it are superseded.

  The new reader has a stated regime. It recovers every weight to within about four
  percent once kappa, the standard deviation of the dequantized base over the coarsest
  grain, exceeds roughly three, and recovers the unit share at every kappa tested. Below
  that the grain split is descriptive rather than a measurement, and the documentation
  says so.

* Documentation corrections, from experiments run for the accompanying paper. `heap_grid()`
  now records that its blind mode returned the imposed grid in 0 of 12 test cases where its
  verifying mode returned it in 12 of 12, so a blind reading is a candidate rather than a
  recovered grid. The README no longer describes `heap_detect()` as blind spectral
  detection; with its shipped settings it abstains in most cells. The package Description
  no longer calls the grid reader blind, and no longer calls the base-R baselines faithful,
  since the Heitjan-Rubin replica has not been checked against a published implementation.

* No change to `deheap_kde()`, `superpose_kde()`, `adkde()`, `heap_fraction()`,
  `heap_detect()`, or any baseline. Density estimates from 1.0.0 are unchanged.

# adheaping 1.0.0

* Initial release. Characteristic-function de-heaping density estimation for heaped and
  rounded data:
  - tuning-free box-deconvolution de-heaping estimator (`deheap_kde`), superposition
    (`superpose_kde`), and the combined band-capacity-gated estimator (`adkde`);
  - blind grid, heaped-fraction, and mixed-grain (subgroup-lattice) readers
    (`heap_grid`, `heap_fraction`, `heap_lattice`);
  - a spectral higher-order comb detector (`heap_detect`);
  - faithful base-R replicas of measurement-error deconvolution (`deconv_kde`) and the
    Heitjan-Rubin multiple-imputation approach (`heitjan_mi`), and a wrapper for the
    Kernelheaping stochastic-EM estimator (`sem_kde`).
  Derived from the spectral-decomposition kernel density estimation of Thornton
  (arXiv:2606.15450).
* Release archived on Zenodo: doi:10.5281/zenodo.21181205.
