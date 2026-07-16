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
