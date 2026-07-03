# Robust Kernel Density Estimation for Heaped and Rounded Data using Algebraic Diversity


Accompanying repository for a manuscript currently under review. Heaped and rounded data defeat kernel
density estimation: rounding aliases probability mass into the high frequencies, so frequency-domain
methods that assume the exact-sampling noise floor of 1/n place their cutoff wrongly and fail by
orders of magnitude. The methods build on the estimator machinery of arXiv:2606.15450 (the
AD-Wiener plug-in, the algebraic-residue floor, and the AD-KDE superposition) and contributes,
against that cited background:

1. The failure analysis: heaping is a floor-elevation phenomenon in the characteristic-function
   domain; the measured (residue) floor tracks it and stays robust under moderate pitch where the
   fixed 1/n floor fails by orders of magnitude.
2. A gated first-order peel for coarse pitch, where the rounding structure is deterministic and no
   white-floor model can describe it: blind pitch detection from the spectrum's exact replica
   periodicity (resurgence above the quietest valley, persistence across harmonics, an n-free
   core-to-floor band-capacity condition), capped deconvolution of the detected kernel with a
   colored-noise residue flow, and strict abstention: on exact samples the pipeline is unchanged,
   and every abstention falls back to the residue-floor estimator. Extended to mixed grain
   lattices: Moebius inversion over the cyclic subgroup lattice recovers grain weights and the
   unrounded-report share from replica center amplitudes, and the peel deconvolves the mixture
   kernel (`scripts/exp_heaping_mixed_v1.py`).
3. Ground-truthed validation on NHANES 2017-2018 survey measurements (controlled coarsening of
   measured adult weights; naturally heaped cigarette counts, qualitative).

## Contents
- `scripts/ad_kde_v31.py` - method module (corrected direct Botev fixed-point ISJ; `ad_kde_v30.py`
  is a compatibility shim kept for imports)
- `scripts/exp_heaping_coarse_v1.py` - the coarse-heaping battery (Table: coarse), gate of record
- `scripts/exp_heaping_mixed_v1.py` - the mixed-grain subgroup-lattice battery (Table: mixed)
- `scripts/exp_nhanes_heaping_v16.py` - the mixed-grain NHANES study (local run, same --data-dir as v15)
- `scripts/exp_nhanes_heaping_v15.py` - the NHANES study (Table: real; Fig. panels), needs the
  public NHANES 2017-2018 files DEMO_J.xpt, BMX_J.xpt, SMQ_J.xpt in --data-dir
- `scripts/exp_datagen_v30.py`, `scripts/adkde_plugins.py` - superposition machinery and plugin hooks
- `results/exp_heaping_coarse_v1.json` - battery results of record (seed 20260627)
- `DATA.md` - data provenance and reproduction commands; `HISTORY.md` - revision record

All synthetic results reproduce with numpy/scipy/matplotlib only; the NHANES study additionally
needs pandas. NHANES files are public-domain U.S. government data and are not redistributed here.

## Reproduction map

Every table and figure in the manuscript regenerates from a single script:

| Result | Script | Command |
|---|---|---|
| Coarse-heaping table (three targets x four pitches, gates and abstention) | `scripts/exp_heaping_coarse_v1.py` | `python3 exp_heaping_coarse_v1.py` |
| Mixed-grain table (subgroup-lattice peel, five mixtures) | `scripts/exp_heaping_mixed_v1.py` | `python3 exp_heaping_mixed_v1.py` |
| Measured-weight figure, both panels, and the imposed-grid study | `scripts/exp_nhanes_heaping_v15.py` | `python3 exp_nhanes_heaping_v15.py --data-dir data/nhanes` |
| Mixed-grain study on measured weights | `scripts/exp_nhanes_heaping_v16.py` | `python3 exp_nhanes_heaping_v16.py --data-dir data/nhanes` |
| Cigarette-count lattice decomposition | `scripts/exp_cig_lattice_probe_v1.py` | `python3 exp_cig_lattice_probe_v1.py --data-dir data/nhanes` |

The NHANES studies need the public 2017-2018 files `DEMO_J.xpt`, `BMX_J.xpt`, and `SMQ_J.xpt`
in the `--data-dir`; all other results are fully synthetic with fixed recorded seeds and run
with no external data.
