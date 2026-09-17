# DATA.md

What every number in this paper was computed from, and how to obtain it. Three datasets
appear. Two are public survey files that are not redistributed here because their providers
ask that they be taken from the source. The third is generated, and its generator and seeds
ship with the bundle.

Nothing in this file is a claim about the data. It is a set of instructions for putting the
same inputs on disk that produced the tables.

## 1. NHANES 2017 to 2018, measured body weight

Used for Table 2 and Figure 4, the controlled coarsening study.

| | |
|---|---|
| Files | `DEMO_J.xpt`, `BMX_J.xpt` |
| Source | National Center for Health Statistics, NHANES 2017-2018 cycle |
| Landing page | https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Cycle=2017-2018 |
| Documentation | `DEMO_J.htm` and `BMX_J.htm` beside the data files |
| Place them in | any directory, then pass it as `NHANES_DIR` |

Both files sit in the 2017-2018 public data directory alongside their documentation pages.
Take them from the landing page above rather than from a deep link, since the file paths on
that site have changed at least once.

Construction, exactly as `scripts/make_figures_v10_7.R` part 3 performs it. `DEMO_J` supplies
`SEQN` and `RIDAGEYR`, `BMX_J` supplies `SEQN` and `BMXWT`, and the two are merged on `SEQN`.
Respondents aged twenty and over with a non-missing weight are retained, which gives
n = 5185. Kilograms are converted to pounds at 2.2046226218. The evaluation grid is M = 4096
cells spanning the measured range extended by forty pounds at each end, so the cell width is
0.1326 lb. The no-heaping reference is a Gaussian kernel estimate of those unrounded weights
at the Silverman bandwidth, 9.681 lb. Grids of D in {5, 10, 20, 30, 40} pounds are imposed by
rounding, and every method is estimated from the rounded values and scored against that
reference.

Two of the five columns draw from the random number generator, the SEM baseline through
`Kernelheaping` and the imputation baseline through its eight draws, so part 3 fixes the seed
of record, 20260627, before it runs. The seed is recorded in the provenance block of
`results/nhanes_v10_7.rds`.

## 2. NHANES 2017 to 2018, self-reported cigarettes per day

Used for the mixed-grain reading in the Results and Supplementary Table S8.

| | |
|---|---|
| File | `SMQ_J.xpt` |
| Source | the same 2017-2018 cycle and landing page as above |
| Place it in | the same `NHANES_DIR` |

Construction, as `scripts/test_lattice_v3.R` part 3 performs it. The first available column
among `SMD650`, `SMD641` and `SMQ020` is taken, and values are kept when they are finite,
strictly positive and below 200, which gives n = 1019 adult smokers reporting a daily count.
The reader is anchored on the unit grid and evaluated on grains 1, 5, 10 and 20. Bootstrap
standard errors come from 1000 resamples under the seed of record.

The grain weights reported in the previous version of this work came from a Python
implementation and are superseded. That history is recorded in Methods. The Python script
that produced them is not shipped, because three of its imports were never public and it
cannot run from the deposit.

## 3. Berlin resident register, December 2015

Used for Table 3 and Figure 5, the interval-grouping study.

| | |
|---|---|
| File | `EWR201512E_Matrix.csv` |
| Source | Amt fuer Statistik Berlin-Brandenburg, Einwohnerregisterstatistik, open data |
| Open data directory | https://www.statistik-berlin-brandenburg.de/opendata/ |
| Landing page | https://www.statistik-berlin-brandenburg.de/meine-region/berlin-statistik/einwohnerbestand/ |
| Place it in | a directory passed as the first argument to `berlin_R.R` |

The open data files follow the pattern `EWR<YYYYMM><letter>_Matrix.csv`, one matrix per
reporting date, so the December 2015 file is `EWR201512E_Matrix.csv`. Read it with
`read.csv2`, since it is semicolon separated with comma decimals and dotted thousands.

Construction, as `berlin_R.R` performs it. Age-band columns are the ones matching
`^E_E[0-9][0-9]_[0-9]+$`, whose name encodes the band edges. Counts are pooled over all 447
planning areas, covering 3,610,156 residents in bands of one to fifteen years. A sample of
40,000 values is then drawn from the pooled counts in proportion to their totals, each drawn
value placed uniformly at random within its native band, under the seed of record 20260627.
That dequantized sample is the native-resolution sample, and its Gaussian kernel estimate at
the Silverman bandwidth is the reference density. The evaluation grid is M = 2048 cells
spanning [-25, 135) years. Coarser uniform grains of D in {5, 10, 15, 20} years are imposed
by regrouping the same pooled counts.

Placing values uniformly within a cell is the operation the imputation baselines themselves
perform, so this construction favors those baselines at fine grains rather than the proposed
method. The paper states this where the study is reported.

## 4. Synthetic targets

Used for Table 1, Supplementary Tables S1 to S6 and S9 to S11, and Figures 1, 2 and 3.

No download. The targets are normal mixtures defined in `scripts/mw15.R`, which encodes the
fifteen densities of Marron and Wand (1992), Table 1, together with two custom targets, a
bimodal and a strongly skewed one, that are not members of that benchmark and are defined in
Supplementary Table S1. `scripts/verify_mw15.R` checks the definitions against the published
parameters.

Samples are drawn at n = 4000 on a grid of M = 2048 cells spanning [-10, 10), so the cell
width is 0.009766. Seeds are not listed. They are generated by the rule of record

    seed = 20260627 + 1000 * replicate + round(100 * D)

so every draw in the paper is reproducible from that one line. The heaped-fraction experiment
extends the rule with a term in the fraction, stated in its own script.

## 5. What is stored rather than recomputed

The bundle ships `results/`, one file per experiment, each carrying a provenance block that
records the generating script, the timestamp, the R version, the platform and the parameters
of the run. Every experimental value in the paper is emitted from those files by
`scripts/make_tables_v10_2.R` and `scripts/make_tables_v10_7.R`. No value is typed into the
source.

This means the paper can be rebuilt with `bash build.sh` from the stored results alone,
without any of the three datasets above. The datasets are needed only to regenerate the
result files themselves.

## 6. Licensing and redistribution

The NHANES files are United States government works in the public domain and may be
redistributed, but they are taken from the source here so that a reader gets the file the
provider currently serves rather than a copy of unknown age.

The Berlin open data files are published by the Amt fuer Statistik Berlin-Brandenburg under
its open data terms. Check the current terms on the landing page above before redistributing
the file.

Neither dataset contains direct identifiers. The NHANES files are public-use releases with
disclosure protection already applied, including the topcoding of age at eighty that the
documentation describes.
