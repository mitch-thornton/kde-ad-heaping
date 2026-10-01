# scripts/

Everything that produced a number in the paper. Each `run_*.sh` wrapper sets the defaults
and calls its `.R` file, writing a result file and a log side by side into `../results/`.
The wrappers are bash, so invoke them as `bash run_whatever.sh` rather than executing them
under zsh.

All of them take `PKG_DIR`, which must point at a checkout that has an `R/` directory, that
is at this repository root rather than at a paper bundle.

## Experiments

| script | produces | what it is |
|---|---|---|
| `run_benchmark_v10_2.sh` | `benchmark_v10_2.rds` | the 50-seed benchmark behind Table 1 and Supplementary Tables S2 and S3, 76 cells, about 28 minutes |
| `run_batch_v10_3.sh` | `batch_v10_3.rds` | the threshold sweep, cutoff ablation, deconvolution bandwidth study and NHANES bandwidth sweep, Supplementary Tables S4 to S7 |
| `run_verify_readers_v10_3.sh` | `verify_readers_v10_3.rds` | grid reader, blind against verifying, and detector firing rates, Supplementary Table S9 |
| `run_endtoend_v10_4.sh` | `endtoend_v10_4.rds` | end-to-end reconstruction from an estimated grid, and the heaped-fraction known-truth test, Supplementary Tables S10 and S11 |
| `run_test_lattice_v3.sh` | `test_lattice_v3.rds` | the mixed-grain known-truth sweep against kappa, and the reading of the real cigarette counts, Supplementary Table S8 |
| `run_make_figures_v12.sh` | `fig1b_v10_7.rds`, `fig2_ise_v10_7.rds`, `nhanes_v10_7.rds` | Figures 1, 2 and 4, and the NHANES study at full precision |
| `run_sensitivity_n_v13.sh` | `sensitivity_n_v13.rds` | the sample-size sweep behind Supplementary Table S12, five sample sizes at one grid, about a second |
| `run_b6_fix_v10_3.sh` | `b6_fix_v10_3.rds` | the corrected attenuation simulation, gated on matching the real digit shares |

`run_sensitivity_n_v13.sh` differs from the others in two ways worth knowing. It reads the
target densities, the evaluation grid, the grid spacing and the error functional out of
`benchmark_v10_2.R` by parsing that file rather than carrying its own copies, so the two
cannot disagree about what the bimodal mixture is. And its draws are nested, each seed
drawing the benchmark's four thousand observations under the same seed rule with the smaller
samples taken as the leading `n`, so its largest sample is a published cell of Table 1 rather
than a fresh run of one. It checks that cell per seed against `benchmark_v10_2.rds` and exits
nonzero if it does not reproduce.

## Emitters

`run_make_tables_v10_2.sh` writes Table 1 and Supplementary Tables S1 to S3.
`run_make_tables_v12.sh` writes Table 2 and Supplementary Tables S4 to S11, and
`tables/numbers_v10_7.tex`, the macro file the paper and the response letter both read.
`run_make_tableS12_v13.sh` writes Supplementary Table S12 and `tables/numbers_v13.tex`, and
refuses to write at all if the stored sensitivity run did not pass its Table 1 check.

No experimental value in the paper is typed. Every one is emitted from a file in
`../results/`, each of which carries a provenance block recording the generating script,
the timestamp, the R version, the platform and the parameters of the run.

## Support

`mw15.R` defines the fifteen normal mixtures of Marron and Wand together with the two
custom targets, and `verify_mw15.R` checks those definitions against the published
parameters. `make_markup_v13.sh` builds the marked-up manuscript against the reviewed
version and is a submission tool rather than an experiment.

## Superseded, kept as the record

`make_tables_v10_7.sh`, `make_tables_v10_7.R`, `make_figures_v10_7.R`,
`run_make_figures_v10_7.sh` and `make_markup_v10_7.sh` are the previous generation of the
emitters and the markup tool. The v12 and v13 files above replace them and are what produced
the published tables and figures. The v10.7 copies are kept because earlier result files and
earlier versions of the manuscript were built with them.

`heap_lattice_v2.R` and `heap_lattice_v3.R` are the second and third mixed-grain readers.
The third is now the implementation of `heap_lattice()` in the package, so the copy here is
historical. The second is the nonnegative least-squares design that was tried and did not
survive its known-truth test; Methods records that history, so it is kept rather than
deleted. Neither should be sourced in preference to the package function.

`test_lattice_v2.R` is the known-truth test of that second reader, kept for the same reason.
