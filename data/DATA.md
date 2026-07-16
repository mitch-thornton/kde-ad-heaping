# Data

The package's own tests and examples use only synthetic data (fixed seed 20260627).

The two real datasets used in the paper are public and are not committed here to keep the
repository light:

* NHANES 2017-2018 (CDC/NCHS): DEMO_J.xpt, BMX_J.xpt, SMQ_J.xpt. Reproduce Table 2 with
  `Rscript nhanes_R.R <dir-with-xpt>`.
* Berlin resident register, December 2015 (Amt fuer Statistik Berlin-Brandenburg, open data):
  https://www.statistik-berlin-brandenburg.de/opendata/RBS_OD_LOR_2015_12.zip . Reproduce
  Table 3 with `Rscript berlin_R.R <dir-with-EWR201512E_Matrix.csv>`.
