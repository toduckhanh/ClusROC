# ClusROC 1.0-2 (17th November 2022)

## Major changes

* Correct the names of variables and functions according to the Tidyverse style.
* Change the names of major functions, as listed below:
  - `lme2()` &rarr; `clus_lme()`
  - `VUS()` &rarr; `clus_vus()`
  - `ROCsurface()` &rarr; `clus_roc_surface()`
  - `TCFs()` &rarr; `clus_tcfs()`
  - `optThres3()` &rarr; `clus_opt_thres3()`
* Adding the test files.
* Adding EnergyEthiopia data file.
* Adding C/C++ code for supporting fast generate all combination when estimate 
covariate-specific VUS.
* Replace `doSNOW` and `snow` packages by `doParallel`.
* Replace `car` package by `ellipse`, in order to plot ellipse-based confidence
interval for optimal pair of thresholds.

## Minor changes

* Adding the COYING file, which contains the terms of GPL-3 license.

## Bug fixes

* `clus_lme()` could not work with `NA` values without `na.action`.
* `clus_lme()` cannot work with option of `subset`
* The argument `newdata` in all major functions were not well-defined.

# ClusROC 1.0-1 (25th May 2022)

*  First release.
