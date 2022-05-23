# ClusROC
The R package ClusROC
This package implements the techniques for ROC surface analysis, in cases of clustered data and in presence of covariates. In particular, the package allows to obtain covariate-specific point and interval estimation for:

  - [x] true class fractions (TCFs) at fixed pairs of thresholds;
  - [x] the ROC surface;
  - [x] the optimal pairs of thresholds;
  - [x] the volume under ROC surface (VUS).
  
Methods considered for the first three quantities are proposed and discussed in To et al. (2022). Referring to the third, three  different selection criteria are implemented: Generalized Youden Index (GYI), Closest to Perfection (CtP) and Maximum Volume (MV). Methods considered for the last are proposed and discussed in Xiong et al. (2018). Visualization tools are also provided. We refer readers to the articles cited above for all details.
