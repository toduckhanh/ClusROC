\name{ClusROC-package}
\alias{ClusROC-package}
\alias{ClusROC}
\docType{package}
\title{
ROC Analysis in Three-Class Classification Problems for Clustered Data
}
\description{
This package implements the techniques for ROC surface analysis, in cases of clustered data and in presence of covariates. In particular, the package allows to obtain covariate-specific point and interval estimation for:
  (i) true class fractions (TCFs) at fixed pairs of thresholds;
  (ii) the ROC surface;
  (iii) the volume under ROC surface (VUS);
  (iv) the optimal pairs of thresholds.
  Methods considered in points (i), (ii) and (iv) are proposed and discussed in To et al. (2022). Referring to point (iv), three  different selection criteria are implemented: Generalized Youden Index (GYI), Closest to Perfection (CtP) and Maximum Volume (MV). Methods considered in point (iii) are proposed and discussed in Xiong et al. (2018). Visualization tools are also provided. We refer readers to the articles cited above for all details.

}
\details{
\tabular{ll}{
Package: \tab ClusROC\cr
Type: \tab Package\cr
Version: \tab 1.0-2\cr
Date: \tab 2022-10-10\cr
License: \tab GPL 2 | GPL 3\cr
Lazy load: \tab yes\cr
}

Major functions are \code{\link{clus_lme}}, \code{\link{clus_roc_surface}}, \code{\link{clus_opt_thres3}}, \code{\link{clus_vus}} and \code{\link{clus_tcfs}}.

}

\author{
Duc-Khanh To, with contributions from Gianfranco Adimari and Monica Chiogna

Maintainer: Duc-Khanh To <toduc@stat.unipd.it>
}

\references{
Bantis, L. E., Nakas, C. T., Reiser, B., Myall, D., and Dalrymple-Alford, J. C. (2017). Construction of joint confidence regions for the optimal true class fractions of Receiver Operating Characteristic (ROC) surfaces and manifolds. \emph{Statistical methods in medical research}, \bold{26}, 3, 1429-1442.

Gurka, M. J., Edwards, L. J. , Muller, K. E., and Kupper, L. L. (2006). Extending the Box-Cox transformation to the linear mixed model. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, \bold{169}, 2, 273-288.

Gurka, M. J. and Edwards, L. J. (2011). Estimating variance components and random effects using the box-cox transformation in the linear mixed model. \emph{Communications in Statistics - Theory and Methods}, \bold{40}, 3, 515-531.

Kauermann, G. and Carroll, R. J. (2001). A note on the efficiency of sandwich covariance matrix estimation. \emph{Journal of the American Statistical Association}, \bold{96}, 456, 1387-1396.

Liang, K. Y. and Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika}, \bold{73}, 1, 13-22.

Mancl, L. A. and DeRouen, T. A. (2001). A covariance estimator for GEE with improved small-sample properties. \emph{Biometrics}, \bold{57}, 1, 126-134.

To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2022). Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data. \emph{Statistical Methods in Medical Research}, \bold{7}, 31, 1325-1341.

Xiong, C., Luo, J., Chen L., Gao, F., Liu, J., Wang, G., Bateman, R. and Morris, J. C. (2018). Estimating diagnostic accuracy for clustered ordinal diagnostic groups in the three-class case -- Application to the early diagnosis of Alzheimer disease. \emph{Statistical Methods in Medical Research}, \bold{27}, 3, 701-714.

}

\keyword{package}

