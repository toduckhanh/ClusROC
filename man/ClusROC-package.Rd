\name{ClusROC-package}
\alias{ClusROC-package}
\alias{ClusROC}
\docType{package}
\title{
ROC analysis in three-class classification problems for clustered data
}
\description{
This package implements the approaches for ROC surface analysis in the setting of clustered data. In particular, the following quantities are estimated:
  (i) the covariate-specific ROC surface;
  (ii) the covariate-specific volume under ROC surface (VUS);
  (iii) the covariate-specific optimal pairs of thresholds.
  The selection methods are: Generalized Youden Index (GYI), Closest to Perfection (CtP) and Max Volume (MV).

See To et al. (2022) for more details.
}
\details{
\tabular{ll}{
Package: \tab ClusROC\cr
Type: \tab Package\cr
Version: \tab 1.0-1\cr
Date: \tab 2022-04-28\cr
License: \tab GPL 2 | GPL 3\cr
Lazy load: \tab yes\cr
}

Major functions are \code{\link{lme2}}, \code{\link{ROCsurface}}, \code{\link{optThres3}}, \code{\link{VUS}} and \code{\link{TCFs}}.

}

\author{
Duc-Khanh To, with contributions from Monica Chiogna and Gianfranco Adimari

Maintainer: Duc-Khanh To <toduc@stat.unipd.it>
}

\references{
Bantis, L. E., Nakas, C. T., Reiser, B., Myall, D., and Dalrymple-Alford, J. C. (2017). Construction of joint confidence regions for the optimal true class fractions of Receiver Operating Characteristic (ROC) surfaces and manifolds. \emph{Statistical methods in medical research}, \bold{26}, 3, 1429-1442.

Gurka, M. J., Edwards, L. J. , Muller, K. E., and Kupper, L. L. (2006). Extending the Box-Cox transformation to the linear mixed model. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, \bold{169}, 2, 273-288.

Gurka, M. J. and Edwards, L. J. (2011). Estimating variance components and random effects using the box-cox transformation in the linear mixed model. \emph{Communications in Statistics - Theory and Methods}, \bold{40}, 3, 515-531.

Kauermann, G. and Carroll, R. J. (2001). A note on the efficiency of sandwich covariance matrix estimation. \emph{Journal of the American Statistical Association}, \bold{96}, 456, 1387-1396.

Liang, K. Y. and Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika}, \bold{73}, 1, 13-22.

Mancl, L. A. and DeRouen, T. A. (2001). A covariance estimator for GEE with improved small-sample properties. \emph{Biometrics}, \bold{57}, 1, 126-134.

To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2022). Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data. \emph{Statistical Methods in Medical Research}, DOI: 10.1177/09622802221089029.

Xiong, C., Luo, J., Chen L., Gao, F., Liu, J., Wang, G., Bateman, R. and Morris, J. C. (2018). Estimating diagnostic accuracy for clustered ordinal diagnostic groups in the three-class case -- Application to the early diagnosis of Alzheimer disease. \emph{Statistical Methods in Medical Research}, \bold{27}, 3, 701-714.

}

\keyword{package}

