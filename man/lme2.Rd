% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_lme2.R
\name{lme2}
\alias{lme2}
\title{Fitting the cluster-effect models for two-class of three-class settings.}
\usage{
lme2(
  name.test,
  name.class,
  name.covars,
  name.clust,
  data,
  levl.class = NULL,
  boxcox = FALSE,
  apVar = TRUE,
  interval_lambda = c(-2, 2),
  trace = TRUE,
  ...
)
}
\arguments{
\item{name.test}{name of variable indicating diagnostic test of biomarker in data. It can be also the transformation, for example, \code{"log(y)"}, where the term \code{log} is for the log-transformation, and \code{y} is the name of test.}

\item{name.class}{name of variable indicating disease classes (diagnostic groups) in data.}

\item{name.covars}{an vector of names of covariates containing in data. The vector can contain also the transformation, for example, \code{c("x1", "sqrt(x2)", "I(x3^2)")}.}

\item{name.clust}{name of variable indicating clusters in data.}

\item{data}{a data frame containing the variables in the model.}

\item{levl.class}{an vector of the unique values (as character strings) that (disease) class might have taken, sorted into increasing order of means of test results corresponding to the disease classes (diagnostic groups). If \code{levl.class = NULL}, the levels will be automatically determined based on data, and sorted into increasing order of means of test results corresponding to the disease classes (diagnostic groups).}

\item{boxcox}{a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a Box-Cox transformation will be applied to the model to guarantee the normally assumptions.}

\item{apVar}{a logical value. Default = \code{TRUE}. If set to \code{TRUE}, the covariance matrix for all estimated parameters in model with be obtained by using the sandwich formula.}

\item{interval_lambda}{a vector containing the end-points of the interval to be searched for the Box-Cox parameter, \code{lambda}. Default = (-2, 2).}

\item{...}{additional arguments for \code{\link{lme}}, such as \code{control}, \code{contrasts}.}
}
\value{
\code{lme2} returns an object of class inheriting from "lme2" class.

The function \code{\link{print.lme2}} can be used to print a summary of the results.

An object of class "lme2" is a list containing at least the following components:

\item{call}{the matched call.}
\item{est_para}{the estimate of all parameters in model.}
\item{se_para}{the standard error, obtained by using the sandwich formula.}
\item{vcov_sand}{the estimated covariance matrix for all estimated parameters, obtained by the sandwich formula.}
\item{residual}{a list of the residuals}
\item{fitted}{a list of the fitted values.}
\item{randf}{a vector of the estimated random effects.}
\item{n_coef}{total numbers of coefficients included in the model.}
\item{n_p}{total numbers of regressors in the model.}
\item{icc}{a estimate of intra-class correlation - ICC}
\item{boxcox}{logical value indicating whether the Box-Cox transformation was implemented or not.}
}
\description{
\code{lme2} fits the cluster-effect models for two-class or three-class setting based on the \code{lme()} routine from \code{nlme}-package.
}
\details{
.....
}
\examples{
## Example for two-class setting
data(data_2class)
head(data_2class)
out1 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
             data = data_2class, apVar = FALSE)
print(out1)
plot(out1)
out2 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
             data = data_2class)
print(out2, digits = 3)

## Example for three-class setting
data(data_3class)
head(data_3class)
out3 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
             data = data_3class, apVar = FALSE)
print(out3)
out4 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
             data = data_3class)
print(out4)
plot(out4)

}
