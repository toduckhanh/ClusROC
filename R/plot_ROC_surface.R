####==========================================================================####
## This file consists of functions for estimating and ploting the ROC surface   ##
## Date: 26/03/2021																															##
####==========================================================================####

#' @import rgl
#' @import utils

shade.ellips <- function(orgi, sig, lev){
  t1 <- sig[2, 2]
  sig[2, 2] <- sig[3, 3]
  sig[3, 3] <- t1
  t1 <- sig[1, 2]
  sig[1, 2] <- sig[1, 3]
  sig[1, 3] <- t1
  sig[lower.tri(sig)] <- sig[upper.tri(sig)]
  ellips <- ellipse3d(sig, centre = orgi[c(1,3,2)], t = sqrt(qchisq(lev, 3)))
  return(ellips)
}

### ---- Plot the ROC surface under normal assumption ----
#' @title Plot the covariate-specific ROC surface for clustered data after fitting the cluster-effect model.
#'
#' @description \code{ROCsurface} estimates and makes a 3D plot of covariate-specific ROC surface for clustered data.
#'
#' @param name.test  name of variable indicating diagnostic test of biomarker in data. It can be also the transformation, for example, \code{"log(y)"}, where the term \code{log} is for the log-transformation, and \code{y} is the name of test.
#' @param name.class  name of variable indicating disease classes (diagnostic groups) in data.
#' @param name.covars  an vector of names of covariates containing in data. The vector can contain also the transformation, for example, \code{c("x1", "sqrt(x2)", "I(x3^2)")}.
#' @param name.clust  name of variable indicating clusters in data.
#' @param data  a data frame containing the variables in the model.
#' @param levl.class  an vector of the unique values (as character strings) that (disease) class might have taken, sorted into increasing order of means of test results corresponding to the disease classes (diagnostic groups). If \code{levl.class = NULL}, the levels will be automatically determined based on data, and sorted into increasing order of means of test results corresponding to the disease classes (diagnostic groups).
#' @param boxcox  a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a Box-Cox transformation will be applied to the model to guarantee the normally assumptions.
#' @param apVar  a logical value. Default = \code{TRUE}. If set to \code{TRUE}, the covariance matrix for all estimated parameters in model with be obtained by using the sandwich formula.
#' @param interval_lambda  a vector containing the end-points of the interval to be searched for the Box-Cox parameter, \code{lambda}. Default = (-2, 2).
#' @param ...  additional arguments for \code{\link{lme}}, such as \code{control}, \code{contrasts}.
#'
#' @details
#' .....
#'
#'
#'
#' @return \code{lme2} returns an object of class inheriting from "lme2" class.
#'
#' The function \code{\link{print.lme2}} can be used to print a summary of the results.
#'
#' An object of class "lme2" is a list containing at least the following components:
#'
#' \item{call}{the matched call.}
#' \item{est_para}{the estimate of all parameters in model.}
#' \item{se_para}{the standard error, obtained by using the sandwich formula.}
#' \item{vcov_sand}{the estimated covariance matrix for all estimated parameters, obtained by the sandwich formula.}
#' \item{residual}{a list of the residuals}
#' \item{fitted}{a list of the fitted values.}
#' \item{randf}{a vector of the estimated random effects.}
#' \item{n_coef}{total numbers of coefficients included in model.}
#' \item{icc}{a estimate of intra-class correlation - ICC}
#' \item{boxcox}{logical value indicating whether the Box-Cox transformation was implemented or not.}
#'
#' @examples
#'

#' @export
ROCsurface <- function(out_lme2, x.val, step.tcf = 0.01, plot = TRUE, main = NULL, file.name = NULL,
                       ellips = FALSE, thresholds = NULL, ci.level = ifelse(ellips, 0.95, NULL)){
  # define parameters
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  if(out_lme2$n_coef/out_lme2$n_p != 3) stop("There is not a case of three-class setting!")
  n_p <- out_lme2$n_p
  par_model <- out_lme2$est_para
  beta_d <- matrix(par_model[1:(3*n_p)], ncol = 3, nrow = n_p, byrow = FALSE)
  sigma_d <- sqrt(par_model[(3*n_p + 2):(3*n_p + 4)]^2 + par_model[(3*n_p + 1)]^2)
  Z <- c(1, x.val)
  mu_d <- Z %*% beta_d
  a12 <- sigma_d[2]/sigma_d[1]
  a32 <- sigma_d[2]/sigma_d[3]
  b12 <- (mu_d[1] - mu_d[2])/sigma_d[1]
  b32 <- (mu_d[3] - mu_d[2])/sigma_d[3]
  # estimate the ROC surface by given p1 = TCF1 and p3 = TCF3
  p1 <- p3 <- seq(0, 1, by = step.tcf)
  dt_p1_p3 <- expand.grid(x = p1, y = p3)
  rocs <- apply(dt_p1_p3, 1, function(x){
    if(qnorm(x[1], mean = mu_d[1], sd = sigma_d[1]) <= qnorm(1 - x[2], mean = mu_d[3], sd = sigma_d[3])){
      res <- pnorm((qnorm(1 - x[2]) + b32)/a32) - pnorm((qnorm(x[1]) + b12)/a12)
    } else res <- 0
    return(res)
  })
  fit <- cbind(dt_p1_p3, rocs)
  colnames(fit) <- c("TCF1", "TCF3", "TCF2")
  out <- list()
  out$fit <- fit
  if(plot){
    tcf1 <- matrix(p1, length(p1), length(p1), byrow = FALSE)
    tcf3 <- matrix(p3, length(p1), length(p1), byrow = TRUE)
    tcf2 <- matrix(rocs, length(p1), length(p1), byrow = FALSE)
    open3d()
    par3d(windowRect = 50 + c(0, 0, 640, 640))
    if(is.null(main)) main <- "Covariate-specific ROC surface"
    plot3d(0, 0, 0, type = "n", box = FALSE, xlab = "TCF 1", ylab = "TCF 3", zlab = "TCF 2",
           xlim = c(0,1), ylim = c(0,1), zlim = c(0,1))
    bgplot3d({
      plot.new()
      title(main = main, line = 1)
    })
    surface3d(tcf1, tcf3, tcf2, col = "gray40", alpha = 0.5)
    if(ellips){
      if(is.null(thresholds)) stop("Need to assign the pair of thresholds!")
      if(inherits(thresholds, "numeric")){
        if(thresholds[1] >= thresholds[2]) stop("The 1st threshold need to less than 2nd threshold!")
        tcfs_ellips <- TCF_normal(par = par_model, x.val = x.val, thresholds = thresholds, n_p = n_p,
                                  boxcox = out_lme2$boxcox)
        vcov_tcfs <- TCF_normal_vcov(par_model = par_model, x.val = x.val, thresholds = thresholds,
                                     vcov_par_model = out_lme2$vcov_sand, n_p = n_p, fixed = TRUE,
                                     boxcox = out_lme2$boxcox)
      }
      ellip.tcf <- shade.ellips(orgi = tcfs_ellips, sig = vcov_tcfs, lev = ci.level)
      plot3d(ellip.tcf, box = FALSE, col = "green", alpha = 0.5, xlim = c(0,1),
             ylim = c(0, 1), zlim = c(0, 1), xlab = " ", ylab = " ", zlab = " ",
             add = TRUE)
      plot3d(tcfs_ellips[1], tcfs_ellips[3], tcfs_ellips[2], type = "s", col = "red",
             radius = 0.01, add = TRUE)
    }
    play3d(spin3d(axis = c(0, 0, 1), rpm = 12.25), duration = 2)
    play3d(spin3d(axis = c(0, 1, 0), rpm = 0.3), duration = 2)
    if(!is.null(file.name)){
      if(!grepl(".png", file.name)) file.name <- paste0(file.name, ".png")
      rgl.snapshot(file.name)
    }
    invisible(out)
  } else return(out)
}
