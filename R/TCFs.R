####==========================================================================####
## This file consists of functions for estimating TCFs at a fixed               ##
## pair of thresholds                                                           ##
## Date: 21/10/2021																															##
####==========================================================================####


TCF_normal <- function(par, z, thresholds, n_p, boxcox = FALSE){
  beta_d <- par[1:(3*n_p)]
  sigma_d <- sqrt(par[(3*n_p + 2):(3*n_p + 4)]^2 + par[(3*n_p + 1)]^2)
  if(boxcox) thresholds <- boxcox_trans(thresholds, par[length(par)])
  mu_est <- z %*% beta_d
  res_tcfs <- numeric(3)
  res_tcfs[1] <- pnorm(thresholds[1], mean = mu_est[1], sd = sigma_d[1])
  res_tcfs[2] <- pnorm(thresholds[2], mean = mu_est[2], sd = sigma_d[2]) -
    pnorm(thresholds[1], mean = mu_est[2], sd = sigma_d[2])
  res_tcfs[3] <- 1 - pnorm(thresholds[2], mean = mu_est[3], sd = sigma_d[3])
  return(res_tcfs)
}

TCF_normal_vcov <- function(par_model, z, thresholds, vcov_par_model, n_p, fixed = FALSE, boxcox = FALSE,
                            type_thresholds = c("GYI", "CtP", "MV")){
  type_thresholds <- match.arg(type_thresholds)
  # if(boxcox) thresholds <- boxcox_trans(thresholds, par_model[length(par_model)])
  if(fixed){
    grad.tcfs <- jacobian(TCF_normal, x = par_model, z = z, thresholds = thresholds, n_p = n_p,
                          boxcox = boxcox)
    vcov.tcfs <- grad.tcfs %*% vcov_par_model %*% t(grad.tcfs)
  } else{
    grad.tcfs_par <- jacobian(TCF_normal, x = par_model, z = z, thresholds = thresholds, n_p = n_p,
                              boxcox = boxcox)
    grad.tcfs_thresholds <- jacobian(TCF_normal, x = thresholds, z = z, par = par_model, n_p = n_p)
    term3 <- switch (type_thresholds,
                     GYI = hessian(func = function(par, z, n_p){
                       GYI_normal(par = par[1:2], par_model = par[-c(1,2)], z = z, n_p = n_p)
                     }, x = c(thresholds, par_model), z = z, n_p = n_p),
                     CtP = hessian(func = function(par, z, n_p){
                       CtP_normal(par = par[1:2], par_model = par[-c(1,2)], z = z, n_p = n_p)
                     }, x = c(thresholds, par_model), z = z, n_p = n_p),
                     MV = hessian(func = function(par, z, n_p){
                       MV_normal(par = par[1:2], par_model = par[-c(1,2)], z = z, n_p = n_p)
                     }, x = c(thresholds, par_model), z = z, n_p = n_p)
    )
    term3.1 <- term3[1:2, 1:2]
    derv.cpt.para <- -solve(term3.1) %*% term3[1:2, -c(1,2)]
    derv.TCF.para <- grad.tcfs_par + grad.tcfs_thresholds %*% derv.cpt.para
    vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  }
  return(vcov.tcfs)
}

## Compute TCFs
### ---- Estimate TCFs under normal assumption ----
#' @title Estimation of the covariate-specific TCFs for clustered data.
#'
#' @description \code{TCFs} estimates covariate-specific True Class Fractions (TCFs), at a specified pair of thresholds, of a continuous diagnostic test in a clustered design with three ordinal groups. This function allows to estimate covariate-specific TCFs at multiple points for covariates.
#'
#' @param out_lme2  an object of class "lme2", a result of \code{\link{lme2}} call.
#' @param x.val  specific value(s) of covariate(s) where the TCFs are estimated. In absence of covariate, no values have to be specified. In case of one covariate, \code{x.val} should be a number. In case of \eqn{p} covariates (\eqn{p > 1}), \code{x.val} should be a vector containing \eqn{p} values; or a matrix with \eqn{p} columns and \eqn{m} rows containing values of the covariates if the user wants to estimate at \eqn{m} points.
#' @param thresholds  specified pair of thresholds.
#' @param apVar  logical value. If set to \code{TRUE}, the variance-covariance matrix of estimated covariate-specific TCFs is estimated.
#' @details
#' This function implements a method in To et al. (2022) for estimating covariate-specific TCFs at a specified pair of thresholds of a continuous diagnostic test in a clustered design with three ordinal groups. The estimator is based on results from  \code{\link{lme2}}, which uses the REML approach. The asymptotic variance-covariance matrix of the estimated covariate-specific TCFs is estimated through the Delta method. Note that, if the Box-Cox transformation is applied for the linear mixed-effect model, the pair of thresholds must be input in the original scale.
#'
#' Before performing estimation, a check for the monotone ordering assumption is performed. This means that, for the fixed values of covariates, three predicted mean values for test results in three diagnostic groups are compared. If the assumption is not meet, the covariate-specific TCFs at the values of covariates are not estimated.
#'
#' @return \code{TCFs} returns an object of class "TCFs", which is a list containing at least the following components:
#'
#' \item{call}{the matched call.}
#' \item{tcfs_est}{a vector or matrix containing the estimated TCFs.}
#' \item{tcf_vcov}{a matrix or list of matrices containing the estimated variance-covariance matrices.}
#' \item{thresholds}{specified pair of thresholds.}
#' \item{mess_order}{a diagnostic message from checking the monontone ordering.}
#' \item{x.val}{value(s) of covariate(s).}
#' \item{n_p}{total number of regressors in the model.}
#'
#' Generic functions such as \code{print} is also used to show the results.
#'
#'
#' @references
#' To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2022)
#' ``Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data''. \emph{Statistical Methods in Medical Research}, DOI: 10.1177/09622802221089029.
#'
#' @examples
#' data(data_3class)
#' ## One covariate
#' out1 <- lme2(fixed.formula = Y ~ X1, name.class = "D", name.clust = "id_Clus",
#'              data = data_3class)
#'
#' ### Estimate TCFs at one single value of X1, (t1, t2) = (1, 4)
#' out_tcfs_1 <- TCFs(out_lme2 = out1, x.val = 1, thresholds = c(1, 4), apVar = TRUE)
#' print(out_tcfs_1)
#'
#' ## Two covariates
#' out2 <- lme2(fixed.formula = Y ~ X1 + X2, name.class = "D", name.clust = "id_Clus",
#'              data = data_3class)
#'
#' ### Estimate covariate-specific TCFs at point (X1, X2) = (1, 0), and (t1, t2) = (1, 4)
#' out_tcfs_2 <- TCFs(out_lme2 = out2, x.val = c(1, 0), thresholds = c(1, 4), apVar = TRUE)
#' print(out_tcfs_2)
#'
#' ### Estimate covariate-specific TCFs at three points and (t1, t2) = (1, 4)
#' out_tcfs_3 <- TCFs(out_lme2 = out2, x.val = rbind(c(-0.5, 0), c(0.5, 0), c(0.5, 1)),
#'                    thresholds = c(1, 4), apVar = TRUE)
#' print(out_tcfs_3)
#'
#' @export
TCFs <- function(out_lme2, x.val, thresholds, apVar = FALSE){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  n_p <- out_lme2$n_p
  if(out_lme2$n_coef/n_p != 3) stop("There is not a case of three-class setting!")
  if(n_p == 1){
    if(!missing(x.val)) {
      if(!is.null(x.val)) warning("Sepecified value(s) of covariate(s) are not used!", call. = FALSE)
    }
    x.val <- NULL
  }
  if(n_p == 2){
    if(missing(x.val)) stop("Please input specific value(s) of covariate.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariate.")
    if(!inherits(x.val, "numeric")) stop("For the case of 1 covariate, please input a number or a vector.")
    if(any(is.na(x.val))) stop("NA value(s) are not allowed!")
  }
  if(n_p > 2){
    if(missing(x.val)) stop("Please input specific value(s) of covariates.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariates.")
    n_vb <- attr(out_lme2$terms, "n_vb")
    if(inherits(x.val, "numeric")){
      if(length(x.val) != n_vb) stop(paste("For case of", n_vb, "covariates, please input a vector of", n_vb, "values of covariates."))
    }
    if(inherits(x.val, "matrix")) {
      if(ncol(x.val) != n_vb) stop(paste("For case of m points of", n_vb, "covariates, please input a matrix with", n_vb, "columns and m rows containing values of covariates."))
    }
    if(any(is.na(x.val))) stop("NA value(s) not allowed!")
    x.val <- matrix(x.val, ncol = n_vb, byrow = FALSE)
  }
  if(apVar){
    if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing!")
    if(any(is.na(out_lme2$vcov_sand))) stop("There are NA values in the estimated covariance matrix of parameters. Unable to estimate standard error of TCFs.")
    # vcov_par_model <- out_lme2$vcov_sand[1:(out_lme2$n_coef + 4), 1:(out_lme2$n_coef + 4)]
  }
  call <- match.call()
  fit <- list()
  fit$call <- call
  ## main
  par_model <- out_lme2$est_para
  ## Check the ordering of means: mu_1 < mu_2 < mu_3
  Z <- make_data(out_lme2, x.val, n_p)
  res_check <- check_mu_order(Z, par_model, n_p)
  if(all(res_check$status == 0))
    stop("The assumption of montone ordering DOES NOT hold for all the value(s) of the covariate(s)")
  if(any(res_check$status == 0)){
    mess_order <- paste("The assumption of montone ordering DOES NOT hold for some points. The points number:",
                        paste(which(res_check$status == 0), collapse = ", "), "are deleted from analysis!")
    fit$mess_order <- mess_order
    message(mess_order)
  }
  Z <- res_check$Z_new
  if(n_p == 1){ # no covariate
    fit$x.val <- x.val
    fit$tcfs_est <- TCF_normal(par = par_model, z = Z, thresholds = thresholds, n_p = n_p,
                               boxcox = out_lme2$boxcox)
    if(apVar){
      fit$tcfs_cov <- list()
      fit$tcfs_cov[[1]] <- TCF_normal_vcov(par_model = par_model, z = Z, thresholds = thresholds,
                                           vcov_par_model = out_lme2$vcov_sand, n_p = n_p, fixed = TRUE,
                                           boxcox = out_lme2$boxcox)
    }
  }
  if(n_p == 2){ # 1 covariate
    fit$x.val <- x.val[res_check$status != 0]
    # tcfs_vector <- Vectorize(TCF_normal, vectorize.args = "x.val")
    fit$tcfs_est <- sapply(Z, function(x){
      TCF_normal(par = par_model, z = x, thresholds = thresholds, n_p = n_p, boxcox = out_lme2$boxcox)
    })
    if(apVar){
      # tcfs_se_vector <- Vectorize(TCF_normal_vcov, vectorize.args = "x.val", SIMPLIFY = FALSE)
      fit$tcfs_cov <- lapply(Z, function(x){
        TCF_normal_vcov(par_model = par_model, z = x, thresholds = thresholds,
                        vcov_par_model = out_lme2$vcov_sand, n_p = n_p, fixed = TRUE,
                        boxcox = out_lme2$boxcox)
      })
    }
  }
  if(n_p > 2){ # multiple covariates
    fit$x.val <- matrix(x.val[res_check$status != 0,], ncol = n_vb, byrow = FALSE)
    fit$tcfs_est <- sapply(Z, function(x){
      TCF_normal(par = par_model, z = x, thresholds = thresholds, n_p = n_p, boxcox = out_lme2$boxcox)
    })
    if(apVar){
      fit$tcfs_cov <- lapply(Z, function(x){
        TCF_normal_vcov(par_model = par_model, z = x, thresholds = thresholds,
                        vcov_par_model = out_lme2$vcov_sand, n_p = n_p, fixed = TRUE,
                        boxcox = out_lme2$boxcox)
      })
    }
  }
  fit$tcfs_est <- t(fit$tcfs_est)
  fit$thresholds <- thresholds
  fit$n_p <- n_p
  class(fit) <- "TCFs"
  return(fit)
}

## ---- The function print.TCFs ----
#' @title Print summary results from TCFs
#'
#' @description \code{print.TCFs} displays the results of the output from \code{\link{TCFs}}.
#'
#' @method print TCFs
#' @param x an object of class "TCFs", a result of \code{\link{TCFs}} call.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param call logical. If \code{TRUE}, the matched call will be printed.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.TCFs} shows a summary table for covariate-specific TCFs estimates.
#'
#' @seealso \code{\link{TCFs}}
#'
#' @export
print.TCFs <- function(x, digits = 3, call = TRUE, ...){
  if(isFALSE(inherits(x, "TCFs"))) stop("The object is not TCFs!")
  cat("\n")
  if(call){
    cat("CALL: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  }
  if(!is.null(x$mess_order)){
    cat("NOTE: ", x$mess_order, "\n \n", sep = "")
  }
  if(x$n_p == 1){
    labels <- "Intercept"
  }
  if(x$n_p == 2) {
    labels <- as.character(x$x.val)
  }
  if(x$n_p > 2) {
    labels <- apply(x$x.val, 1, function(y) paste0("(", paste(y, collapse = ", "), ")"))
  }
  infer_tab <- data.frame(labels, x$tcfs_est[,1], x$tcfs_est[,2], x$tcfs_est[,3])
  infer_tab[,2:4] <- signif(infer_tab[,2:4], digits = digits)
  colnames(infer_tab) <- c("Covariate(s) Values", "TCF 1", "TCF 2", "TCF 3")
  if(!is.null(x$tcfs_cov)){
    tcfs_se <- t(sapply(x$tcfs_cov, function(y) sqrt(diag(y))))
    infer_tab <- data.frame(infer_tab, tcfs_se[,1], tcfs_se[,2], tcfs_se[,3])
    infer_tab[,5:7] <- signif(infer_tab[,5:7], digits = digits)
    colnames(infer_tab) <- c("Covariate(s) Values", "TCF 1", "TCF 2", "TCF 3", "Se.TCF 1", "Se.TCF 2", "Se.TCF 3")
  }
  cat(paste0("Covariate-specific TCFs at (", x$thresholds[1], ",", x$thresholds[2], ")"), ": \n")
  print(infer_tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, ...)
  cat("\n")
  invisible(x)
}
