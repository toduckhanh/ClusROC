####========================================================================####
## This file consists of functions for estimating TCFs at a fixed             ##
## pair of thresholds                                                         ##
####========================================================================####

tcf_normal <- function(par, z, thresholds, n_p, boxcox = FALSE) {
  beta_d <- par[1:(3 * n_p)]
  sigma_d <- sqrt(par[(3 * n_p + 2):(3 * n_p + 4)]^2 + par[(3 * n_p + 1)]^2)
  if (boxcox) {
    thresholds <- boxcox_trans(thresholds, par[length(par)])
  }
  mu_est <- z %*% beta_d
  res_tcfs <- numeric(3)
  res_tcfs[1] <- pnorm(thresholds[1], mean = mu_est[1], sd = sigma_d[1])
  res_tcfs[2] <- pnorm(thresholds[2], mean = mu_est[2], sd = sigma_d[2]) -
    pnorm(thresholds[1], mean = mu_est[2], sd = sigma_d[2])
  res_tcfs[3] <- 1 - pnorm(thresholds[2], mean = mu_est[3], sd = sigma_d[3])
  return(res_tcfs)
}

tcf_normal_vcov <- function(par_model, z, thresholds, vcov_par_model, n_p,
                            fixed = FALSE, boxcox = FALSE,
                            type_thresholds = c("GYI", "CtP", "MV")) {
  type_thresholds <- match.arg(type_thresholds)
  if (fixed) {
    grad_tcfs <- jacobian(tcf_normal, x = par_model, z = z,
                          thresholds = thresholds, n_p = n_p, boxcox = boxcox)
    vcov_tcfs <- grad_tcfs %*% vcov_par_model %*% t(grad_tcfs)
  } else {
    grad_tcfs_par <- jacobian(tcf_normal, x = par_model, z = z,
                              thresholds = thresholds, n_p = n_p,
                              boxcox = boxcox)
    grad_tcfs_thresholds <- jacobian(tcf_normal, x = thresholds, z = z,
                                     par = par_model, n_p = n_p)
    term3 <- switch(type_thresholds,
                    GYI = hessian(func = function(par, z, n_p) {
                      gyi_normal(par = par[1:2], par_model = par[-c(1, 2)],
                                 z = z, n_p = n_p)
                     }, x = c(thresholds, par_model), z = z, n_p = n_p),
                    CtP = hessian(func = function(par, z, n_p) {
                      ctp_normal(par = par[1:2], par_model = par[-c(1, 2)],
                                 z = z, n_p = n_p)
                     }, x = c(thresholds, par_model), z = z, n_p = n_p),
                    MV = hessian(func = function(par, z, n_p) {
                      mv_normal(par = par[1:2], par_model = par[-c(1, 2)],
                                z = z, n_p = n_p)
                     }, x = c(thresholds, par_model), z = z, n_p = n_p)
    )
    term3_1 <- term3[1:2, 1:2]
    derv_cpt_para <- -solve(term3_1) %*% term3[1:2, -c(1, 2)]
    derv_tcf_para <- grad_tcfs_par + grad_tcfs_thresholds %*% derv_cpt_para
    vcov_tcfs <- derv_tcf_para %*% vcov_par_model %*% t(derv_tcf_para)
  }
  return(vcov_tcfs)
}

## Compute TCFs
### ---- Estimate TCFs under normal assumption ----
#' @title Estimation of the covariate-specific TCFs for clustered data.
#'
#' @description \code{clus_tcfs} estimates covariate-specific True Class Fractions (TCFs), at a specified pair of thresholds, of a continuous diagnostic test in a clustered design with three ordinal groups. This function allows to estimate covariate-specific TCFs at multiple points for covariates.
#'
#' @param out_clus_lme  an object of class "clus_lme", a result of \code{\link{clus_lme}} call.
#' @param newdata  a data frame (containing specific value(s) of covariate(s)) in which to look for variables with which to estimate covariate-specific TCFs. In absence of covariate, no values have to be specified.
#' @param thresholds  a specified pair of thresholds.
#' @param ap_var  logical value. If set to \code{TRUE}, the variance-covariance matrix of estimated covariate-specific TCFs is estimated.
#' @details
#' This function implements a method in To et al. (2022) for estimating covariate-specific TCFs at a specified pair of thresholds of a continuous diagnostic test in a clustered design with three ordinal groups. The estimator is based on results from  \code{\link{clus_lme}}, which uses the REML approach. The asymptotic variance-covariance matrix of the estimated covariate-specific TCFs is estimated through the Delta method. Note that, if the Box-Cox transformation is applied for the linear mixed-effect model, the pair of thresholds must be input in the original scale.
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
#' \item{newdata}{value(s) of covariate(s).}
#' \item{n_p}{total number of regressors in the model.}
#'
#' Generic functions such as \code{print} is also used to show the results.
#'
#'
#' @references
#' To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2022)
#' ``Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data''. \emph{Statistical Methods in Medical Research}, \bold{7}, 31, 1325-1341.
#'
#' @examples
#' data(data_3class)
#' ## One covariate
#' out1 <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
#'                  name_clust = "id_Clus", data = data_3class)
#'
#' ### Estimate TCFs at one single value of X1, (t1, t2) = (1, 4)
#' out_tcfs_1 <- clus_tcfs(out_clus_lme = out1, newdata = data.frame(X1 = 1),
#'                         thresholds = c(1, 4), ap_var = TRUE)
#' print(out_tcfs_1)
#'
#' ## Two covariates
#' out2 <- clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
#'                  name_clust = "id_Clus", data = data_3class)
#'
#' ### Estimate covariate-specific TCFs at point (X1, X2) = (1, 0), and (t1, t2) = (1, 4)
#' out_tcfs_2 <- clus_tcfs(out_clus_lme = out2,
#'                         newdata = data.frame(X1 = 1, X2 = 0),
#'                         thresholds = c(1, 4), ap_var = TRUE)
#' print(out_tcfs_2)
#'
#' ### Estimate covariate-specific TCFs at three points and (t1, t2) = (1, 4)
#' out_tcfs_3 <- clus_tcfs(out_clus_lme = out2,
#'                         newdata = data.frame(X1 = c(-0.5, 0.5, 0.5),
#'                                              X2 = c(0, 0, 1)),
#'                         thresholds = c(1, 4), ap_var = TRUE)
#' print(out_tcfs_3)
#'
#' @export
clus_tcfs <- function(out_clus_lme, newdata, thresholds, ap_var = FALSE) {
  if (isFALSE(inherits(out_clus_lme, "clus_lme"))) {
    stop("out_clus_lme was not from lme2()!")
  }
  n_p <- out_clus_lme$n_p
  out_check_newdata <- check_newdata_vus(out_clus_lme$fixed_formula, newdata,
                                         n_p)
  newdata <- out_check_newdata$newdata
  if (ap_var) {
    if (is.null(out_clus_lme$vcov_sand)) {
      stop("The estimated covariance matrix of parameters was missing!")
    }
    if (any(is.na(out_clus_lme$vcov_sand))) {
      stop("There are NA values in the estimated covariance matrix of parameters. Unable to estimate standard error of TCFs.")
    }
  }
  call <- match.call()
  fit <- list()
  fit$call <- call
  ## main
  par_model <- out_clus_lme$est_para
  ## Check the ordering of means: mu_1 < mu_2 < mu_3
  z <- make_data(out_clus_lme, newdata, n_p)
  res_check <- check_mu_order(z, par_model, n_p)
  if (all(res_check$status == 0)) {
    stop("The assumption of montone ordering DOES NOT hold for all the value(s) of the covariate(s)")
  }
  if (any(res_check$status == 0)) {
    mess_order <- paste("The assumption of montone ordering DOES NOT hold for some points. The points number:",
                        paste(which(res_check$status == 0), collapse = ", "),
                        "are deleted from analysis!")
    fit$mess_order <- mess_order
    message(mess_order)
  }
  z <- res_check$z_new
  if (n_p == 1) { # with out covariate
    fit$newdata <- newdata
    fit$tcfs_est <- tcf_normal(par = par_model, z = z[[1]],
                               thresholds = thresholds, n_p = n_p,
                               boxcox = out_clus_lme$boxcox)
    if (ap_var) {
      fit$tcfs_cov <- list()
      fit$tcfs_cov[[1]] <- tcf_normal_vcov(
        par_model = par_model, z = z[[1]], thresholds = thresholds,
        vcov_par_model = out_clus_lme$vcov_sand, n_p = n_p, fixed = TRUE,
        boxcox = out_clus_lme$boxcox
        )
    }
  } else { # with covariate
    fit$newdata <- as.data.frame(newdata[res_check$status != 0, ])
    names(fit$newdata) <- names(newdata)
    fit$tcfs_est <- sapply(z, function(x) {
      tcf_normal(par = par_model, z = x, thresholds = thresholds, n_p = n_p,
                 boxcox = out_clus_lme$boxcox)
    })
    if (ap_var) {
      fit$tcfs_cov <- lapply(z, function(x) {
        tcf_normal_vcov(par_model = par_model, z = x, thresholds = thresholds,
                        vcov_par_model = out_clus_lme$vcov_sand, n_p = n_p,
                        fixed = TRUE, boxcox = out_clus_lme$boxcox)
      })
    }
  }
  fit$tcfs_est <- t(fit$tcfs_est)
  fit$thresholds <- thresholds
  fit$n_p <- n_p
  class(fit) <- "clus_tcfs"
  return(fit)
}

## ---- The function print.clus_tcfs ----
#' @title Print summary results from clus_tcfs
#'
#' @description \code{print.clus_tcfs} displays the results of the output from \code{\link{clus_tcfs}}.
#'
#' @method print clus_tcfs
#' @param x an object of class "clus_tcfs", a result of \code{\link{clus_tcfs}} call.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param call logical. If set to \code{TRUE}, the matched call will be printed.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.clus_tcfs} shows a summary table for covariate-specific TCFs estimates.
#'
#' @return \code{print.clus_tcfs} returns a summary table for covariate-specific TCFs estimates.
#'
#' @seealso \code{\link{clus_tcfs}}
#'
#' @export
print.clus_tcfs <- function(x, digits = 3, call = TRUE, ...) {
  if (isFALSE(inherits(x, "clus_tcfs"))) {
    stop("The object is not clus_tcfs!")
  }
  cat("\n")
  if (call) {
    cat("CALL: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n",
        sep = "")
  }
  if (!is.null(x$mess_order)) {
    cat("NOTE: ", x$mess_order, "\n \n", sep = "")
  }
  if (x$n_p == 1) {
    labels <- "Intercept"
  }
  if (x$n_p == 2) {
    labels <- apply(x$newdata, 1, function(y) paste0(y))
  }
  if (x$n_p > 2) {
    labels <- apply(x$newdata, 1, function(y) {
      paste0("(", paste(y, collapse = ", "), ")")
    })
  }
  infer_tab <- data.frame(labels, x$tcfs_est[, 1], x$tcfs_est[, 2],
                          x$tcfs_est[, 3])
  infer_tab[, 2:4] <- signif(infer_tab[, 2:4], digits = digits)
  colnames(infer_tab) <- c("Covariate(s) Values", "TCF 1", "TCF 2", "TCF 3")
  if (!is.null(x$tcfs_cov)) {
    tcfs_se <- t(sapply(x$tcfs_cov, function(y) sqrt(diag(y))))
    infer_tab <- data.frame(infer_tab, tcfs_se[, 1], tcfs_se[, 2], tcfs_se[, 3])
    infer_tab[, 5:7] <- signif(infer_tab[, 5:7], digits = digits)
    colnames(infer_tab) <- c("Covariate(s) Values", "TCF 1", "TCF 2", "TCF 3",
                             "Se.TCF 1", "Se.TCF 2", "Se.TCF 3")
  }
  cat(paste0("Covariate-specific TCFs at (", x$thresholds[1], ",",
             x$thresholds[2], ")"), ": \n")
  print(infer_tab, quote = FALSE, right = TRUE, na.print = "--",
        row.names = FALSE, ...)
  cat("\n")
  invisible(x)
}
