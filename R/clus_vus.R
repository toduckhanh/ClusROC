####========================================================================####
## This file consists of functions for estimating covariate-specific VUS      ##
####========================================================================####

#' @importFrom Rcpp evalCpp
#' @useDynLib ClusROC, .registration = TRUE
#' @import utils
#' @import numDeriv
#' @import stats

theta_1 <- function(par, z, n_p, n_c, n_k, n, subdivisions = 1000, ...) {
  beta_d <- par[1:(3 * n_p)]
  sigma_e <- par[(3 * n_p + 2):length(par)]
  mu_est <- z %*% beta_d
  f01 <- function(s) {
    f1 <- pnorm((s - mu_est[1]) / sigma_e[1])
    f2 <- dnorm((s - mu_est[2]) / sigma_e[2]) / sigma_e[2]
    f3 <- pnorm((s - mu_est[3]) / sigma_e[3], lower.tail = FALSE)
    return(f1 * f2 * f3)
  }
  rr <- mu_est[2] + c(-1, 1) * 6 * sigma_e[2]
  return(integrate(f01, lower = rr[1], upper = rr[2],
                   subdivisions = subdivisions, ...)$value)
}

theta_2 <- function(par, z, n_p, n_c, n_k, n, subdivisions = 1000, ...) {
  beta_d <- par[1:(3 * n_p)]
  sigma_e <- par[(3 * n_p + 2):length(par)]
  sigma_c <- par[(3 * n_p + 1)]
  mu_est <- z %*% beta_d
  f02 <- function(s) {
    f1 <- pnorm((s - mu_est[1]) / sigma_e[1])
    f2 <- dnorm((s - mu_est[2]) / sigma_e[2]) / sigma_e[2]
    f3 <- pnorm((s - mu_est[3]) / sqrt(2 * sigma_c^2 + sigma_e[3]^2),
                lower.tail = FALSE)
    return(f1 * f2 * f3)
  }
  rr <- mu_est[2] + c(-1, 1) * 6 * sigma_e[2]
  return(integrate(f02, lower = rr[1], upper = rr[2],
                   subdivisions = subdivisions, ...)$value)
}

theta_3 <- function(par, z, n_p, n_c, n_k, n, subdivisions = 1000, ...) {
  beta_d <- par[1:(3 * n_p)]
  sigma_e <- par[(3 * n_p + 2):length(par)]
  sigma_c <- par[(3 * n_p + 1)]
  mu_est <- z %*% beta_d
  f03 <- function(s) {
    f1 <- pnorm((s - mu_est[1]) / sigma_e[1])
    f2 <- dnorm((s - mu_est[2]) / sqrt(2 * sigma_c^2 + sigma_e[2]^2)) /
      sqrt(2 * sigma_c^2 + sigma_e[2]^2)
    f3 <- pnorm((s - mu_est[3]) / sigma_e[3], lower.tail = FALSE)
    return(f1 * f2 * f3)
  }
  rr <- mu_est[2] + c(-1, 1) * 6 * sqrt(2 * sigma_c^2 + sigma_e[2]^2)
  return(integrate(f03, lower = rr[1], upper = rr[2],
                   subdivisions = subdivisions, ...)$value)
}

theta_4 <- function(par, z, n_p, n_c, n_k, n, subdivisions = 1000, ...) {
  beta_d <- par[1:(3 * n_p)]
  sigma_e <- par[(3 * n_p + 2):length(par)]
  sigma_c <- par[(3 * n_p + 1)]
  mu_est <- z %*% beta_d
  f04 <- function(s) {
    f1 <- pnorm((s - mu_est[1]) / sqrt(2 * sigma_c^2 + sigma_e[1]^2))
    f2 <- dnorm((s - mu_est[2]) / sigma_e[2]) / sigma_e[2]
    f3 <- pnorm((s - mu_est[3]) / sigma_e[3], lower.tail = FALSE)
    return(f1 * f2 * f3)
  }
  rr <- mu_est[2] + c(-1, 1) * 6 * sigma_e[2]
  return(integrate(f04, lower = rr[1], upper = rr[2],
                   subdivisions = subdivisions, ...)$value)
}

theta_5 <- function(par, z, n_p, n_c, n_k, n, subdivisions = 1000, ...) {
  beta_d <- par[1:(3 * n_p)]
  sigma_e <- par[(3 * n_p + 2):length(par)]
  sigma_c <- par[(3 * n_p + 1)]
  mu_est <- z %*% beta_d
  f05 <- function(s) {
    f1 <- pnorm((s - mu_est[1]) / sqrt(sigma_c^2 + sigma_e[1]^2))
    f2 <- dnorm((s - mu_est[2]) / sqrt(sigma_c^2 + sigma_e[2]^2)) /
      sqrt(sigma_c^2 + sigma_e[2]^2)
    f3 <- pnorm((s - mu_est[3]) / sqrt(sigma_c^2 + sigma_e[3]^2),
                lower.tail = FALSE)
    return(f1 * f2 * f3)
  }
  rr <- mu_est[2] + c(-1, 1) * 6 * sqrt(sigma_c^2 + sigma_e[2]^2)
  return(integrate(f05, lower = rr[1], upper = rr[2],
                   subdivisions = subdivisions, ...)$value)
}

## vus_core function
vus_core <- function(par, z, n_p, n_c, n_k, n, p_sss, p_ssk, p_sks, p_skk,
                     p_ijk, subdivisions = 1000, ...) {
  theta_est_1 <- theta_1(par = par, z = z, n_p = n_p, n_c = n_c, n_k = n_k,
                         n = n, subdivisions = subdivisions, ...)
  theta_est_2 <- theta_2(par = par, z = z, n_p = n_p, n_c = n_c, n_k = n_k,
                         n = n, subdivisions = subdivisions, ...)
  theta_est_3 <- theta_3(par = par, z = z, n_p = n_p, n_c = n_c, n_k = n_k,
                         n = n, subdivisions = subdivisions, ...)
  theta_est_4 <- theta_4(par = par, z = z, n_p = n_p, n_c = n_c, n_k = n_k,
                         n = n, subdivisions = subdivisions, ...)
  theta_est_5 <- theta_5(par = par, z = z, n_p = n_p, n_c = n_c, n_k = n_k,
                         n = n, subdivisions = subdivisions, ...)
  vus_est <- theta_est_1 * p_sss + theta_est_2 * p_ssk + theta_est_3 * p_sks +
    theta_est_4 * p_skk + theta_est_5 * p_ijk
  return(vus_est)
}

vus_se <- function(par, vcov_par_model, z, n_p, n_c, n_k, n, p_sss, p_ssk,
                   p_sks, p_skk, p_ijk) {
  jac_vus <- rbind(jacobian(theta_1, x = par, z = z, n_p = n_p, n_c = n_c,
                            n_k = n_k, n = n),
                   jacobian(theta_2, x = par, z = z, n_p = n_p, n_c = n_c,
                            n_k = n_k, n = n),
                   jacobian(theta_3, x = par, z = z, n_p = n_p, n_c = n_c,
                            n_k = n_k, n = n),
                   jacobian(theta_4, x = par, z = z, n_p = n_p, n_c = n_c,
                            n_k = n_k, n = n),
                   jacobian(theta_5, x = par, z = z, n_p = n_p, n_c = n_c,
                            n_k = n_k, n = n))
  sig_1 <- jac_vus %*% vcov_par_model %*% t(jac_vus)
  pp <- c(p_sss, p_ssk, p_sks, p_skk, p_ijk)
  vus_sd <- as.numeric(sqrt(pp %*% sig_1 %*% pp))
  return(vus_sd)
}

### ---- Main function to estimate the covariate-specific VUS ----
#' @title Estimation of the covariate-specific VUS for clustered data.
#'
#' @description This function estimates the covariate-specific VUS of a continuous diagnostic test in the setting of clustered data as described in Xiong et al. (2018). This function allows to estimate covariate-specific VUS at multiple points for covariates.
#'
#' @param out_clus_lme  an object of class "clus_lme", a result of \code{\link{clus_lme}} call.
#' @param newdata   a data frame (containing specific value(s) of covariate(s)) in which to look for variables with which to estimate covariate-specific VUS. In absence of covariate, no values have to be specified.
#' @param ap_var  logical value. If set to \code{TRUE} (default), the standard error for (estimated) covariate-specific VUS are estimated.
#' @param subdivisions  the maximum number of subintervals used to approximate integral. Default is 1000.
#' @param ...  additional arguments to be passed to \code{\link[stats]{integrate}}.
#'
#' @details
#' This function implements a method in Xiong et al. (2018) for estimating covariate-specific VUS of a continuous diagnostic test in a clustered design with three ordinal groups. The estimator is based on results from \code{\link{clus_lme}}, which uses the REML approach. The standard error of the estimated covariate-specific VUS is approximated through the Delta method.
#'
#' Before performing estimation, a check for the monotone ordering assumption is performed. This means that, for the fixed values of covariates, three predicted mean values for test results in three diagnostic groups are compared. If the assumption is not meet, the covariate-specific VUS at the values of covariates are not estimated. In addition, this function also performs the statistical test, \eqn{H_0: VUS = 1/6} versus an alternative of interest.
#'
#'
#' @return \code{clus_vus} returns an object of class "VUS" which is a list containing at least the following components:
#'
#' \item{call}{the matched call.}
#' \item{vus_est}{a vector containing the estimated covariate-specific VUS.}
#' \item{vus_se}{a vector containing the standard errors.}
#' \item{mess_order}{a diagnostic message from checking the monontone ordering.}
#' \item{newdata}{value(s) of covariate(s).}
#' \item{n_p}{total number of regressors in the model.}
#'
#' Generic functions such as \code{print} is also used to show the results.
#'
#' @references
#' Xiong, C., Luo, J., Chen L., Gao, F., Liu, J., Wang, G., Bateman, R. and Morris, J. C. (2018)
#' ``Estimating diagnostic accuracy for clustered ordinal diagnostic groups in the three-class case -- Application to the early diagnosis of Alzheimer disease''.
#' \emph{Statistical Methods in Medical Research}, \bold{27}, 3, 701-714.
#'
#'
#' @examples
#' data(data_3class)
#' ## One covariate
#' out1 <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
#'                  name_clust = "id_Clus", data = data_3class)
#'
#' ### Estimate covariate-specific VUS at one value of one covariate
#' out_vus1 <- clus_vus(out1, newdata = data.frame(X1 = 0.5))
#' ci_clus_vus(out_vus1, ci_level = 0.95)
#'
#' ### Estimate covariate-specific VUS at multiple values of one covariate
#' out_vus2 <- clus_vus(out1, newdata = data.frame(X1 = c(-0.5, 0, 0.5)))
#' ci_clus_vus(out_vus2, ci_level = 0.95)
#'
#' ## Two covariates
#' out2 <- clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
#'                  name_clust = "id_Clus", data = data_3class)
#'
#' ### Estimate covariate-specific VUS at one point
#' out_vus3 <- clus_vus(out2, newdata = data.frame(X1 = 1.5, X2 = 1))
#' ci_clus_vus(out_vus3, ci_level = 0.95)
#'
#' ### Estimate covariate-specific VUS at three points
#' out_vus4 <- clus_vus(out2, newdata = data.frame(X1 = c(-0.5, 0.5, 0.5),
#'                                                 X2 = c(0, 0, 1)))
#' ci_clus_vus(out_vus4, ci_level = 0.95)
#'
#' @export
clus_vus <- function(out_clus_lme, newdata, ap_var = TRUE,
                     subdivisions = 1000, ...) {
  ## Check all conditions
  if (isFALSE(inherits(out_clus_lme, "clus_lme"))) {
    stop("out_clus_lme was not from clus_lme()!")
  }
  n_p <- out_clus_lme$n_p
  out_check_newdata <- check_newdata_vus(out_clus_lme$fixed_formula, newdata,
                                         n_p)
  newdata <- out_check_newdata$newdata
  ##
  if (ap_var) {
    if (is.null(out_clus_lme$vcov_sand)) {
      stop("The estimated covariance matrix of parameters was missing!")
    }
    if (any(is.na(out_clus_lme$vcov_sand))) {
      stop("There are NA values in the estimated covariance matrix of parameters. Unable to estimate standard error.")
    }
    vcov_par_model <- out_clus_lme$vcov_sand[1:(out_clus_lme$n_coef + 4),
                                             1:(out_clus_lme$n_coef + 4)]
  }
  call <- match.call()
  fit <- list()
  fit$call <- call
  ## check if all clusters/families have the same cluster size, if yes, the weight calculation_k is simplified
  n_c <- out_clus_lme$cls
  n_k <- out_clus_lme$n_c
  n <- out_clus_lme$n
  unique_n_k <- unique(n_k)
  equal_n_k_flag <- ifelse(length(unique_n_k) == 1, TRUE, FALSE)
  if (equal_n_k_flag) {
    n_k_div_n <- (unique_n_k / n)^3
    p_sss <- n_c * n_k_div_n
    p_ssk <- n_c * (n_c - 1) * n_k_div_n
    p_ijk <- n_c * (n_c - 1) * (n_c - 2) * n_k_div_n
  } else {
    n_k3 <- n_k^3
    p_sss <- sum(n_k3) / n^3
    n_k_sqr <- n_k^2
    outer_n_ksqr_n_k <- outer(n_k_sqr, n_k)
    diag(outer_n_ksqr_n_k) <- NA
    p_ssk <- sum(outer_n_ksqr_n_k, na.rm = TRUE) / (n^3)
    p_ijk <- fast_combn_sum(vals = 1:n_c, n = 3, n_k = n_k) / (n^3)
    p_ijk <- 6 * p_ijk
  }
  p_skk <- p_sks <- p_ssk
  par <- out_clus_lme$est_para[1:(out_clus_lme$n_coef + 4)]
  z <- make_data(out_clus_lme, newdata, n_p)
  ## Check the ordering of means: mu_1 < mu_2 < mu_3
  res_check <- check_mu_order(z, par, n_p)
  if (all(res_check$status == 0)) {
    stop("The assumption of montone ordering DOES NOT hold for all the value(s) of the covariate(s)")
  }
  if (any(res_check$status == 0)) {
    mess_order <- paste("The assumption of montone ordering DOES NOT hold for some points. The points number:",
                        paste(which(res_check$status == 0), collapse = ", "),
                        "are excluded from analysis!")
    fit$mess_order <- mess_order
    message(mess_order)
  }
  z <- res_check$z_new
  ##
  if (n_p == 1) { # without covariate
    fit$newdata <- newdata
    fit$vus_est <- vus_core(par = par, z = z[[1]], n_p = n_p, n_c = n_c,
                            n_k = n_k, n = n, p_sss = p_sss, p_ssk = p_ssk,
                            p_sks = p_sks, p_skk = p_skk, p_ijk = p_ijk,
                            subdivisions = subdivisions, ...)
    if (ap_var) {
      fit$vus_se <- vus_se(par = par, vcov_par_model = vcov_par_model,
                           z = z[[1]], n_p = n_p, n_c = n_c,
                           n_k = n_k, n = n, p_sss = p_sss, p_ssk = p_ssk,
                           p_sks = p_sks, p_skk = p_skk, p_ijk = p_ijk)
    }
  } else { # with covariate
    fit$newdata <- as.data.frame(newdata[res_check$status != 0, ])
    names(fit$newdata) <- names(newdata)
    fit$vus_est <- sapply(z, function(x) {
      vus_core(par = par, z = x, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
               p_sss = p_sss, p_ssk = p_ssk, p_sks = p_sks, p_skk = p_skk,
               p_ijk = p_ijk, subdivisions = subdivisions, ...)
      })
    if (ap_var) {
      fit$vus_se <- sapply(z, function(x) {
        vus_se(par = par, vcov_par_model = vcov_par_model, z = x, n_p = n_p,
               n_c = n_c, n_k = n_k, n = n, p_sss = p_sss, p_ssk = p_ssk,
               p_sks = p_sks, p_skk = p_skk, p_ijk = p_ijk)
        })
    }
  }
  fit$n_p <- n_p
  class(fit) <- "clus_vus"
  return(fit)
}

trunc_ci_vus <- function(ci) {
  if (ci[1] < 0) {
    ci[1] <- 0
  }
  if (ci[2] > 1) {
    ci[2] <- 1
  }
  return(ci)
}

## ---- The function ci_clus_vus ----
#' @title Confidence Intervals for Covariate-specific VUS
#'
#' @description Computes confidence intervals for covariate-specific VUS.
#'
#' @param x an object of class "VUS", a result of \code{\link{clus_vus}} call.
#' @param ci_level a confidence level to be used for constructing the confidence interval; default is 0.95.
#'
#' @details A confidence interval for covariate-specific VUS is given based on normal approximation. If the lower bound (or the upper bound) of the confidence interval is smaller than 0 (or greater than 1), it will be set as 0 (or 1). Also, logit and probit transformations are available if one wants guarantees that confidence limits are inside (0, 1).
#'
#' @return \code{ci_clus_vus} returns an object of class inheriting from "ci_VUS" class. An object of class "ci_VUS" is a list, containing at least the following components:
#'
#' \item{ci_vus_norm}{the normal approximation-based confidence interval for covariate-specific VUS.}
#' \item{ci_vus_log}{the confidence interval for covariate-specific VUS, after using logit-transformation.}
#' \item{ci_vus_prob}{the confidence interval for covariate-specific VUS, after using probit-transformation.}
#' \item{ci_level}{fixed confidence level.}
#' \item{newdata}{value(s) of covariate(s).}
#' \item{n_p}{total numbers of the regressors in the model.}
#'
#' @seealso \code{\link{clus_vus}}
#'
#' @export
ci_clus_vus <- function(x, ci_level = 0.95) {
  if (isFALSE(inherits(x, "clus_vus"))) {
    stop("The object is not clus_vus!")
  }
  if (is.null(x$vus_se)) {
    stop("Can not compute CI without standard error!")
  }
  n_p <- x$n_p
  fit <- list()
  c_alp <- qnorm((1 + ci_level) / 2)
  if (n_p == 1) { # no covariate
    ## Normal-approach with truncated boundary
    temp <- trunc_ci_vus(x$vus_est + c(-1, 1) * c_alp * x$vus_se)
    fit$ci_vus_norm <- matrix(temp, ncol = 2)
    ## logit-transform
    logit_vus <- qlogis(x$vus_est)
    logit_vus_sd <- x$vus_se / (x$vus_est * (1 - x$vus_est))
    fit$ci_vus_log <- matrix(plogis(logit_vus + c(-1, 1) * c_alp *
                                      logit_vus_sd), ncol = 2)
    ## probit-transform
    probit_vus <- qnorm(x$vus_est)
    probit_vus_sd <- grad(function(x) qnorm(x), x = x$vus_est) * x$vus_se
    fit$ci_vus_prob <- matrix(pnorm(probit_vus + c(-1, 1) * c_alp *
                                      probit_vus_sd), ncol = 2)
  }
  if (n_p == 2) { # 1 covariate
    fit$ci_vus_norm <- t(mapply(FUN = function(x, y) {
      trunc_ci_vus(x + c(-1, 1) * c_alp * y)
    }, x = x$vus_est, y = x$vus_se))
    ## logit-transform
    logit_vus <- qlogis(x$vus_est)
    logit_vus_sd <- x$vus_se / (x$vus_est * (1 - x$vus_est))
    fit$ci_vus_log <- t(mapply(FUN = function(x, y) {
      plogis(x + c(-1, 1) * c_alp * y)
    }, x = logit_vus, y = logit_vus_sd))
    ## probit-transform
    probit_vus <- qnorm(x$vus_est)
    probit_vus_sd <- grad(function(x) qnorm(x), x = x$vus_est) * x$vus_se
    fit$ci_vus_prob <- t(mapply(FUN = function(x, y) {
      pnorm(x + c(-1, 1) * c_alp * y)
    }, x = probit_vus, y = probit_vus_sd))
  }
  if (n_p > 2) { # multiple covariates
    fit$ci_vus_norm <- t(mapply(FUN = function(x, y) {
      trunc_ci_vus(x + c(-1, 1) * c_alp * y)
    }, x = x$vus_est, y = x$vus_se))
    ## logit-transform
    logit_vus <- qlogis(x$vus_est)
    logit_vus_sd <- x$vus_se / (x$vus_est * (1 - x$vus_est))
    fit$ci_vus_log <- t(mapply(FUN = function(x, y) {
      plogis(x + c(-1, 1) * c_alp * y)
      }, x = logit_vus, y = logit_vus_sd))
    ## probit-transform
    probit_vus <- qnorm(x$vus_est)
    probit_vus_sd <- grad(function(x) qnorm(x), x = x$vus_est) * x$vus_se
    fit$ci_vus_prob <- t(mapply(FUN = function(x, y) {
      pnorm(x + c(-1, 1) * c_alp * y)
    }, x = probit_vus, y = probit_vus_sd))
  }
  fit$ci_level <- ci_level
  fit$n_p <- n_p
  fit$newdata <- x$newdata
  class(fit) <- "ci_clus_vus"
  return(fit)
}

## ---- The function print.clus_vus ----
#' @title Print summary results from clus_vus
#'
#' @description \code{print.clus_vus} displays the results of the output from \code{\link{clus_vus}}.
#'
#' @method print clus_vus
#' @param x an object of class "VUS", a result of \code{\link{clus_vus}} call.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param call logical. If set to \code{TRUE}, the matched call will be printed.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.clus_vus} shows a summary table for covariate-specific VUS estimates, containing estimates, standard errors, z-values and p-values for the hypothesis testing \eqn{H_0: VUS = 1/6} versus an alternative \eqn{H_A: VUS > 1/6}.
#'
#' @return \code{print.clus_vus} returns a summary table for covariate-specific VUS estimates.
#'
#' @seealso \code{\link{clus_vus}}
#'
#' @export
print.clus_vus <- function(x, digits = 3, call = TRUE, ...) {
  if (isFALSE(inherits(x, "clus_vus"))) {
    stop("The object is not clus_vus!")
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
  if (!is.null(x$vus_se)) {
    z <- (x$vus_est - rep(1 / 6, length(x$vus_est))) / x$vus_se
    p_val <- pnorm(z, lower.tail = FALSE)
    infer_tab <- data.frame(labels, x$vus_est, x$vus_se, z, p_val)
    infer_tab[, 2:4] <- signif(infer_tab[, 2:4], digits = digits)
    pv <- as.vector(infer_tab[, 5])
    dig_tst <- max(1, min(5, digits - 1))
    infer_tab[, 5] <- format.pval(pv, digits = dig_tst, eps = 0.001)
    signif <- symnum(pv, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " "))
    sleg <- attr(signif, "legend")
    sleg <- strwrap(sleg, width = getOption("width") - 2, prefix = "  ")
    infer_tab <- cbind(infer_tab, format(signif))
    colnames(infer_tab) <- c("Covariates Values", "Est.", "Std.Error",
                             "z-value", "p-value", "")
    cat("Covariate-specific VUS: \n")
    print(infer_tab, quote = FALSE, right = TRUE, na.print = "--",
          row.names = FALSE, ...)
    cat("---\nSignif. codes:  ", sleg, sep = "",
        fill = getOption("width") + 4 + max(nchar(sleg, "bytes") - nchar(sleg)))
    cat("z-value and p-value are for testing the null hypothesis H0: VUS = 1/6 vs H1: VUS > 1/6 \n")
  } else {
    infer_tab <- data.frame(labels, x$vus_est)
    infer_tab[, 2] <- signif(infer_tab[, 2], digits = digits)
    colnames(infer_tab) <- c("Covariates Values", "Est.")
    cat("Covariate-specific VUS: \n")
    print(infer_tab, quote = FALSE, right = TRUE, na.print = "--",
          row.names = FALSE, ...)
  }
  cat("\n")
  invisible(x)
}

## ---- The function print.ci_clus_vus ----
#' @title Print summary results from ci_clus_vus
#'
#' @description \code{print.ci_vus} displays the results of the output from \code{\link{ci_clus_vus}}.
#'
#' @method print ci_clus_vus
#' @param x an object of class "ci_clus_vus", a result of \code{\link{ci_clus_vus}} call.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.ci_clus_vus} shows a summary table for confidence interval limits for covariate-specific VUS.
#'
#' @return \code{print.ci_clus_vus} shows a summary table for confidence intervals for covariate-specific VUS.
#'
#' @seealso \code{\link{clus_vus}}
#'
#' @export
print.ci_clus_vus <- function(x, digits = 3, ...) {
  if (isFALSE(inherits(x, "ci_clus_vus"))) {
    stop("The object is not ci_clus_vus!")
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
  ci_tab <- cbind(x$ci_vus_norm, x$ci_vus_log, x$ci_vus_prob)
  ci_tab <- format(round(ci_tab, digits = digits))
  res_ci_tab <- data.frame(labels,
                           apply(matrix(ci_tab[, 1:2], ncol = 2, byrow = FALSE),
                                 1, function(y) {
                                   paste0("(", paste(y, collapse = ", "), ")")
                                 }),
                           apply(matrix(ci_tab[, 3:4], ncol = 2, byrow = FALSE),
                                 1, function(y) {
                                   paste0("(", paste(y, collapse = ", "), ")")
                                 }),
                           apply(matrix(ci_tab[, 5:6], ncol = 2, byrow = FALSE),
                                 1, function(y) {
                                   paste0("(", paste(y, collapse = ", "), ")")
                                 }))
  colnames(res_ci_tab) <- c("Covariates Values", "Normal approximation",
                            "Logit transformation", "Probit transformation")
  cat(paste0("The ", x$ci.level * 100,
             "% confidence intervals for covariate-specific VUS:\n"))
  print(res_ci_tab, quote = FALSE, right = TRUE, na.print = "--",
        row.names = FALSE, print.gap = 3, ...)
  cat("\n")
  invisible(x)
}
