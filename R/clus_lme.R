####========================================================================####
## This file consists of functions for fitting the cluster-effect models      ##
## REML approach, based on the lme() routine                                  ##
####========================================================================####
#' @import numDeriv
#' @import nlme
#' @import stats
#' @import utils

### ---- reml log-likelihood for normal distribution ----
reml_loglik_vec <- function(par, d_list, y_list, z_list, v_list, cls, n_p,
                            n_class) {
  beta_fit <- par[1:(n_class * n_p)]
  sigma_c <- par[(n_class * n_p + 1)]
  sigma_e <- par[(n_class * n_p + 2):length(par)]
  sigma_all <- lapply(1:cls, function(i) {
    tem <- tcrossprod(v_list[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(d_list[[i]] %*% sigma_e^2),
                                  ncol = ncol(tem), nrow = ncol(tem)))
  })
  sig_chol <- lapply(sigma_all, chol)
  sig_chol_inv <- lapply(sig_chol, FUN = function(x) {
    backsolve(x, diag(1, ncol(x)))
  })
  sig_chol_inv_z <- mapply(FUN = function(x, y) crossprod(x, y),
                           x = sig_chol_inv, y = z_list, SIMPLIFY = FALSE)
  term_1 <- numeric(cls)
  for (i in 1:cls) {
    sig_chol_inv_y <- crossprod(sig_chol_inv[[i]], y_list[[i]])
    term_1[i] <- 0.5 * crossprod(sig_chol_inv_y -
                                   sig_chol_inv_z[[i]] %*% beta_fit) +
      sum(log(diag(sig_chol[[i]])))
  }
  tem <- Reduce("+", lapply(sig_chol_inv_z, FUN = crossprod))
  term_2 <- sum(log(diag(chol(tem)))) / cls
  return(term_1 + term_2 / 2)
}

reml_loglik_vec_sum <- function(par, d_list, y_list, z_list, v_list, cls, n_p,
                                n_class) {
  beta_fit <- par[1:(n_class * n_p)]
  sigma_c <- par[(n_class * n_p + 1)]
  sigma_e <- par[(n_class * n_p + 2):length(par)]
  sigma_all <- lapply(1:cls, function(i) {
    tem <- tcrossprod(v_list[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(d_list[[i]] %*% sigma_e^2),
                                  ncol = ncol(tem), nrow = ncol(tem)))
  })
  sig_chol <- lapply(sigma_all, chol)
  sig_chol_inv <- lapply(sig_chol, FUN = function(x) {
    backsolve(x, diag(1, ncol(x)))
  })
  sig_chol_inv_z <- mapply(FUN = function(x, y) crossprod(x, y),
                           x = sig_chol_inv, y = z_list, SIMPLIFY = FALSE)
  term_1 <- numeric(cls)
  for (i in 1:cls) {
    sig_chol_inv_y <- crossprod(sig_chol_inv[[i]], y_list[[i]])
    term_1[i] <- 0.5 * crossprod(sig_chol_inv_y -
                                   sig_chol_inv_z[[i]] %*% beta_fit) +
      sum(log(diag(sig_chol[[i]])))
  }
  tem <- Reduce("+", lapply(sig_chol_inv_z, FUN = crossprod))
  term_2 <- sum(log(diag(chol(tem))))
  return(sum(term_1) + term_2 / 2)
}

### ---- reml log-likelihood for Box-Cox transformation ----
reml_bcx_loglik_vec <- function(par, d_list, y_list, z_list, v_list, cls, n_p,
                                n_class) {
  beta_fit <- par[1:(n_class * n_p)]
  sigma_c <- par[(n_class * n_p + 1)]
  sigma_e <- par[(n_class * n_p + 2):(length(par) - 1)]
  lambda <- par[length(par)]
  y_boxcox <- lapply(y_list, function(x) boxcox_trans(x, lambda))
  sigma_all <- lapply(1:cls, function(i) {
    tem <- tcrossprod(v_list[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(d_list[[i]] %*% sigma_e^2),
                                  ncol = ncol(tem), nrow = ncol(tem)))
  })
  sig_chol <- lapply(sigma_all, chol)
  sig_chol_inv <- lapply(sig_chol, FUN = function(x) {
    backsolve(x, diag(1, ncol(x)))
  })
  sig_chol_inv_z <- mapply(FUN = function(x, y) crossprod(x, y),
                           x = sig_chol_inv, y = z_list, SIMPLIFY = FALSE)
  term_1 <- numeric(cls)
  term_3 <- numeric(cls)
  for (i in 1:cls) {
    sig_chol_inv_y <- crossprod(sig_chol_inv[[i]], y_boxcox[[i]])
    term_1[i] <- 0.5 * crossprod(sig_chol_inv_y -
                                   sig_chol_inv_z[[i]] %*% beta_fit) +
      sum(log(diag(sig_chol[[i]])))
    term_3[i] <- sum(log(y_list[[i]]))
  }
  tem <- Reduce("+", lapply(sig_chol_inv_z, FUN = crossprod))
  term_2 <- sum(log(diag(chol(tem)))) / cls
  return(term_1 + term_2 / 2 - (lambda - 1) * term_3)
}

reml_bcx_loglik_vec_sum <- function(par, d_list, y_list, z_list, v_list, cls,
                                    n_p, n_class) {
  beta_fit <- par[1:(n_class * n_p)]
  sigma_c <- par[(n_class * n_p + 1)]
  sigma_e <- par[(n_class * n_p + 2):(length(par) - 1)]
  lambda <- par[length(par)]
  y_boxcox <- lapply(y_list, function(x) boxcox_trans(x, lambda))
  sigma_all <- lapply(1:cls, function(i) {
    tem <- tcrossprod(v_list[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(d_list[[i]] %*% sigma_e^2),
                                  ncol = ncol(tem), nrow = ncol(tem)))
  })
  sig_chol <- lapply(sigma_all, chol)
  sig_chol_inv <- lapply(sig_chol, FUN = function(x) {
    backsolve(x, diag(1, ncol(x)))
  })
  sig_chol_inv_z <- mapply(FUN = function(x, y) crossprod(x, y),
                           x = sig_chol_inv, y = z_list, SIMPLIFY = FALSE)
  term_1 <- numeric(cls)
  term_3 <- numeric(cls)
  for (i in 1:cls) {
    sig_chol_inv_y <- crossprod(sig_chol_inv[[i]], y_boxcox[[i]])
    term_1[i] <- 0.5 * crossprod(sig_chol_inv_y -
                                   sig_chol_inv_z[[i]] %*% beta_fit) +
      sum(log(diag(sig_chol[[i]])))
    term_3[i] <- sum(log(y_list[[i]]))
  }
  tem <- Reduce("+", lapply(sig_chol_inv_z, FUN = crossprod))
  term_2 <- sum(log(diag(chol(tem))))
  return(sum(term_1) + term_2 / 2 - (lambda - 1) * sum(term_3))
}

llike_bcx_fun <- function(par, fixed, random, weights, data, y_tit, ...) {
  y <- model.response(model.frame(fixed, data = data))
  y_boxcox <- boxcox_trans(y, par)
  w_boxcox <- y_boxcox / (y_tit^(par - 1))
  fixed_new <- update(fixed, w_boxcox  ~ .)
  data$w_boxcox <- w_boxcox
  out_model <- lme(fixed = fixed_new, random = random, weights = weights,
                   method = "REML", data = data, ...)
  return(as.numeric(out_model$logLik))
}

### ---- Fitting cluster-effect models ----

#' @title Linear Mixed-Effects Models for a continuous diagnostic test or a biomarker (or a classifier).
#'
#' @description \code{clus_lme} fits the cluster-effect model for a continuous diagnostic test in a three-class setting as described in Xiong et al. (2018) and To et al. (2022).
#'
#' @param fixed_formula  a two-sided linear formula object, describing the fixed-effects part of the model for three classes, with the response on the left of ~ operator and the terms, separated by + operators, on the right. For example, \code{Y ~ X1 + X2}, \code{Y ~ X1 + X2 + X1:X2} or \code{log(Y) ~ X1 + X2 + I(X1^2)}.
#' @param name_class  name of variable indicating disease classes (or diagnostic groups) in the data.
#' @param name_clust  name of variable indicating clusters in the data.
#' @param data  a data frame containing the variables in the model.
#' @param levl_class  a vector (of strings) containing the ordered name chosen for the disease classes. The ordering is intended to be ``increasing'' with respect to the disease severity. If \code{levl_class = NULL} (default), the elements of the vector will be automatically determined from data, by considering the order of the means of the test values for each disease class (diagnostic group).
#' @param ap_var  a logical value. Default = \code{TRUE}. If set to \code{TRUE}, the estimated covariance matrix for all estimated parameters in the model will be obtained (by using the sandwich formula).
#' @param boxcox  a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a Box-Cox transformation will be applied to the model.
#' @param interval_lambda  a vector containing the end-points of the interval for searching the Box-Cox parameter, \code{lambda}. Default = (-2, 2).
#' @param trace  a logical value. Default = \code{TRUE}. If set to \code{TRUE}, the information about the check for the monotonic ordering of test values will be provided.
#' @param ...  additional arguments for \code{\link[nlme]{lme}}, such as \code{control}, \code{contrasts}.
#'
#' @details
#' This function fits a linear mixed-effect model for a continuous diagnostic test in a three-class setting in order to account for the cluster and covariates effects on the test result. See Xiong et al. (2018) and To et al. (2022) for more details.
#' \itemize{
#' \item Estimation is done by using \code{\link[nlme]{lme}} with the restricted maximum log-likelihood (REML) method.
#' \item Box-Cox transformation for the model can be used when the distributions of test results are skewed (Gurka et al. 2006). The estimation procedure is described in To et al. (2022). The Box-Cox parameter \eqn{\lambda} is estimated by a grid search on the interval (-2, 2), as discussed in Gurka and Edwards (2011).
#' \item The estimated variance-covariance matrix for the estimated parameters are obtained by sandwich formula (see, Liang and Zeger, 1986; Kauermann and Carroll, 2001; Mancl and DeRouen, 2001) as discussed in To et al. (2022).
#' }
#'
#'
#' @return \code{clus_lme} returns an object of class "clus_lme" class, i.e., a list containing at least the following components:
#'
#' \item{call}{the matched call.}
#' \item{est_para}{a vector containing the estimated parameters.}
#' \item{se_para}{a vector containing the standard errors.}
#' \item{vcov_sand}{the estimated covariance matrix for all estimated parameters.}
#' \item{residual}{a list of residuals.}
#' \item{fitted}{a list of fitted values.}
#' \item{randf}{a vector of estimated random effects for each cluster level.}
#' \item{n_coef}{total number of coefficients included in the model.}
#' \item{n_p}{total numbers of regressors in the model.}
#' \item{icc}{an estimate of intra-class correlation - ICC}
#' \item{terms}{the \code{\link[stats]{terms}} object used.}
#' \item{boxcox}{logical value indicating whether the Box-Cox transformation was applied or not.}
#'
#' Generic functions such as \code{print} and \code{plot} are also used to show results of the fit.
#'
#' @references
#' Gurka, M. J., Edwards, L. J. , Muller, K. E., and Kupper, L. L. (2006) ``Extending the Box-Cox transformation to the linear mixed model''. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, \bold{169}, 2, 273-288.
#'
#' Gurka, M. J. and Edwards, L. J. (2011) ``Estimating variance components and random effects using the box-cox transformation in the linear mixed model''. \emph{Communications in Statistics - Theory and Methods}, \bold{40}, 3, 515-531.
#'
#' Kauermann, G. and Carroll, R. J. (2001)
#' ``A note on the efficiency of sandwich covariance matrix estimation''.
#' \emph{Journal of the American Statistical Association}, \bold{96}, 456, 1387-1396.
#'
#' Liang, K. Y. and Zeger, S. L. (1986)
#' ``Longitudinal data analysis using generalized linear models''. \emph{Biometrika}, \bold{73}, 1, 13-22.
#'
#' Mancl, L. A. and DeRouen, T. A. (2001) ``A covariance estimator for GEE with improved small-sample properties''.
#'  \emph{Biometrics}, \bold{57}, 1, 126-134.
#'
#' To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2022)
#' ``Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data''. \emph{Statistical Methods in Medical Research}, \bold{7}, 31, 1325-1341.
#'
#' Xiong, C., Luo, J., Chen L., Gao, F., Liu, J., Wang, G., Bateman, R. and Morris, J. C. (2018)
#' ``Estimating diagnostic accuracy for clustered ordinal diagnostic groups in the three-class case -- Application to the early diagnosis of Alzheimer disease''.
#' \emph{Statistical Methods in Medical Research}, \bold{27}, 3, 701-714.
#'
#'
#' @examples
#' ## Example 1:
#' data(data_3class)
#' head(data_3class)
#' ## A model with two covariate: X1 + X2
#' out1 <- clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D", name_clust = "id_Clus",
#'                  data = data_3class)
#' print(out1)
#' plot(out1)
#'
#' \donttest{
#' ## Example 2: Box-Cox transformation
#' data(data_3class_bcx)
#' out2 <- clus_lme(fixed_formula = Y ~ X, name_class = "D", name_clust = "id_Clus",
#'                  data = data_3class_bcx, boxcox = TRUE)
#' print(out2)
#' plot(out2)
#' }
#'
#' @export
clus_lme <- function(fixed_formula, name_class, name_clust, data,
                     levl_class = NULL, ap_var = TRUE, boxcox = FALSE,
                     interval_lambda = c(-2, 2), trace = TRUE, ...) {
  if (missing(data)) {
    data <- .GlobalEnv
    cat("Warning: the data is missing, the global environment is used!\n")
  }
  data <- as.data.frame(data)
  if (missing(fixed_formula)) {
    stop("argument fixed_formula is missing with no default")
  }
  if (!inherits(fixed_formula, "formula") || length(fixed_formula) != 3) {
    stop("\nfixed-effects model must be a formula of the form \"resp ~ pred\"")
  }
  if (missing(name_class)) {
    stop("argument name_class is missing with no default")
  }
  if (!inherits(name_class, "character") || length(name_class) != 1) {
    stop("agrument name_class must be a character vector with length 1.")
  }
  if (!is.element(name_class, names(data))) {
    stop(paste("Could not find name_class:", name_class, "in the input data."))
  }
  if (missing(name_clust)) {
    stop("argument name_clust is missing with no default")
  }
  if (!inherits(name_clust, "character") || length(name_clust) != 1) {
    stop("agrument name_clust must be a character vector with length 1.")
  }
  if (!is.element(name_clust, names(data))) {
    stop(paste("Could not find name_clust:", name_clust, "in the input data."))
  }
  ##
  n_class <- length(table(data[, name_class]))
  if (n_class != 3) {
    stop("There is not a case of three-class setting!")
  }
  mf <- unlist(strsplit(as.character(fixed_formula), "~"))[-1]
  name_test <- mf[1]
  if (boxcox) {
    if (any(model.extract(model.frame(fixed_formula, data), "response") < 0)) {
      stop("Cannot apply Box-Cox transform for negative values.")
    }
  }
  form_mean <- as.formula(paste(name_test, "~", name_class))
  mean_temp <- aggregate(form_mean, FUN = mean, data = data)
  temp_levl <- mean_temp[order(mean_temp[, 2]), 1]
  if (is.null(levl_class)) {
    if (trace) {
      cat("The ordered levels of disease classes are specified by the order of \n the means of the test values for each disease class:\n")
      cat(paste(temp_levl, collapse = " < "), "\n")
    }
    levl_class <- temp_levl
  } else {
    if (!inherits(levl_class, "character") || length(levl_class) != 3) {
      stop("agrument levl_class must be a character vector with length 3.")
    }
    if (all(levl_class == temp_levl)) {
      if (trace) {
        cat("The orders of inputed levels of disease classes are the same as \n the one obtained by the orders of averages of tests results:\n")
        cat(paste(levl_class, collapse = " < "), "\n")
      }
    } else {
      if (trace) {
        cat("The orders of inputed levels of disease classes are not the same as \n the one obtained by the orders of averages of tests results:\n")
        cat("The correct one should be:\n")
        cat(paste(temp_levl, collapse = " < "), "\n")
      }
      levl_class <- temp_levl
    }
  }
  data[, name_class] <- factor(data[, name_class], levels = levl_class)
  name_covars <- mf[2]
  ## define the formulas for fixed, random and weights
  call <- match.call()
  fit <- list()
  fit$call <- call
  fit$boxcox <- boxcox
  fit$name_test <- name_test
  fit$name_class <- name_class
  fit$name_clust <- name_clust
  fixed <- as.formula(paste(name_test, "~", name_class, "+",
                            "(", name_covars, ")", ":", name_class, "-1"))
  fit$name_covars <- name_covars
  random <- as.formula(paste("~", "1|", name_clust))
  form_weights <- as.formula(paste("~", "1|", name_class))
  weights <- varIdent(form = form_weights)
  ##
  md_frame <- model.frame(fixed, data)
  md_frame[, name_class] <- factor(md_frame[, name_class], levels = levl_class)
  fit$terms <- terms(fixed)
  xx <- model.matrix(fit$terms, md_frame, contrasts.arg = NULL)
  attr(fit$terms, "levl_class") <- levl_class
  attr(fit$terms, "n_vb") <- length(as.character(attr(fit$terms,
                                                      "variables"))[-c(1:3)])
  attr(fit$terms, "xlevels") <- .getXlevels(fit$terms, md_frame)
  attr(fit$terms, "contrasts") <- attr(xx, "contrasts")
  n <- nrow(data)
  clus <- model.frame(getGroupsFormula(random), data = data)[, 1]
  n_c <- table(clus)
  cls <- length(unique(clus))
  if (boxcox) {
    all_y <- model.response(model.frame(fixed, data = data))
    ## check ok!
    list_y <- split(all_y, clus)
    y_tit <- prod(sapply(list_y, function(x) prod(x^(1 / n))))
    lambda_est <- optimize(llike_bcx_fun, interval = interval_lambda,
                           fixed = fixed, random = random, weights = weights,
                           data = data, y_tit = y_tit,
                           maximum = TRUE, ...)$maximum
    all_y_boxcox <- boxcox_trans(all_y, lambda_est)
    data$y_boxcox <- all_y_boxcox
    fixed_new <- update(fixed, y_boxcox  ~ .)
    out_model <- lme(fixed = fixed_new, random = random, weights = weights,
                     method = "REML", data = data, ...)
  } else {
    out_model <- lme(fixed = fixed, random = random, weights = weights,
                     method = "REML", data = data, ...)
  }
  ## collecting results
  n_coef <- length(out_model$coefficients$fixed)
  n_p <- n_coef / n_class
  if (n_class == 2) {
    id_coef <- c(seq(1, n_coef - 1, by = 2), seq(2, n_coef, by = 2))
  }
  if (n_class == 3) {
    id_coef <- c(seq(1, n_coef - 2, by = 3), seq(2, n_coef - 1, by = 3),
                 seq(3, n_coef, by = 3))
  }
  sigma_e_est <- coef(out_model$modelStruct$varStruct, unconstrained = FALSE,
                      allCoef = TRUE) * out_model$sigma
  sigma_e_est <- sigma_e_est[as.character(levl_class)]
  sigma_c_est <- sqrt(as.numeric(getVarCov(out_model)))
  if (boxcox) {
    par_est <- c(out_model$coefficients$fixed[id_coef], sigma_c_est,
                 sigma_e_est, lambda_est)
    names(par_est) <- c(names(par_est)[1:n_coef], "sigma_c",
                        paste0("sigma_", c(1:n_class)), "lambda")
  } else {
    par_est <- c(out_model$coefficients$fixed[id_coef], sigma_c_est,
                 sigma_e_est)
    names(par_est) <- c(names(par_est)[1:n_coef], "sigma_c",
                        paste0("sigma_", c(1:n_class)))
  }
  icc <- sigma_c_est^2 / (sigma_c_est^2 + mean(sigma_e_est)^2)
  fit$n_coef <- n_coef
  fit$n_p <- n_p
  fit$n <- n
  fit$cls <- cls
  fit$n_c <- n_c
  fit$est_para <- par_est
  fit$icc <- icc
  ##
  resid <- out_model$residual[, 2]
  fit$residual <- split(resid, data[, name_class])
  fitted <- out_model$fitted[, 2]
  fit$fitted <- split(fitted, data[, name_class])
  fit$randf <- ranef(out_model)$`(Intercept)`
  if (ap_var) {
    if (boxcox) {
      data_matrix <- model.matrix(out_model$terms, data = data)
      all_d <- data_matrix[, 1:n_class]
      all_z <- data_matrix[, id_coef]
      y_list <- split(all_y, clus)
      d_list <- lapply(split(as.data.frame(all_d), clus), as.matrix)
      z_list <- lapply(split(as.data.frame(all_z), clus), as.matrix)
      v_list <- lapply(n_c, function(x) rep(1, x))
      jac <- jacobian(reml_bcx_loglik_vec, x = par_est, d_list = d_list,
                      y_list = y_list, z_list = z_list, v_list = v_list,
                      cls = cls, n_p = n_p, n_class = n_class)
      hes <- hessian(reml_bcx_loglik_vec_sum, x = par_est, d_list = d_list,
                     y_list = y_list, z_list = z_list, v_list = v_list,
                     cls = cls, n_p = n_p, n_class = n_class)
      vcov_sand <- solve(hes) %*% matrix(rowSums(apply(jac, 1, tcrossprod)),
                                         n_coef + n_class + 2,
                                         n_coef + n_class + 2) %*% solve(hes)
    } else {
      all_y <- model.response(model.frame(out_model$terms, data = data))
      data_matrix <- model.matrix(out_model$terms, data = data)
      all_d <- data_matrix[, 1:n_class]
      all_z <- data_matrix[, id_coef]
      y_list <- split(all_y, clus)
      d_list <- lapply(split(as.data.frame(all_d), clus), as.matrix)
      z_list <- lapply(split(as.data.frame(all_z), clus), as.matrix)
      v_list <- lapply(n_c, function(x) rep(1, x))
      jac <- jacobian(reml_loglik_vec, x = par_est, d_list = d_list,
                      y_list = y_list, z_list = z_list, v_list = v_list,
                      cls = cls, n_p = n_p, n_class = n_class)
      hes <- hessian(reml_loglik_vec_sum, x = par_est, d_list = d_list,
                     y_list = y_list, z_list = z_list, v_list = v_list,
                     cls = cls, n_p = n_p, n_class = n_class)
      vcov_sand <- solve(hes) %*% matrix(rowSums(apply(jac, 1, tcrossprod)),
                                         n_coef + n_class + 1,
                                         n_coef + n_class + 1) %*% solve(hes)
    }
    fit$vcov_sand <- vcov_sand
    fit$se_para <- sqrt(diag(fit$vcov_sand))
    names(fit$se_para) <- names(fit$est_para)
  }
  class(fit) <- "clus_lme"
  return(fit)
}

## ---- The function print.clus_lme ----
#' @title Print summary results of an clus_lme object
#'
#' @description \code{print.clus_lme} displays results of the output from \code{\link{clus_lme}}.
#'
#' @method print clus_lme
#' @param x an object of class "clus_lme", a result of \code{\link{clus_lme}} call.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param call logical. If set to \code{TRUE}, the matched call will be printed.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.clus_lme} shows a summary table for the estimated parameters in the cluster-effect model (continuous diagnostic test in three-class setting).
#'
#' @return \code{print.clus_lme} returns a summary table for the estimated parameters in the cluster-effect model.
#'
#' @seealso \code{\link{clus_lme}}
#'
#' @export
print.clus_lme <- function(x, digits = max(3L, getOption("digits") - 3L),
                           call = TRUE, ...) {
  cat("\n")
  if (call) {
    cat("CALL: ",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  }
  if (!is.null(x$se_para)) {
    z <- (x$est_para[1:x$n_coef] - rep(0, x$n_coef)) / x$se_para[1:x$n_coef]
    p_val <- 2 * pnorm(abs(z), lower.tail = FALSE)
    infer_tab <- cbind(c(x$est_para, x$icc),
                       c(x$se_para, rep(NA, 1)),
                       c(z, rep(NA, length(x$est_para) - x$n_coef + 1)),
                       c(p_val, rep(NA, length(x$est_para) - x$n_coef + 1)))
    rownames(infer_tab) <- c(names(x$est_para), "ICC")
    colnames(infer_tab) <- c("Est.", "Std.Error", "z-value", "p-value")
    cat("Coefficients:\n")
    printCoefmat(x = infer_tab, has.Pvalue = TRUE, digits = digits,
                 na.print = "--", ...)
  } else {
    infer_tab <- cbind(c(x$est_para, x$icc))
    rownames(infer_tab) <- c(names(x$est_para), "ICC")
    colnames(infer_tab) <- c("Est.")
    cat("Coefficients:\n")
    printCoefmat(x = infer_tab, has.Pvalue = FALSE, digits = digits, ...)
  }
  cat("\n")
  cat("Number of clusters:", x$cls, "\n")
  cat("Sample size within cluster:\n")
  print(c(Min = min(x$n_c), Max = max(x$n_c), Average = mean(x$n_c)))
  cat("Box-Cox transformation:", x$boxcox, "\n")
  cat("\n")
  invisible(x)
}

## ---- The function plot.clus_lme ----
#' @title Plot an clus_lme object.
#'
#' @description Diagnostic plots for the linear mixed-effect model, fitted by clus_lme.
#'
#' @method plot clus_lme
#'
#' @param x an object of class "clus_lme", i.e., a result of \code{\link{clus_lme}} call.
#' @param file_name File name to create on disk.
#' @param ... further arguments used with \code{\link{ggexport}} function, for example, \code{width}, \code{height}.
#'
#' @details \code{plot.clus_lme} provides three diagnostic plots: Q-Q plots for residuals, Fitted vs. Residuals values, and Q-Q plot for cluster effects, based on \code{ggplot()}.
#'
#' @return \code{plot.clus_lme} returns the diagnostic plots for the linear mixed-effect model, fitted by clus_lme.
#'
#' @seealso \code{\link{clus_lme}}
#'
#' @import ggplot2
#' @import ggpubr
#' @export
plot.clus_lme <- function(x, file_name = NULL, ...) {
  levl <- paste("Class:", names(x$residual))
  df_fitted_resid <- data.frame(fitted = unlist(x$fitted, use.names = FALSE),
                                resid = unlist(x$residual, use.names = FALSE),
                                Groups = factor(rep(levl,
                                                    sapply(x$residual, length)),
                                                levels = levl))
  p1 <- ggplot(df_fitted_resid, aes(sample = resid)) +
    stat_qq() +
    stat_qq_line() +
    facet_grid(~ Groups) +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles") +
    labs(subtitle = "Q-Q plots of residuals")
  p2 <- ggplot(df_fitted_resid, aes(x = fitted, y = resid)) +
    geom_point() +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    facet_grid(~ Groups) +
    xlab("Fitted") +
    ylab("Residuals") +
    labs(subtitle = "Fitted vs. Residuals")
  data_plot <- data.frame(randf = x$randf,
                          name = rep("Cluster effects", length(x$randf)))
  p3 <- ggplot(data_plot, aes_string(sample = "randf")) +
    stat_qq() +
    stat_qq_line() +
    facet_grid(~ name) +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles") +
    labs(subtitle = "Q-Q plot of Cluster effects")
  res <- ggarrange(ggarrange(p1, p2, ncol = 2, nrow = 1), p3, nrow = 2)
  if (!is.null(file_name)) {
    ggexport(res, filename = file_name, ...)
  }
  res
}
