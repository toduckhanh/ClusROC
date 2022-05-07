####==========================================================================####
## This file consists of functions for fitting the cluster-effect models        ##
## REML approach, based on the lme() routine                                    ##
## Date: 23/03/2021																															##
####==========================================================================####
#' @import numDeriv
#' @import nlme
#' @import stats
#' @import utils

### ---- reml log-likelihood for normal distribution ----
reml_loglik_vec <- function(par, D, Y, Z, V, cls, n_p, n_class){
  beta_fit <- par[1:(n_class*n_p)]
  sigma_c <- par[(n_class*n_p + 1)]
  sigma_e <- par[(n_class*n_p + 2):length(par)]
  Sigma <- lapply(1:cls, function(i){
    tem <- tcrossprod(V[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(D[[i]] %*% sigma_e^2), ncol = ncol(tem), nrow = ncol(tem)))
  })
  Sig_chol <- lapply(Sigma, chol)
  Sig_chol_inv <- lapply(Sig_chol, FUN = function(x) backsolve(x, diag(1, ncol(x))))
  Sig_chol_inv_Z <- mapply(FUN = function(x, y) crossprod(x, y), x = Sig_chol_inv, y = Z, SIMPLIFY = FALSE)
  term_1 <- numeric(cls)
  for(i in 1:cls){
    Sig_chol_inv_y <- crossprod(Sig_chol_inv[[i]], Y[[i]])
    term_1[i] <- 0.5*crossprod(Sig_chol_inv_y - Sig_chol_inv_Z[[i]] %*% beta_fit) +
      sum(log(diag(Sig_chol[[i]])))
  }
  tem <- Reduce("+", lapply(Sig_chol_inv_Z, FUN = crossprod))
  term_2 <- sum(log(diag(chol(tem))))/cls
  return(term_1 + term_2/2)
}

reml_loglik_vec_sum <- function(par, D, Y, Z, V, cls, n_p, n_class){
  beta_fit <- par[1:(n_class*n_p)]
  sigma_c <- par[(n_class*n_p + 1)]
  sigma_e <- par[(n_class*n_p + 2):length(par)]
  Sigma <- lapply(1:cls, function(i){
    tem <- tcrossprod(V[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(D[[i]] %*% sigma_e^2), ncol = ncol(tem), nrow = ncol(tem)))
  })
  Sig_chol <- lapply(Sigma, chol)
  Sig_chol_inv <- lapply(Sig_chol, FUN = function(x) backsolve(x, diag(1, ncol(x))))
  Sig_chol_inv_Z <- mapply(FUN = function(x,y) crossprod(x, y), x = Sig_chol_inv, y = Z, SIMPLIFY = FALSE)
  term_1 <- numeric(cls)
  for(i in 1:cls){
    Sig_chol_inv_y <- crossprod(Sig_chol_inv[[i]], Y[[i]])
    term_1[i] <- 0.5*crossprod(Sig_chol_inv_y - Sig_chol_inv_Z[[i]] %*% beta_fit) +
      sum(log(diag(Sig_chol[[i]])))
  }
  tem <- Reduce("+", lapply(Sig_chol_inv_Z, FUN = crossprod))
  term_2 <- sum(log(diag(chol(tem))))
  return(sum(term_1) + term_2/2)
}

### ---- reml log-likelihood for Box-Cox transformation ----
reml_bcx_loglik_vec <- function(par, D, Y, Z, V, cls, n_p, n_class){
  beta_fit <- par[1:(n_class*n_p)]
  sigma_c <- par[(n_class*n_p + 1)]
  sigma_e <- par[(n_class*n_p + 2):(length(par) - 1)]
  lambda <- par[length(par)]
  y_boxcox <- lapply(Y, function(x) boxcox_trans(x, lambda))
  Sigma <- lapply(1:cls, function(i){
    tem <- tcrossprod(V[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(D[[i]] %*% sigma_e^2), ncol = ncol(tem), nrow = ncol(tem)))
  })
  Sig_chol <- lapply(Sigma, chol)
  Sig_chol_inv <- lapply(Sig_chol, FUN = function(x) backsolve(x, diag(1, ncol(x))))
  Sig_chol_inv_Z <- mapply(FUN = function(x,y) crossprod(x, y), x = Sig_chol_inv, y = Z, SIMPLIFY = FALSE)
  term_1 <- numeric(cls)
  term_3 <- numeric(cls)
  for(i in 1:cls){
    Sig_chol_inv_y <- crossprod(Sig_chol_inv[[i]], y_boxcox[[i]])
    term_1[i] <- 0.5*crossprod(Sig_chol_inv_y - Sig_chol_inv_Z[[i]] %*% beta_fit) +
      sum(log(diag(Sig_chol[[i]])))
    term_3[i] <- sum(log(Y[[i]]))
  }
  tem <- Reduce("+", lapply(Sig_chol_inv_Z, FUN = crossprod))
  term_2 <- sum(log(diag(chol(tem))))/cls
  return(term_1 + term_2/2 - (lambda - 1)*term_3)
}

reml_bcx_loglik_vec_sum <- function(par, D, Y, Z, V, cls, n_p, n_class){
  beta_fit <- par[1:(n_class*n_p)]
  sigma_c <- par[(n_class*n_p + 1)]
  sigma_e <- par[(n_class*n_p + 2):(length(par) - 1)]
  lambda <- par[length(par)]
  y_boxcox <- lapply(Y, function(x) boxcox_trans(x, lambda))
  Sigma <- lapply(1:cls, function(i){
    tem <- tcrossprod(V[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(D[[i]] %*% sigma_e^2), ncol = ncol(tem), nrow = ncol(tem)))
  })
  Sig_chol <- lapply(Sigma, chol)
  Sig_chol_inv <- lapply(Sig_chol, FUN = function(x) backsolve(x, diag(1, ncol(x))))
  Sig_chol_inv_Z <- mapply(FUN = function(x,y) crossprod(x, y), x = Sig_chol_inv, y = Z, SIMPLIFY = FALSE)
  term_1 <- numeric(cls)
  term_3 <- numeric(cls)
  for(i in 1:cls){
    Sig_chol_inv_y <- crossprod(Sig_chol_inv[[i]], y_boxcox[[i]])
    term_1[i] <- 0.5*crossprod(Sig_chol_inv_y - Sig_chol_inv_Z[[i]] %*% beta_fit) +
      sum(log(diag(Sig_chol[[i]])))
    term_3[i] <- sum(log(Y[[i]]))
  }
  tem <- Reduce("+", lapply(Sig_chol_inv_Z, FUN = crossprod))
  term_2 <- sum(log(diag(chol(tem))))
  return(sum(term_1) + term_2/2 - (lambda - 1)*sum(term_3))
}

llike_bcx_fun <- function(par, fixed, random, weights, data, y_tit, ...){
  y <- model.response(model.frame(fixed, data = data))
  y_boxcox <- boxcox_trans(y, par)
  W <- y_boxcox/(y_tit^(par - 1))
  fixed_new <- update(fixed, W  ~ . )
  data$W <- W
  out_model <- lme(fixed = fixed_new, random = random, weights = weights, method = "REML", data = data, ...)
  return(as.numeric(out_model$logLik))
}

### ---- Fitting cluster-effect models ----

#' @title Linear Mixed-Effects Models for a continuous diagnostic test.
#'
#' @description \code{lme2} fits the cluster-effect model for a continuous diagnostic test in a three-class setting as described in Xiong et al. (2018) and To et al. (2022).
#'
#' @param fixed.formula  a two-sided linear formula object, describing the fixed-effects part of the model for three classes, with the response on the left of ~ operator and the terms, separated by + operators, on the right. For example, \code{Y ~ X1 + X2}, \code{Y ~ X1 + X2 + X1:X2} or \code{log(Y) ~ X1 + X2 + I(X1^2)}.
#' @param name.class  name of variable indicating disease classes (or diagnostic groups) in the data.
#' @param name.clust  name of variable indicating clusters in the data.
#' @param data  a data frame containing the variables in the model.
#' @param levl.class  a vector (of strings) containing the ordered name chosen for the disease classes. The ordering is intended to be ``increasing'' with respect to the disease severity. If \code{levl.class = NULL} (default), the elements of the vector will be automatically determined from data, by considering the order of the means of the test values for each disease class (diagnostic group).
#' @param apVar  a logical value. Default = \code{TRUE}. If set to \code{TRUE}, the estimated covariance matrix for all estimated parameters in the model will be obtained (by using the sandwich formula).
#' @param boxcox  a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a Box-Cox transformation will be applied to the model.
#' @param interval_lambda  a vector containing the end-points of the interval for searching the Box-Cox parameter, \code{lambda}. Default = (-2, 2).
#' @param trace  a logical value. Default = \code{TRUE}. If set to \code{TRUE}, the information about the check for the monotonic ordering of test values will be provided.
#' @param ...  additional arguments for \code{\link[nlme]{lme}}, such as \code{control}, \code{contrasts}.
#'
#' @details
#' This function fits a linear mixed-effect model for a continuous diagnostic test in a three-class setting in order to account for the cluster and covariates effects on the test result. See Xiong et al. (2018) and To et al. (2022) for more details.
#' \itemize{
#' \item Estimation is done by using \code{\link[nlme]{lme}} with the restricted maximum log-likelihood (REML) method.
#' \item Box-Cox transformation for the model can be used when the distributions of test results are skewed (Gurka et al. 2006). The estimation procedure is described in To et al. (2022). The Box-Cox parameter \eqn{\lambda} is estimated by a grid search on the interval [-2, 2], as discussed in Gurka and Edwards (2011).
#' \item The estimated variance-covariance matrix for the estimated parameters are obtained by sandwich formula (see, Liang and Zeger, 1986; Kauermann and Carroll, 2001; Mancl and DeRouen, 2001) as discussed in To et al. (2022).
#' }
#'
#'
#' @return \code{lme2} returns an object of class "lme2" class, i.e., a list containing at least the following components:
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
#' ``Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data''. \emph{Statistical Methods in Medical Research}, DOI: 10.1177/09622802221089029.
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
#' out1 <- lme2(fixed.formula = Y ~ X1 + X2, name.class = "D", name.clust = "id_Clus",
#'              data = data_3class)
#' print(out1)
#' plot(out1)
#'
#' ## Example 2: Box-Cox transformation
#' data(data_3class_bcx)
#' out2 <- lme2(fixed.formula = Y ~ X, name.class = "D", name.clust = "id_Clus",
#'              data = data_3class_bcx, boxcox = TRUE)
#' print(out2)
#' plot(out2)
#'
#'
#' @export
lme2 <- function(fixed.formula, name.class, name.clust, data, levl.class = NULL, apVar = TRUE,
                 boxcox = FALSE, interval_lambda = c(-2, 2), trace = TRUE, ...){
  if(missing(data)){
    data <- .GlobalEnv
    cat("Warning: the data is missing, the global environment is used!\n")
  }
  if(missing(fixed.formula)) stop("argument fixed.formula is missing with no default")
  if(!inherits(fixed.formula, "formula") || length(fixed.formula) != 3) {
    stop("\nfixed-effects model must be a formula of the form \"resp ~ pred\"")
  }
  if(missing(name.class)) stop("argument name.class is missing with no default")
  if(!inherits(name.class, "character") || length(name.class) != 1)
    stop("agrument name.class must be a character vector with length 1.")
  if(missing(name.clust)) stop("argument name.clust is missing with no default")
  if(!inherits(name.clust, "character") || length(name.clust) != 1)
    stop("agrument name.clust must be a character vector with length 1.")
  n_class <- length(table(data[, name.class]))
  if(n_class != 3) stop("There is not a case of three-class setting!")
  mf <- unlist(strsplit(as.character(fixed.formula), "~"))[-1]
  name.test <- mf[1]
  if(boxcox){
    if(any(data[, name.test] < 0)) stop("Cannot apply Box-Cox transform for negative values.")
  }
  form.mean <- as.formula(paste(name.test, "~", name.class))
  mean.temp <- aggregate(form.mean, FUN = mean, data = data)
  temp.levl <- mean.temp[order(mean.temp[,2]), 1]
  if(is.null(levl.class)){
    if(trace){
      cat("The ordered levels of classes were not assigned by user!\n")
      cat("The ordered levels of classes are now determined by the orders of averages of tests results:\n")
      cat(paste(temp.levl, collapse = " < "), "\n")
    }
    levl.class <- temp.levl
  } else{
    if(!inherits(levl.class, "character") || length(levl.class) != 3)
      stop("agrument levl.class must be a character vector with length 3.")
    if(all(levl.class == temp.levl)){
      if(trace){
        cat("The orders of inputed levels of classes are the same as the one obtained by the orders of averages of tests results:\n")
        cat(paste(levl.class, collapse = " < "),"\n")
      }
    } else{
      if(trace){
        cat("The orders of inputed levels of classes are not the same as the one obtained by the orders of averages of tests results:\n")
        cat("The correct one should be:\n")
        cat(paste(temp.levl, collapse = " < "),"\n")
      }
      levl.class <- temp.levl
    }
  }
  data[, name.class] <- factor(data[, name.class], levels = levl.class)
  name.covars <- mf[2]
  ## define the formulas for fixed, random and weights
  call <- match.call()
  fit <- list()
  fit$call <- call
  fit$boxcox <- boxcox
  fit$name.test <- name.test
  fit$name.class <- name.class
  fit$name.clust <- name.clust
  fixed <- as.formula(paste(name.test, "~", name.class, "+", "(", name.covars, ")", ":", name.class, "-1"))
                            # paste0(name.covars, ":", name.class, collapse = " + "), "- 1"))
  fit$name.covars <- name.covars
  random <- as.formula(paste("~", "1|", name.clust))
  form.weights <- as.formula(paste("~", "1|", name.class))
  weights <- varIdent(form = form.weights)
  fit$terms <- terms(fixed)
  attr(fit$terms, "levl.class") <- levl.class
  attr(fit$terms, "n_vb") <- length(as.character(attr(fit$terms, "variables"))[-c(1:3)])
  n <- nrow(data)
  Clus <- model.frame(getGroupsFormula(random), data = data)[,1]
  n_c <- table(Clus)
  cls <- length(unique(Clus))
  if(isFALSE(boxcox)){
    out_model <- lme(fixed = fixed, random = random, weights = weights, method = "REML", data = data)
  } else{
    all.Y <- model.response(model.frame(fixed, data = data))
    ## check ok!
    list.Y <- split(all.Y, Clus)
    y_tit <- prod(sapply(list.Y, function(x) prod(x^(1/n))))
    lambda_est <- optimize(llike_bcx_fun, interval = interval_lambda, fixed = fixed, random = random,
                           weights = weights, data = data, y_tit = y_tit, maximum = TRUE, ...)$maximum
    all.Y_boxcox <- boxcox_trans(all.Y, lambda_est)
    data$Y_boxcox <- all.Y_boxcox
    fixed_new <- update(fixed, Y_boxcox  ~ . )
    out_model <- lme(fixed = fixed_new, random = random, weights = weights, method = "REML", data = data, ...)
  }
  ## collecting results
  n_coef <- length(out_model$coefficients$fixed)
  n_p <- n_coef/n_class
  if(n_class == 2){
    id_coef <- c(seq(1, n_coef - 1, by = 2), seq(2, n_coef, by = 2))
  }
  if(n_class == 3){
    id_coef <- c(seq(1, n_coef - 2, by = 3), seq(2, n_coef - 1, by = 3), seq(3, n_coef, by = 3))
  }
  sigma_e_est <- coef(out_model$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)*out_model$sigma
  sigma_e_est <- sigma_e_est[as.character(levl.class)]
  sigma_c_est <- sqrt(as.numeric(getVarCov(out_model)))
  if(boxcox){
    par_est <- c(out_model$coefficients$fixed[id_coef], sigma_c_est, sigma_e_est, lambda_est)
    names(par_est) <- c(names(par_est)[1:n_coef], "sigma_c", paste0("sigma_", c(1:n_class)), "lambda")
    #fit$lambda <- lambda_est
  } else{
    par_est <- c(out_model$coefficients$fixed[id_coef], sigma_c_est, sigma_e_est)
    names(par_est) <- c(names(par_est)[1:n_coef], "sigma_c", paste0("sigma_", c(1:n_class)))
  }
  icc <- sigma_c_est^2/(sigma_c_est^2 + mean(sigma_e_est)^2)
  fit$n_coef <- n_coef
  fit$n_p <- n_p
  fit$n <- n
  fit$cls <- cls
  fit$n_c <- n_c
  fit$est_para <- par_est #[1:(n_class*n_p + 4)]
  fit$icc <- icc
  ##
  resid <- out_model$residual[,2]
  fit$residual <- split(resid, data[, name.class])
  fitted <- out_model$fitted[,2]
  fit$fitted <- split(fitted, data[, name.class])
  fit$randf <- ranef(out_model)$`(Intercept)`
  if(apVar){
    if(isFALSE(boxcox)){
      all.Y <- model.response(model.frame(out_model$terms, data = data))
      data_matrix <- model.matrix(out_model$terms, data = data)
      all.D <- data_matrix[, 1:n_class]
      all.Z <- data_matrix[, id_coef]
      Y <- split(all.Y, Clus)
      D <- lapply(split(as.data.frame(all.D), Clus), as.matrix)
      Z <- lapply(split(as.data.frame(all.Z), Clus), as.matrix)
      V <- lapply(n_c, function(x) rep(1, x))
      jac <- jacobian(reml_loglik_vec, x = par_est, D = D, Y = Y, Z = Z, V = V,
                      cls = cls, n_p = n_p, n_class = n_class)
      hes <- hessian(reml_loglik_vec_sum, x = par_est, D = D, Y = Y, Z = Z, V = V,
                     cls = cls, n_p = n_p, n_class = n_class)
      vcov_sand <- solve(hes) %*% matrix(rowSums(apply(jac, 1, tcrossprod)), n_coef + n_class + 1,
                                         n_coef + n_class + 1) %*% solve(hes)
    } else{
      data_matrix <- model.matrix(out_model$terms, data = data)
      all.D <- data_matrix[, 1:n_class]
      all.Z <- data_matrix[, id_coef]
      Y <- split(all.Y, Clus)
      D <- lapply(split(as.data.frame(all.D), Clus), as.matrix)
      Z <- lapply(split(as.data.frame(all.Z), Clus), as.matrix)
      V <- lapply(n_c, function(x) rep(1, x))
      jac <- jacobian(reml_bcx_loglik_vec, x = par_est, D = D, Y = Y, Z = Z, V = V, cls = cls, n_p = n_p,
                      n_class = n_class)
      hes <- hessian(reml_bcx_loglik_vec_sum, x = par_est, D = D, Y = Y, Z = Z, V = V, cls = cls,
                     n_p = n_p, n_class = n_class)
      vcov_sand <- solve(hes) %*% matrix(rowSums(apply(jac, 1, tcrossprod)), n_coef + n_class + 2,
                                         n_coef + n_class + 2) %*% solve(hes)
    }
    fit$vcov_sand <- vcov_sand
    fit$se_para <- sqrt(diag(fit$vcov_sand))
    names(fit$se_para) <- names(fit$est_para)
  }
  class(fit) <- "lme2"
  return(fit)
}

## ---- The function print.lme2 ----
#' @title Print summary results of an lme2 object
#'
#' @description \code{print.lme2} displays results of the output from \code{\link{lme2}}.
#'
#' @method print lme2
#' @param x an object of class "lme2", a result of \code{\link{lme2}} call.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param call logical. If \code{TRUE}, the matched call will be printed.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.lme2} shows a summary table for the estimated parameters in the cluster-effect model (continuous diagnostic test in three-class setting).
#'
#' @seealso \code{\link{lme2}}
#'
#' @export
print.lme2 <- function(x, digits = max(3L, getOption("digits") - 3L), call = TRUE, ...){
  cat("\n")
  if(call){
    cat("CALL: ",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  }
  if(!is.null(x$se_para)){
    z <- (x$est_para[1:x$n_coef] - rep(0, x$n_coef))/x$se_para[1:x$n_coef]
    p_val <- 2*pnorm(abs(z), lower.tail = FALSE)
    infer_tab <- cbind(c(x$est_para, x$icc),
                       c(x$se_para[1:x$n_coef], rep(NA, length(x$est_para) - x$n_coef + 1)),
                       c(z, rep(NA, length(x$est_para) - x$n_coef + 1)),
                       c(p_val, rep(NA, length(x$est_para) - x$n_coef + 1)))
    rownames(infer_tab) <- c(names(x$est_para), "ICC")
    colnames(infer_tab) <- c("Est.", "Std.Error", "z-value", "p-value")
    cat("Coefficients:\n")
    printCoefmat(infer_tab, has.Pvalue = TRUE, digits = digits, na.print = "--")
  }
  else{
    infer_tab <- cbind(c(x$est_para, x$icc))
    rownames(infer_tab) <- c(names(x$est_para), "ICC")
    colnames(infer_tab) <- c("Est.")
    cat("Coefficients:\n")
    printCoefmat(infer_tab, has.Pvalue = FALSE, digits = digits)
  }
  cat("\n")
  cat("Number of clusters:", x$cls, "\n")
  cat("Sample size within cluster:\n")
  print(c(Min = min(x$n_c), Max = max(x$n_c), Average = mean(x$n_c)))
  cat("Box-Cox transformation:", x$boxcox, "\n")
  cat("\n")
  invisible(x)
}

## ---- The function plot.lme2 ----
#' @title Plot an lme2 object.
#'
#' @description Diagnostic plots for the linear mixed-effect model, fitted by lme2.
#'
#' @method plot lme2
#'
#' @param x an object of class "lme2", i.e., a result of \code{\link{lme2}} call.
#' @param file.name File name to create on disk.
#' @param ... further arguments used with \code{\link{ggexport}} function, for example, \code{width}, \code{height}.
#'
#' @details \code{plot.lme2} provides three diagnostic plots: Q-Q plots for residuals, Fitted vs. Residuals values, and Q-Q plot for cluster effects, based on \code{ggplot()}.
#'
#' @seealso \code{\link{lme2}}
#'
#' @import ggplot2
#' @import ggpubr
#' @export
plot.lme2 <- function(x, file.name = NULL, ...){
  levl <- paste("Class:", names(x$residual))
  df_fitted_resid <- data.frame(fitted = unlist(x$fitted, use.names = FALSE),
                                resid = unlist(x$residual, use.names = FALSE),
                                Groups = factor(rep(levl, sapply(x$residual, length)), levels = levl))
  p1 <- ggplot(df_fitted_resid, aes(sample = resid)) + stat_qq() + stat_qq_line() +
    facet_grid( ~ Groups) + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") +
    labs(subtitle = "Q-Q plots of residuals")
  p2 <- ggplot(df_fitted_resid, aes(x = fitted, y = resid)) +
    geom_point() + geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    facet_grid( ~ Groups) + xlab("Fitted") + ylab("Residuals") +
    labs(subtitle = "Fitted vs. Residuals")
  data_plot <- data.frame(randf = x$randf, name = rep("Cluster effects", length(x$randf)))
  p3 <- ggplot(data_plot, aes_string(sample = "randf")) +
    stat_qq() + stat_qq_line() +
    facet_grid(~ name) +
    xlab("Theoretical Quantiles") + ylab("Sample Quantiles") +
    labs(subtitle = "Q-Q plot of Cluster effects")
  res <- ggarrange(ggarrange(p1, p2, ncol = 2, nrow = 1), p3, nrow = 2)
  if(!is.null(file.name)){
    ggexport(res, filename = file.name, ...)
  }
  res
}
