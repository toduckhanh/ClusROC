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

#' @title Fitting the cluster-effect models for two-class of three-class settings.
#'
#' @description \code{lme2} fits the cluster-effect models for two-class or three-class setting based on the \code{lme()} routine from \code{nlme}-package.
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
#' \item{n_coef}{total numbers of coefficients included in the model.}
#' \item{n_p}{total numbers of regressors in the model.}
#' \item{icc}{a estimate of intra-class correlation - ICC}
#' \item{boxcox}{logical value indicating whether the Box-Cox transformation was implemented or not.}
#'
#' @examples
#' ## Example for two-class setting
#' data(data_2class)
#' head(data_2class)
#' out1 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
#'              data = data_2class, apVar = FALSE)
#' print(out1)
#' plot(out1)
#' out2 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
#'              data = data_2class)
#' print(out2, digits = 3)
#'
#' ## Example for three-class setting
#' data(data_3class)
#' head(data_3class)
#' out3 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
#'              data = data_3class, apVar = FALSE)
#' print(out3)
#' out4 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
#'              data = data_3class)
#' print(out4)
#' plot(out4)
#'
#' @export
lme2 <- function(name.test, name.class, name.covars, name.clust, data, levl.class = NULL,
                 boxcox = FALSE, apVar = TRUE, interval_lambda = c(-2, 2), trace = TRUE, ...){
  form.mean <- as.formula(paste(name.test, "~", name.class))
  mean.temp <- aggregate(form.mean, FUN = mean, data = data)
  temp.levl <- mean.temp[order(mean.temp[,2]), 1]
  if(is.null(levl.class)){
    if(trace){
      cat("The ordered levels of classes were not assigned by user!\n")
      cat("The ordered levels of classes are now determined by the orders of averages of diagnostic tests:\n")
      cat(paste(temp.levl, collapse = " < "), "\n")
    }
    levl.class <- temp.levl
  } else{
    if(all(levl.class == temp.levl)){
      if(trace){
        cat("The orders of inputed levels of classes are the same as the one obtained by the orders of averages of diagnostic tests:\n")
        cat(paste(levl.class, collapse = " < "),"\n")
      }
    } else{
      if(trace){
        cat("The orders of inputed levels of classes are not the same as the one obtained by the orders of averages of diagnostic tests:\n")
        cat("The correct one should be:\n")
        cat(paste(temp.levl, collapse = " < "),"\n")
      }
      levl.class <- temp.levl
    }
  }
  data[, name.class] <- factor(data[, name.class], levels = levl.class)
  ## define the formulas for fixed, random and weights
  call <- match.call()
  fit <- list()
  fit$call <- call
  fit$boxcox <- boxcox
  fit$name.test <- name.test
  fit$name.class <- name.class
  fit$name.covars <- name.covars
  fit$name.clust <- name.clust
  if(!missing(name.covars)){
    fixed <- as.formula(paste(name.test, "~", name.class, "+",
                              paste0(name.covars, ":", name.class, collapse = " + "), "- 1"))
  } else{
    fixed <- as.formula(paste(name.test, "~", name.class, "- 1"))
  }
  random <- as.formula(paste("~", "1|", name.clust))
  form.weights <- as.formula(paste("~", "1|", name.class))
  weights <- varIdent(form = form.weights)
  n <- nrow(data)
  Clus <- model.frame(getGroupsFormula(random), data = data)[,1]
  n_c <- table(Clus)
  cls <- length(unique(Clus))
  if(isFALSE(boxcox)){
    out_model <- lme(fixed = fixed, random = random, weights = weights, method = "REML", data = data)
  } else{
    all.Y <- model.response(model.frame(fixed, data = data))
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
  n_class <- length(table(data[, name.class]))
  n_p <- n_coef/n_class
  if(n_class == 2){
    id_coef <- c(seq(1, n_coef - 1, by = 2), seq(2, n_coef, by = 2))
  }
  if(n_class == 3){
    id_coef <- c(seq(1, n_coef - 2, by = 3), seq(2, n_coef - 1, by = 3), seq(3, n_coef, by = 3))
  }
  sigma_e_est <- coef(out_model$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)*out_model$sigma
  sigma_e_est <- sigma_e_est[levl.class]
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
#' @title Print summary results of lme2
#'
#' @description \code{print.lme2} prints the results for the output of function \code{\link{lme2}}.
#'
#' @method print lme2
#' @param x an object of class "lme2", a result of a call to \code{\link{lme2}}.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.lme2} shows a nice format of the summary table for fitting the cluster-effect models.
#'
#' @seealso \code{\link{lme2}}
#'
#' @export
print.lme2 <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\n")
  cat("CALL: ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  if(!is.null(x$se_para)){
    z <- (x$est_para[1:x$n_coef] - rep(0, x$n_coef))/x$se_para[1:x$n_coef]
    p_val <- 2*pnorm(abs(z), lower.tail = FALSE)
    infer_tab <- cbind(c(x$est_para, x$icc),
                       c(x$se_para[1:x$n_coef], rep(NA, length(x$est_para) - x$n_coef + 1)),
                       c(z, rep(NA, length(x$est_para) - x$n_coef + 1)),
                       c(p_val, rep(NA, length(x$est_para) - x$n_coef + 1)))
    rownames(infer_tab) <- c(names(x$est_para), "ICC")
    colnames(infer_tab) <- c("Est.", "Std.Error", "z-value", "p-value")
    printCoefmat(infer_tab, has.Pvalue = TRUE, digits = digits, na.print = "--")
  }
  else{
    infer_tab <- cbind(c(x$est_para, x$icc))
    rownames(infer_tab) <- c(names(x$est_para), "ICC")
    colnames(infer_tab) <- c("Est.")
    printCoefmat(infer_tab, has.Pvalue = FALSE, digits = digits)
  }
  invisible(x)
}

## ---- The function plot.lme2 ----
#' @title Diagnostic plots for model fitted by lme2
#'
#' @description \code{plot.lme2} provides diagnostic plots (based on \code{ggplot()}) for the residuals and cluster effects.
#'
#' @method plot lme2
#' @param x an object of class "lme2", a result of a call to \code{\link{lme2}}.
#' @param file.name File name to create on disk.
#'
#' @details \code{plot.lme2} shows three diagnostic plots: Q-Q plots of residuals, Fitted vs. Residuals, and Q-Q plot of cluster effects.
#'
#' @seealso \code{\link{lme2}}
#'
#' @import ggplot2
#' @import ggpubr
#' @export
plot.lme2 <- function(x, file.name = NULL){
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
  p3 <- ggplot(data = data.frame(randf = x$randf, name = rep("Cluster effects", length(x$randf))),
               aes(sample = randf)) +
    stat_qq() + stat_qq_line() +
    facet_grid(~ name) +
    xlab("Theoretical Quantiles") + ylab("Sample Quantiles") +
    labs(subtitle = "Q-Q plot of Cluster effects")
  res <- ggarrange(ggarrange(p1, p2, ncol = 2, nrow = 1), p3, nrow = 2)
  if(!is.null(file.name)){
    ggexport(res, filename = file.name)
  }
  res
}


