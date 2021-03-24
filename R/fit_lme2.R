####==========================================================================####
## This file consists of functions for fitting the cluster-effect models        ##
## REML approach, based on the lme() routine                                    ##
## Date: 23/03/2021																															##
####==========================================================================####

#' @import numDeriv
#' @import nlme
#' @import stats
#' @import utils

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

#' @title Fitting the cluster-effect models for two-class of three-class settings.
#'
#' @description \code{lme2} the cluster-effect models for two-class or three-class setting based on the \code{lme()} routine from \code{nlme}-package.
#'
#' @param name.test  name of variable indicating diagnostic test of biomarker in data. It can be also the transformation, for example, \code{"log(y)"}, where the term \code{log} is for the log-transformation, and \code{y} is the name of test.
#' @param name.class  name of variable indicating disease classes (diagnostic groups) in data.
#' @param name.covars  an vector of names of covariates containing in data. The vector can contain also the transformation, for example, \code{c("x1", "sqrt(x2)", "I(x3^2)")}.
#' @param name.clust  name of variable indicating clusters in data.
#' @param data  a data frame containing the variables in the model.
#' @param levl.class  an vector of the unique values (as character strings) that (disease) class might have taken, sorted into increasing order of means of test results corresponding to the disease classes (diagnostic groups). If \code{levl.class = NULL}, the levels will be automatically determined based on data, and sorted into increasing order of means of test results corresponding to the disease classes (diagnostic groups).
#' @param apVar  a logical value. Default = \code{TRUE}. If set to \code{TRUE}, the covariance matrix for all estimated parameters in model with be obtained by using the sandwich formula.
#' @param ...  additional arguments for \code{\link{lme}}.
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
lme2 <- function(name.test, name.class, name.covars, name.clust, data, levl.class = NULL, apVar = TRUE, ...){
  form.mean <- as.formula(paste(name.test, "~", name.class))
  mean.temp <- aggregate(form.mean, FUN = mean, data = data)
  temp.levl <- mean.temp[order(mean.temp[,2]), 1]
  if(is.null(levl.class)){
    cat("The ordered levels of classes were not assigned by user!\n")
    cat("The ordered levels of classes are now determined by the orders of averages of diagnostic tests:\n")
    cat(paste(temp.levl, collapse = " < "), "\n")
    levl.class <- temp.levl
  } else{
    if(all(levl.class == temp.levl)){
      cat("The orders of inputed levels of classes are the same as the one obtained by the orders of averages of diagnostic tests:\n")
      cat(paste(levl.class, collapse = " < "),"\n")
    } else{
      cat("The orders of inputed levels of classes are not the same as the one obtained by the orders of averages of diagnostic tests:\n")
      cat("The correct one should be:\n")
      cat(paste(temp.levl, collapse = " < "),"\n")
      levl.class <- temp.levl
    }
  }
  data[, name.class] <- factor(data[, name.class], levels = levl.class)
  ## define the formulas for fixed, random and weights
  call <- match.call()
  fixed <- as.formula(paste(name.test, "~", name.class, "+",
                            paste0(name.covars, ":", name.class, collapse = " + "), "- 1"))
  random <- as.formula(paste("~", "1|", name.clust))
  form.weights <- as.formula(paste("~", "1|", name.class))
  weights <- varIdent(form = form.weights)
  ## call lme() with REML
  out_lme <- lme(fixed = fixed, random = random, weights = weights, method = "REML", data = data)
  ## collecting results
  n_coef <- length(out_lme$coefficients$fixed)
  n_class <- length(table(data[, name.class]))
  n_p <- n_coef/n_class
  if(n_class == 2){
    id_coef <- c(seq(1, n_coef - 1, by = 2), seq(2, n_coef, by = 2))
  }
  if(n_class == 3){
    id_coef <- c(seq(1, n_coef - 2, by = 3), seq(2, n_coef - 1, by = 3), seq(3, n_coef, by = 3))
  }
  sigma_e_est <- coef(out_lme$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)*out_lme$sigma
  sigma_e_est <- sigma_e_est[levl.class]
  sigma_c_est <- sqrt(as.numeric(getVarCov(out_lme)))
  par_lme_est <- c(out_lme$coefficients$fixed[id_coef], sigma_c_est, sigma_e_est)
  icc <- sigma_c_est^2/(sigma_c_est^2 + mean(sigma_e_est)^2)
  fit <- list()
  fit$call <- call
  fit$n_coef <- n_coef
  fit$est_para <- par_lme_est
  names(fit$est_para) <- c(names(par_lme_est)[1:n_coef], "sigma_c", paste0("sigma_", c(1:n_class)))
  fit$icc <- icc
  ##
  resid <- out_lme$residual[,2]
  fit$residual <- split(resid, data[, name.class])
  fitted <- out_lme$fitted[,2]
  fit$fitted <- split(fitted, data[, name.class])
  fit$randf <- ranef(out_lme)$`(Intercept)`
  if(apVar){
    ## preparing input for sandwich covariance matrix
    Clus <- out_lme$groups[,1]
    cls <- length(unique(Clus))
    n_c <- table(Clus)
    all.Y <- model.response(model.frame(out_lme$terms, data = data))
    data_matrix <- model.matrix(out_lme$terms, data = data)
    all.D <- data_matrix[, 1:n_class]
    all.Z <- data_matrix[, id_coef]
    Y <- split(all.Y, Clus)
    D <- lapply(split(as.data.frame(all.D), Clus), as.matrix)
    Z <- lapply(split(as.data.frame(all.Z), Clus), as.matrix)
    V <- lapply(n_c, function(x) rep(1, x))
    jac <- jacobian(reml_loglik_vec, x = par_lme_est, D = D, Y = Y, Z = Z, V = V,
                    cls = cls, n_p = n_p, n_class = n_class)
    hes <- hessian(reml_loglik_vec_sum, x = par_lme_est, D = D, Y = Y, Z = Z, V = V,
                   cls = cls, n_p = n_p, n_class = n_class)
    vcov_sand <- solve(hes) %*% matrix(rowSums(apply(jac, 1, tcrossprod)), n_coef + 1 + n_class,
                                       n_coef + 1 + n_class) %*% solve(hes)
    se_lme_est <- sqrt(diag(vcov_sand))
    fit$se_para <- se_lme_est
    names(fit$se_para) <- c(names(par_lme_est)[1:n_coef], "sigma_c", paste0("sigma_", c(1:n_class)))
    fit$vcov_sand <- vcov_sand
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
    colnames(infer_tab) <- c("Est.", "Std.Error", "z-value", "p-value")
    rownames(infer_tab) <- c(names(x$est_para), "ICC")
    printCoefmat(infer_tab, has.Pvalue = TRUE, digits = digits, na.print = "--")
  }
  else{
    infer_tab <- cbind(c(x$est_para, x$icc))
    colnames(infer_tab) <- c("Est.")
    rownames(infer_tab) <- c(names(x$est_para), "ICC")
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


