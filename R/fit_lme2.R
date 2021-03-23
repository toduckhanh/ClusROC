####==========================================================================####
## This file consists of functions for fitting the cluster-effect models        ##
## REML approach, based on the lme() routine                                    ##
## Date: 23/03/2021																															##
####==========================================================================####

#' @import numDeriv
#' @import nlme
#' @import stats
#' @import utils

reml_loglik_vec <- function(par, D, Y, Z, V, Clus, cls, n_p, n_class){
  beta_fit <- par[1:(n_class*n_p)]
  sigma_c <- par[(n_class*n_p + 1)]
  sigma_e <- par[(n_class*n_p + 2):length(par)]
  Sigma <- lapply(1:cls, function(i){
    tem <- tcrossprod(V[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(matrix(D[Clus == i,], ncol = n_class,
                                                    byrow = FALSE) %*% sigma_e^2),
                                  ncol = ncol(tem), nrow = ncol(tem)))
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

reml_loglik_vec_sum <- function(par, D, Y, Z, V, Clus, cls, n_p, n_class){
  beta_fit <- par[1:(n_class*n_p)]
  sigma_c <- par[(n_class*n_p + 1)]
  sigma_e <- par[(n_class*n_p + 2):length(par)]
  Sigma <- lapply(1:cls, function(i){
    tem <- tcrossprod(V[[i]])
    return(sigma_c^2 * tem + diag(as.numeric(matrix(D[Clus == i,], ncol = n_class,
                                                    byrow = FALSE) %*% sigma_e^2),
                                  ncol = ncol(tem), nrow = ncol(tem)))
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
#' @param name.test  ....
#' @param name.class  ....
#' @param name.covars  ...
#' @param name.clust  ...
#' @param data  ...
#' @param apVar  a logical value ...
#' @param residual  a logical value ...
#' @param fitted a logical value. Default = \code{FALSE}. If set to \code{TRUE}, ...
#' @param randf  a logical value ...
#' @param ...  additional arguments for \code{lme()}.
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
#' \item{est_para}{the estimate of VUS.}
#' \item{se_para}{the standard error, obtained by using asymptotic theory or bootstrap resampling method.}
#' \item{vcov_sand}{the matched call.}
#' \item{residual}{...}
#' \item{fitted}{...}
#' \item{randf}{...}
#' \item{n_coef}{...}
#'
#' @examples
#' ## Example for two-class setting
#' data(data_2class)
#' head(data_2class)
#' out1 <- lme2(name.test = "Y", name.class = "D", name.covars = c("X1", "X2"), name.clust = "id_Clus",
#'              data = data_2class, apVar = FALSE)
#' print(out1)
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
#'
#' @export
lme2 <- function(name.test, name.class, name.covars, name.clust, data, apVar = TRUE,
                 residual = FALSE, fitted = FALSE, randf = FALSE, ...){
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
  sigma_e_est <- sigma_e_est[order(names(sigma_e_est))]
  sigma_c_est <- sqrt(as.numeric(getVarCov(out_lme)))
  par_lme_est <- c(out_lme$coefficients$fixed[id_coef], sigma_c_est, sigma_e_est)
  fit <- list()
  fit$call <- call
  fit$n_coef <- n_coef
  fit$est_para <- par_lme_est
  names(fit$est_para) <- c(names(par_lme_est)[1:n_coef], "sigma_c", paste0("sigma_", c(1:n_class)))
  ## preparing input for sandwich covariance matrix
  Clus <- out_lme$groups[,1]
  cls <- length(unique(Clus))
  all.Y <- model.response(model.frame(out_lme$terms, data = data))
  data_matrix <- model.matrix(out_lme$terms, data = data)
  D <- data_matrix[, 1:n_class]
  all.Z <- data_matrix[,id_coef]
  Y <- split(all.Y, Clus)
  Z <- lapply(split(as.data.frame(all.Z), Clus), as.matrix)
  V <- lapply(1:cls, function(x) rep(1, sum(Clus == x)))
  D_f <- apply(D, 1, function(x) which(x == 1))
  if(residual){
    resid <- out_lme$residual[,2]
    fit$residual <- split(resid, D_f)
  }
  if(fitted){
    fitted <- out_lme$fitted[,2]
    fit$fitted <- split(fitted, D_f)
  }
  if(randf){
    fit$randf <- ranef(out_lme)$`(Intercept)`
  }
  if(apVar){
    jac <- jacobian(reml_loglik_vec, x = par_lme_est, D = D, Y = Y, Z = Z, V = V, Clus = Clus,
                    cls = cls, n_p = n_p, n_class = n_class)
    hes <- hessian(reml_loglik_vec_sum, x = par_lme_est, D = D, Y = Y, Z = Z, V = V, Clus = Clus,
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
    infer_tab <- cbind(x$est_para, c(x$se_para[1:x$n_coef], rep(NA, length(x$est_para) - x$n_coef)),
                       c(z, rep(NA, length(x$est_para) - x$n_coef)),
                       c(p_val, rep(NA, length(x$est_para) - x$n_coef)))
    colnames(infer_tab) <- c("Est.", "Std.Error", "z-value", "p-value")
    rownames(infer_tab) <- names(x$est_para)
    printCoefmat(infer_tab, has.Pvalue = TRUE, digits = digits, na.print = "--")
  }
  else{
    infer_tab <- cbind(x$est_para)
    colnames(infer_tab) <- c("Est.")
    rownames(infer_tab) <- names(x$est_para)
    printCoefmat(infer_tab, has.Pvalue = FALSE, digits = digits)
  }
  invisible(x)
}

