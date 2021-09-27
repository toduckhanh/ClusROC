####==========================================================================####
## This file consists of functions for estimating covariate-specific VUS        ##
## Date: 26/03/2021																															##
####==========================================================================####

#' @import utils
#' @import numDeriv
#' @import stats

theta_1 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f01 <- function(s){
    dnorm((s - mu_2)/sigma_e[2])*pnorm((s - mu_1)/sigma_e[1])*
      pnorm((s - mu_3)/sigma_e[3], lower.tail = FALSE)/sigma_e[2]
  }
  return(integrate(f01, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

theta_2 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f02 <- function(s){
    dnorm((s - mu_2)/sigma_e[2])*pnorm((s - mu_1)/sigma_e[1])*pnorm((s - mu_3)/sqrt(2*sigma_c^2 + sigma_e[3]^2),
                                                                    lower.tail = FALSE)/sigma_e[2]
  }
  return(integrate(f02, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

theta_3 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f03 <- function(s){
    dnorm((s - mu_2)/sqrt(2*sigma_c^2 + sigma_e[2]^2))*pnorm((s - mu_1)/sigma_e[1])*
      pnorm((s - mu_3)/sigma_e[3], lower.tail = FALSE)/sqrt(2*sigma_c^2 + sigma_e[2]^2)
  }
  return(integrate(f03, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

theta_4 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f04 <- function(s){
    dnorm((s - mu_2)/sigma_e[2])*pnorm((s - mu_1)/sqrt(2*sigma_c^2 + sigma_e[1]^2))*
      pnorm((s - mu_3)/sigma_e[3], lower.tail = FALSE)/sigma_e[2]
  }
  return(integrate(f04, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

theta_5 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f05 <- function(s){
    dnorm((s - mu_2)/sqrt(sigma_c^2 + sigma_e[2]^2))*pnorm((s - mu_1)/sqrt(sigma_c^2 + sigma_e[1]^2))*
      pnorm((s - mu_3)/sqrt(sigma_c^2 + sigma_e[3]^2), lower.tail = FALSE)/sqrt(sigma_c^2 + sigma_e[2]^2)
  }
  return(integrate(f05, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

## vus_core function
vus_core <- function(par, x.val, n_p, n_c, n_k, n, p.sss, p.ssk, p.sks, p.skk, p.ijk,
                     subdivisions = 1000, ...){
  theta_est_1 <- theta_1(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  theta_est_2 <- theta_2(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  theta_est_3 <- theta_3(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  theta_est_4 <- theta_4(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  theta_est_5 <- theta_5(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  vus_est <- theta_est_1*p.sss + theta_est_2*p.ssk + theta_est_3*p.sks + theta_est_4*p.skk + theta_est_5*p.ijk
  return(vus_est)
}

vus_se <- function(par, vcov_par_model, x.val, n_p, n_c, n_k, n, p.sss, p.ssk, p.sks, p.skk, p.ijk){
  jac_vus <- rbind(jacobian(theta_1, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n),
                   jacobian(theta_2, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n),
                   jacobian(theta_3, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n),
                   jacobian(theta_4, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n),
                   jacobian(theta_5, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n))
  Sig_1 <- jac_vus %*% vcov_par_model %*% t(jac_vus)
  pp <- c(p.sss, p.ssk, p.sks, p.skk, p.ijk)
  vus_sd <- as.numeric(sqrt(pp %*% Sig_1 %*% pp))
  return(vus_sd)
}

#' @export
VUS <- function(out_lme2, x.val, apVar = FALSE, ci = FALSE, ci.level = ifelse(ci, 0.95, NULL),
                subdivisions = 1000, ...){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  if(out_lme2$n_coef/out_lme2$n_p != 3) stop("There is not a case of three-class setting!")
  if(isFALSE(apVar) & ci) stop("Confidence intervals cannot computed without option of apVar!")
  call <- match.call()
  fit <- list()
  fit$call <- call
  ## check if all clusters/families have the same cluster size, if yes, the weight calculation_k is simplified
  n_p <- out_lme2$n_p
  n_c <- out_lme2$cls
  n_k <- out_lme2$n_c
  n <- out_lme2$n
  unique.n_k <- unique(n_k)
  equal.n_k.flag <- ifelse(length(unique.n_k) == 1, TRUE, FALSE)
  if(equal.n_k.flag){
    n_k.div.N <- (unique.n_k/n)^3;
    p.sss <- n_c*n_k.div.N;
    p.ssk <- n_c*(n_c - 1)*n_k.div.N
    p.ijk <- n_c*(n_c - 1)*(n_c - 2)*n_k.div.N
  } else{
    n_k3 <- n_k^3
    p.sss <- sum(n_k3)/n^3
    n_k.sqr <- n_k^2
    outer.n_kSQR.n_k <- outer(n_k.sqr, n_k)
    diag(outer.n_kSQR.n_k) <- NA
    p.ssk <- sum(outer.n_kSQR.n_k, na.rm = TRUE)/(n^3)
    combo <- combn(x = 1:n_c, m = 3)
    p.ijk <- sum(apply(combo, 2, function(idx) prod(n_k[idx])))/(n^3)
    p.ijk <- 6*p.ijk
  }
  p.skk <- p.sks <- p.ssk
  par <- out_lme2$est_para[1:(out_lme2$n_coef + 4)]
  ##
  if(n_p == 1){ # no covariate
    if(!missing(x.val)) {
      if(!is.null(x.val)) warning("Sepecified value(s) of covariate(s) are not used!", call. = FALSE)
    }
    x.val <- NULL
    fit$vus_est <- vus_core(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n, p.sss = p.sss,
                            p.ssk = p.ssk, p.sks = p.sks, p.skk = p.skk, p.ijk = p.ijk,
                            subdivisions = subdivisions, ...)
    if(apVar){
      if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing! You may try with bootstrap = TRUE!")
      if(all(is.na(out_lme2$vcov_sand))) stop("Unable to estimate se of VUS. You may try with bootstrap procedure!")
      vcov_par_model <- out_lme2$vcov_sand[1:(out_lme2$n_coef + 4), 1:(out_lme2$n_coef + 4)]
      fit$vus_se <- vus_se(par = par, vcov_par_model = vcov_par_model, x.val = x.val, n_p = n_p, n_c = n_c,
                           n_k = n_k, n = n, p.sss = p.sss, p.ssk = p.ssk, p.sks = p.sks, p.skk = p.skk,
                           p.ijk = p.ijk)
      if(ci){
        ## Normal-approach with truncated boundary
        temp <- fit$vus_est + c(-1, 1)*qnorm((1 + ci.level)/2)* fit$vus_se
        if(temp[1] < 0) temp[1] <- 0
        if(temp[2] > 1) temp[2] <- 1
        fit$vus_ci_norm <- matrix(temp, ncol = 2)
        ## logit-transform
        logit_vus <- qlogis(fit$vus_est)
        logit_vus_sd <- fit$vus_se/(fit$vus_est*(1 - fit$vus_est))
        fit$vus_ci_log <- matrix(plogis(logit_vus + c(-1, 1)*qnorm((1 + ci.level)/2)*logit_vus_sd), ncol = 2)
        ## probit-transform
        probit_vus <- qnorm(fit$vus_est)
        probit_vus_sd <- grad(function(x) qnorm(x), x = fit$vus_est)*fit$vus_se
        fit$vus_ci_prob <- matrix(pnorm(probit_vus + c(-1, 1)*qnorm((1 + ci.level)/2)*probit_vus_sd), ncol = 2)
        fit$ci.level <- ci.level
      }
    }
  }
  if(n_p == 2){ # 1 covariate
    if(missing(x.val)) stop("Please input specific value(s) of covariate.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariate.")
    if(!inherits(x.val, "numeric")) stop("For case of 1 covariate, please input a number or a vector.")
    if(any(is.na(x.val))) stop("NA value(s) not allowed!")
    vus_core_vector <- Vectorize(vus_core, vectorize.args = "x.val")
    fit$vus_est <- vus_core_vector(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                                   p.sss = p.sss, p.ssk = p.ssk, p.sks = p.sks, p.skk = p.skk, p.ijk = p.ijk,
                                   subdivisions = subdivisions, ...)
    if(apVar){
      if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing! You may try with bootstrap = TRUE!")
      if(all(is.na(out_lme2$vcov_sand))) stop("Unable to estimate se of VUS. You may try with bootstrap procedure!")
      vcov_par_model <- out_lme2$vcov_sand[1:(out_lme2$n_coef + 4), 1:(out_lme2$n_coef + 4)]
      vus_se_vector <- Vectorize(vus_se, vectorize.args = "x.val")
      fit$vus_se <- vus_se_vector(par = par, vcov_par_model = vcov_par_model, x.val = x.val, n_p = n_p,
                                  n_c = n_c, n_k = n_k, n = n, p.sss = p.sss, p.ssk = p.ssk, p.sks = p.sks,
                                  p.skk = p.skk, p.ijk = p.ijk)
      if(ci){
        ## Normal-approach with truncated boundary
        fit$vus_ci_norm <- t(mapply(FUN = function(x, y) {
          res <- x + c(-1, 1)*qnorm((1 + ci.level)/2)*y
          if(res[1] < 0) res[1] <- 0
          if(res[2] > 1) res[2] <- 1
          return(res)
        }, x = fit$vus_est, y = fit$vus_se))
        ## logit-transform
        logit_vus <- qlogis(fit$vus_est)
        logit_vus_sd <- fit$vus_se/(fit$vus_est*(1 - fit$vus_est))
        fit$vus_ci_log <- t(mapply(FUN = function(x, y) plogis(x + c(-1, 1)*qnorm((1 + ci.level)/2)*y),
                                   x = logit_vus, y = logit_vus_sd))
        ## probit-transform
        probit_vus <- qnorm(fit$vus_est)
        probit_vus_sd <- grad(function(x) qnorm(x), x = fit$vus_est)*fit$vus_se
        fit$vus_ci_prob <- t(mapply(FUN = function(x, y) pnorm(x + c(-1, 1)*qnorm((1 + ci.level)/2)*y),
                                   x = probit_vus, y = probit_vus_sd))
        fit$ci.level <- ci.level
      }
    }
  }
  if(n_p > 2){ # multiple covariates
    if(missing(x.val)) stop("Please input specific value(s) of covariates.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariates.")
    if(inherits(x.val, "numeric")){
      if(length(x.val) != n_p - 1) stop(paste("For case of", n_p - 1, "covariates, please input a vector of", n_p - 1, "values of covariates."))
    }
    if(inherits(x.val, "matrix")) {
      if(ncol(x.val) != n_p - 1) stop(paste("For case of m points of", n_p - 1, "covariates, please input a matrix with", n_p - 1, "columns and m rows containing values of covariates."))
    }
    if(any(is.na(x.val))) stop("NA value(s) not allowed!")
    x.val <- matrix(x.val, ncol = n_p - 1, byrow = FALSE)
    fit$vus_est <- apply(x.val, 1, function(x){
      vus_core(par = par, x.val = x, n_p = n_p, n_c = n_c, n_k = n_k, n = n, p.sss = p.sss, p.ssk = p.ssk,
               p.sks = p.sks, p.skk = p.skk, p.ijk = p.ijk, subdivisions = subdivisions, ...)
    })
    if(apVar){
      if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing! You may try with bootstrap = TRUE!")
      if(all(is.na(out_lme2$vcov_sand))) stop("Unable to estimate se of VUS. You may try with bootstrap procedure!")
      vcov_par_model <- out_lme2$vcov_sand[1:(out_lme2$n_coef + 4), 1:(out_lme2$n_coef + 4)]
      vus_se_vector <- Vectorize(vus_se, vectorize.args = "x.val")
      fit$vus_se <- apply(x.val, 1, function(x){
        vus_se(par = par, vcov_par_model = vcov_par_model, x.val = x, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
               p.sss = p.sss, p.ssk = p.ssk, p.sks = p.sks, p.skk = p.skk, p.ijk = p.ijk)
      })
      if(ci){
        ## Normal-approach with truncated boundary,
        fit$vus_ci_norm <- t(mapply(FUN = function(x, y) {
          res <- x + c(-1, 1)*qnorm((1 + ci.level)/2)*y
          if(res[1] < 0) res[1] <- 0
          if(res[2] > 1) res[2] <- 1
          return(res)
          }, x = fit$vus_est, y = fit$vus_se))
        ## logit-transform
        logit_vus <- qlogis(fit$vus_est)
        logit_vus_sd <- fit$vus_se/(fit$vus_est*(1 - fit$vus_est))
        fit$vus_ci_log <- t(mapply(FUN = function(x, y) plogis(x + c(-1, 1)*qnorm((1 + ci.level)/2)*y),
                                   x = logit_vus, y = logit_vus_sd))
        ## probit-transform
        probit_vus <- qnorm(fit$vus_est)
        probit_vus_sd <- grad(function(x) qnorm(x), x = fit$vus_est)*fit$vus_se
        fit$vus_ci_prob <- t(mapply(FUN = function(x, y) pnorm(x + c(-1, 1)*qnorm((1 + ci.level)/2)*y),
                                    x = probit_vus, y = probit_vus_sd))
        fit$ci.level <- ci.level
      }
    }
  }
  fit$x.val <- x.val
  fit$n_p <- n_p
  class(fit) <- "VUS"
  return(fit)
}

## ---- The function print.VUS ----
#' @title Print summary results of VUS
#'
#' @description \code{print.VUS} prints the results for the output of function \code{\link{VUS}}.
#'
#' @method print VUS
#' @param x an object of class "VUS", a result of a call to \code{\link{VUS}}.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.vus} shows a nice format of the summary table for covariate-specific VUS estimates.
#'
#' @seealso \code{\link{VUS}}
#'
#' @export
print.VUS <- function(x, digits = max(3, getOption("digits") - 2), dig.tst = max(1, min(5, digits - 1)), ...){
  if(isFALSE(inherits(x, "VUS"))) stop("The object is not VUS!")
  cat("\n")
  cat("CALL: ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  if(x$n_p == 1){
    labels <- "Intercept"
  }
  if(x$n_p == 2) {
    labels <- as.character(x$x.val)
  }
  if(x$n_p > 2) {
    labels <- apply(x$x.val, 1, function(y) paste0("(", paste(y, collapse = ", "), ")"))
  }
  if(!is.null(x$vus_se)){
    z <- (x$vus_est - rep(1/6, length(x$vus_est)))/x$vus_se
    p_val <- 2*pnorm(abs(z), lower.tail = FALSE)
    infer_tab <- data.frame(labels, x$vus_est, x$vus_se, z, p_val) # as.factor(labels),
    infer_tab[,2:4] <- signif(infer_tab[,2:4], digits = digits)
    pv <- as.vector(infer_tab[,5])
    infer_tab[,5] <- format.pval(pv, digits = dig.tst, eps = 0.001)
    Signif <- symnum(pv, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " "))
    sleg <- attr(Signif, "legend")
    sleg <- strwrap(sleg, width = getOption("width") - 2, prefix = "  ")
    infer_tab <- cbind(infer_tab, format(Signif))
    colnames(infer_tab) <- c("Covariates Values", "Est.", "Std.Error", "z-value", "p-value", "")
    cat("Covariate-specific VUS: \n")
    print(infer_tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, ...)
    cat("---\nSignif. codes:  ", sleg, sep = "",
        fill = getOption("width") + 4 + max(nchar(sleg, "bytes") - nchar(sleg)))
    cat("z-value and p-value are for testing the null hypothesis H0: VUS = 1/6 \n")
    # rownames(infer_tab) <- labels
    # printCoefmat(infer_tab, has.Pvalue = TRUE, digits = digits, na.print = "--") # cs.ind = 2:3, tst.ind = 4,
    if(!is.null(x$vus_ci_norm)){
      ci.tab <- cbind(x$vus_ci_norm, x$vus_ci_log, x$vus_ci_prob)
      ci.tab <- format(round(ci.tab, digits = digits))
      res.ci.tab <- data.frame(labels,
                               apply(matrix(ci.tab[,1:2], ncol = 2, byrow = FALSE), 1,
                                     function(y) paste0("(", paste(y, collapse = ", "), ")")),
                               apply(matrix(ci.tab[,3:4], ncol = 2, byrow = FALSE), 1,
                                     function(y) paste0("(", paste(y, collapse = ", "), ")")),
                               apply(matrix(ci.tab[,5:6], ncol = 2, byrow = FALSE), 1,
                                     function(y) paste0("(", paste(y, collapse = ", "), ")")))
      colnames(res.ci.tab) <- c("Covariates Values", "Normal approximation", "Logit transformation",
                                "Probit transformation")
      cat("---\n")
      cat(paste0("The ", x$ci.level*100, "% confidence Intervals for covariate-specific VUS:\n"))
      print(res.ci.tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, print.gap = 3)
    }
  }
  else{
    infer_tab <- data.frame(labels, x$vus_est)
    infer_tab[,2] <- signif(infer_tab[,2], digits = digits)
    colnames(infer_tab) <- c("Covariates Values", "Est.")
    # rownames(infer_tab) <- labels
    cat("Covariate-specific VUS: \n")
    print(infer_tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, ...)
  }
  cat("\n")
  invisible(x)
}



