####==========================================================================####
## This file consists of functions for estimating covariate-specific AUC        ##
## Date: 26/05/2021																															##
####==========================================================================####

auc_norm <- function(mu, sigma2){
  pnorm((mu[2] - mu[1])/sqrt(sum(sigma2)))
}

theta_auc_1 <- function(par, x.val, n_p){
  beta_d <- matrix(par[1:(2*n_p)], ncol = n_p, nrow = 2, byrow = TRUE)
  sigma_e <- par[(2*n_p + 2):length(par)]
  sigma_c <- par[(2*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  return(pnorm((mu_2 - mu_1)/sqrt(sum(sigma_e^2))))
}

theta_auc_2 <- function(par, x.val, n_p){
  beta_d <- matrix(par[1:(2*n_p)], ncol = n_p, nrow = 2, byrow = TRUE)
  sigma_e <- par[(2*n_p + 2):length(par)]
  sigma_c <- par[(2*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  return(pnorm((mu_2 - mu_1)/sqrt(sum(sigma_e^2) + 2*sigma_c^2)))
}


## auc_core function
auc_core <- function(par, x.val, n_p, p.ss, p.sk){
  theta_est_1 <- theta_auc_1(par = par, x.val = x.val, n_p = n_p)
  theta_est_2 <- theta_auc_2(par = par, x.val = x.val, n_p = n_p)
  auc_est <- theta_est_1*p.ss + theta_est_2*p.sk
  return(auc_est)
}

auc_se <- function(par, vcov_par_model, x.val, n_p, p.ss, p.sk){
  jac_auc <- rbind(jacobian(theta_auc_1, x = par, x.val = x.val, n_p = n_p),
                   jacobian(theta_auc_2, x = par, x.val = x.val, n_p = n_p))
  Sig_1 <- jac_auc %*% vcov_par_model %*% t(jac_auc)
  pp <- c(p.ss, p.sk)
  auc_sd <- as.numeric(sqrt(pp %*% Sig_1 %*% pp))
  return(auc_sd)
}

#' @export
AUC <- function(out_lme2, x.val, apVar = FALSE, ci = FALSE, ci.level = ifelse(ci, 0.95, NULL), ...){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  if(out_lme2$n_coef/out_lme2$n_p != 2) stop("There is not a case of two-class setting!")
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
    n_k.div.N <- (unique.n_k/n)^2;
    p.ss <- n_c*n_k.div.N;
    p.sk <- n_c*(n_c - 1)*n_k.div.N
  } else{
    n_k2 <- n_k^2
    p.ss <- sum(n_k2)/n^2
    combo <- combn(x = 1:n_c, m = 2)
    p.sk <- sum(apply(combo, 2, function(idx) prod(n_k[idx])))/(n^2)
    p.sk <- 2*p.sk
  }
  par <- out_lme2$est_para[1:(out_lme2$n_coef + 3)]
  names(par) <- NULL
  if(n_p == 1){ # no covariate
    if(!missing(x.val)) {
      if(!is.null(x.val)) warning("Sepecified value(s) of covariate(s) are not used!", call. = FALSE)
    }
    x.val <- NULL
    fit$auc_est <- auc_core(par = par, x.val = x.val, n_p = n_p, p.ss = p.ss, p.sk = p.sk)
    if(apVar){
      if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing! You may try with bootstrap = TRUE!")
      if(all(is.na(out_lme2$vcov_sand))) stop("Unable to estimate se of AUC. You may try with bootstrap procedure!")
      vcov_par_model <- out_lme2$vcov_sand[1:(out_lme2$n_coef + 3), 1:(out_lme2$n_coef + 3)]
      fit$auc_se <- auc_se(par = par, vcov_par_model = vcov_par_model, x.val = x.val, n_p = n_p, p.ss = p.ss,
                           p.sk = p.sk)
      if(ci){
        ## Normal-approach with truncated boundary
        temp <- fit$auc_est + c(-1, 1)*qnorm((1 + ci.level)/2)* fit$auc_se
        if(temp[1] < 0) temp[1] <- 0
        if(temp[2] > 1) temp[2] <- 1
        fit$auc_ci_norm <- matrix(temp, ncol = 2)
        ## logit-transform
        logit_auc <- qlogis(fit$auc_est)
        logit_auc_sd <- fit$auc_se/(fit$auc_est*(1 - fit$auc_est))
        fit$auc_ci_log <- matrix(plogis(logit_auc + c(-1, 1)*qnorm((1 + ci.level)/2)*logit_auc_sd), ncol = 2)
        ## probit-transform
        probit_auc <- qnorm(fit$auc_est)
        probit_auc_sd <- grad(function(x) qnorm(x), x = fit$auc_est)*fit$auc_se
        fit$auc_ci_prob <- matrix(pnorm(probit_auc + c(-1, 1)*qnorm((1 + ci.level)/2)*probit_auc_sd), ncol = 2)
        fit$ci.level <- ci.level
      }
    }
  }
  if(n_p == 2){ # 1 covariate
    if(missing(x.val)) stop("Please input specific value(s) of covariate.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariate.")
    if(!inherits(x.val, "numeric")) stop("For case of 1 covariate, please input a number or a vector.")
    if(any(is.na(x.val))) stop("NA value(s) not allowed!")
    auc_core_vector <- Vectorize(auc_core, vectorize.args = "x.val")
    fit$auc_est <- auc_core_vector(par = par, x.val = x.val, n_p = n_p, p.ss = p.ss, p.sk = p.sk)
    if(apVar){
      if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing! You may try with bootstrap = TRUE!")
      if(all(is.na(out_lme2$vcov_sand))) stop("Unable to estimate se of AUC. You may try with bootstrap procedure!")
      vcov_par_model <- out_lme2$vcov_sand[1:(out_lme2$n_coef + 3), 1:(out_lme2$n_coef + 3)]
      auc_se_vector <- Vectorize(auc_se, vectorize.args = "x.val")
      fit$auc_se <- auc_se_vector(par = par, vcov_par_model = vcov_par_model, x.val = x.val, n_p = n_p,
                                  p.ss = p.ss, p.sk = p.sk)
      if(ci){
        ## Normal-approach with truncated boundary
        fit$auc_ci_norm <- t(mapply(FUN = function(x, y) {
          res <- x + c(-1, 1)*qnorm((1 + ci.level)/2)*y
          if(res[1] < 0) res[1] <- 0
          if(res[2] > 1) res[2] <- 1
          return(res)
        }, x = fit$auc_est, y = fit$auc_se))
        ## logit-transform
        logit_auc <- qlogis(fit$auc_est)
        logit_auc_sd <- fit$auc_se/(fit$auc_est*(1 - fit$auc_est))
        fit$auc_ci_log <- t(mapply(FUN = function(x, y) plogis(x + c(-1, 1)*qnorm((1 + ci.level)/2)*y),
                                   x = logit_auc, y = logit_auc_sd))
        ## probit-transform
        probit_auc <- qnorm(fit$auc_est)
        probit_auc_sd <- fit$auc_se/(dnorm(probit_auc))
        fit$auc_ci_prob <- t(mapply(FUN = function(x, y) pnorm(x + c(-1, 1)*qnorm((1 + ci.level)/2)*y),
                                    x = probit_auc, y = probit_auc_sd))
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
    fit$auc_est <- apply(x.val, 1, function(x){
      auc_core(par = par, x.val = x, n_p = n_p, p.ss = p.ss, p.sk = p.sk)
    })
    if(apVar){
      if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing! You may try with bootstrap = TRUE!")
      if(all(is.na(out_lme2$vcov_sand))) stop("Unable to estimate se of AUC. You may try with bootstrap procedure!")
      vcov_par_model <- out_lme2$vcov_sand[1:(out_lme2$n_coef + 3), 1:(out_lme2$n_coef + 3)]
      fit$auc_se <- apply(x.val, 1, function(x){
        auc_se(par = par, vcov_par_model = vcov_par_model, x.val = x, n_p = n_p, p.ss = p.ss, p.sk = p.sk)
      })
      if(ci){
        ## Normal-approach with truncated boundary,
        fit$auc_ci_norm <- t(mapply(FUN = function(x, y) {
          res <- x + c(-1, 1)*qnorm((1 + ci.level)/2)*y
          if(res[1] < 0) res[1] <- 0
          if(res[2] > 1) res[2] <- 1
          return(res)
        }, x = fit$auc_est, y = fit$auc_se))
        ## logit-transform
        logit_auc <- qlogis(fit$auc_est)
        logit_auc_sd <- fit$auc_se/(fit$auc_est*(1 - fit$auc_est))
        fit$auc_ci_log <- t(mapply(FUN = function(x, y) plogis(x + c(-1, 1)*qnorm((1 + ci.level)/2)*y),
                                   x = logit_auc, y = logit_auc_sd))
        ## probit-transform
        probit_auc <- qnorm(fit$auc_est)
        probit_auc_sd <- grad(function(x) qnorm(x), x = fit$auc_est)*fit$auc_se
        fit$auc_ci_prob <- t(mapply(FUN = function(x, y) pnorm(x + c(-1, 1)*qnorm((1 + ci.level)/2)*y),
                                    x = probit_auc, y = probit_auc_sd))
        fit$ci.level <- ci.level
      }
    }
  }
  fit$x.val <- x.val
  fit$n_p <- n_p
  class(fit) <- "AUC"
  return(fit)
}

## ---- The function print.AUC ----
#' @title Print summary results of AUC
#'
#' @description \code{print.AUC} prints the results for the output of function \code{\link{AUC}}.
#'
#' @method print AUC
#' @param x an object of class "AUC", a result of a call to \code{\link{AUC}}.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.auc} shows a nice format of the summary table for covariate-specific AUC estimates.
#'
#' @seealso \code{\link{AUC}}
#'
#' @export
print.AUC <- function(x, digits = max(3, getOption("digits") - 2), dig.tst = max(1, min(5, digits - 1)), ...){
  if(isFALSE(inherits(x, "AUC"))) stop("The object is not AUC!")
  cat("\n")
  cat("CALL: ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  if(x$n_p == 1){
    labels <- "Intercept"
  }
  if(x$n_p == 2) {
    labels <- format(x$x.val, digits = digits)
  }
  if(x$n_p > 2) {
    labels <- apply(x$x.val, 1, function(y) paste0("(", paste(y, collapse = ", "), ")"))
  }
  if(!is.null(x$auc_se)){
    z <- (x$auc_est - rep(1/2, length(x$auc_est)))/x$auc_se
    p_val <- 2*pnorm(abs(z), lower.tail = FALSE)
    infer_tab <- data.frame(labels, x$auc_est, x$auc_se, z, p_val) # as.factor(labels),
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
    cat("Covariate-specific AUC: \n")
    print(infer_tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, ...)
    cat("---\nSignif. codes:  ", sleg, sep = "",
        fill = getOption("width") + 4 + max(nchar(sleg, "bytes") - nchar(sleg)))
    cat("z-value and p-value are for testing the null hypothesis H0: AUC = 1/2 \n")
    # rownames(infer_tab) <- labels
    # printCoefmat(infer_tab, has.Pvalue = TRUE, digits = digits, na.print = "--") # cs.ind = 2:3, tst.ind = 4,
    if(!is.null(x$auc_ci_norm)){
      ci.tab <- cbind(x$auc_ci_norm, x$auc_ci_log, x$auc_ci_prob)
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
      cat(paste0("The ", x$ci.level*100, "% confidence Intervals for covariate-specific AUC:\n"))
      print(res.ci.tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, print.gap = 3)
    }
  }
  else{
    infer_tab <- data.frame(labels, x$auc_est)
    infer_tab[,2] <- signif(infer_tab[,2], digits = digits)
    colnames(infer_tab) <- c("Covariates Values", "Est.")
    # rownames(infer_tab) <- labels
    cat("Covariate-specific AUC: \n")
    print(infer_tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, ...)
  }
  cat("\n")
  invisible(x)
}
