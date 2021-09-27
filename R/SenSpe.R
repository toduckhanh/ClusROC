####==========================================================================####
## This file consists of functions for estimating Sensitivity and Specificity   ##
## at given fixed threshold or optimal threshold                                ##
## Date: 03/06/2021																															##
####==========================================================================####

### Sp and Se
TCF2_normal <- function(par, x.val, threshold, n_p, boxcox = FALSE){
  beta_d <- matrix(par[1:(2*n_p)], ncol = n_p, nrow = 2, byrow = TRUE)
  sigma_d <- sqrt(par[(2*n_p + 2):(2*n_p + 3)]^2 + par[(2*n_p + 1)]^2)
  if(boxcox) threshold <- boxcox_trans(thresholds, par[length(par)])
  Z <- c(1, x.val)
  res_tcfs <- numeric(2)
  res_tcfs[1] <- pnorm(threshold, mean = Z %*% beta_d[1,], sd = sigma_d[1])
  res_tcfs[2] <- 1 - pnorm(threshold, mean = Z %*% beta_d[2,], sd = sigma_d[2])
  return(res_tcfs)
}


#' @export
SenSpe <- function(out_lme2, out_optThres2, apVar = FALSE, ci = FALSE, ci.level = ifelse(ci, 0.95, NULL), ...){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  if(out_lme2$n_coef/out_lme2$n_p != 2) stop("There is not a case of two-class setting!")
  if(isFALSE(inherits(out_optThres2, "optThres2"))) stop("out_lme2 was not from optThres2()!")
  if(isFALSE(apVar) & ci) stop("Confidence intervals cannot computed without option of apVar!")
  call <- match.call()
  fit <- list()
  fit$call <- call
  ##
  par <- out_lme2$est_para
  dt_thres2 <- out_optThres2$thres2
  sp_se_mt <- mapply(function(x, y){
    TCF2_normal(par = par, x.val = x, threshold = y, n_p = out_lme2$n_p, boxcox = out_lme2$boxcox)
  }, x = dt_thres2$x, y = dt_thres2$threshold)
  fit$sp_se <- data.frame(threshold = dt_thres2$threshold, specificity = sp_se_mt[1,],
                          sensitivity = sp_se_mt[2,], Method = dt_thres2$Method, x = dt_thres2$x)
  if(apVar){
    if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing!")
    dt_thres2_list <- split(dt_thres2, dt_thres2$Method)
    fit$vcov.sp.se <- list()
    if("YI" %in% method){
      YI_normal_var(opt_ctps = -1.370565, par_model = out_2class_md1$est_para, vcov_par_model = out_2class_md1$vcov_sand, x.val = -2.5, n_p = out_2class_md1$n_p)
    }
  }
  class(fit) <- "SenSpe"
  return(fit)
}


