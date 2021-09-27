####==========================================================================####
## This file consists of functions for estimating optimal threshold             ##
## Based on NC method (Northwest Corner)                                        ##
## Date: 28/05/2021																															##
####==========================================================================####

NC_normal <- function(par, par_model, x.val, n_p){
  beta_d <- matrix(par_model[1:(2*n_p)], ncol = n_p, nrow = 2, byrow = TRUE)
  sigma_d <- sqrt(par_model[(2*n_p + 2):length(par_model)]^2 + par_model[(2*n_p + 1)]^2)
  Z <- c(1, x.val)
  res_tcfs <- numeric(2)
  res_tcfs[1] <- pnorm(par, mean = Z %*% beta_d[1,], sd = sigma_d[1])
  res_tcfs[2] <- 1 - pnorm(par, mean = Z %*% beta_d[2,], sd = sigma_d[2])
  return(sqrt(sum((1 - res_tcfs)^2)))
}

NC_normal_var <- function(opt_ctps, par_model, vcov_par_model, x.val, n_p){
  term1 <- grad(NC_normal, x = par_model, par = opt_ctps, x.val = x.val, n_p = n_p)
  term2 <- grad(NC_normal, x = opt_ctps, par_model = par_model, x.val = x.val, n_p = n_p)
  grad.tcfs_par <- jacobian(TCF2_normal, x = par_model, x.val = x.val, threshold = opt_ctps, n_p = n_p)
  term3 <- hessian(func = function(par, x.val, n_p){
    NC_normal(par = par[1], par_model = par[-c(1)], x.val = x.val, n_p = n_p)
  }, x = c(opt_ctps, par_model), x.val = x.val, n_p = n_p)
  term3.1 <- term3[1, 1]
  term4 <- jacobian(TCF2_normal, x = opt_ctps, x.val = x.val, par = par_model, n_p = n_p)
  derv.cpt.para <- -solve(term3.1) %*% term3[1, -c(1)]
  derv.CtP.para <- term1 + term2 * derv.cpt.para
  derv.TCF.para <- grad.tcfs_par + term4 %*% derv.cpt.para
  vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  var_D <- as.numeric(derv.CtP.para %*% vcov_par_model %*% t(derv.CtP.para))
  var_cpts <- as.numeric(derv.cpt.para %*% vcov_par_model %*% t(derv.cpt.para))
  return(list(var_D = var_D, var_cpts = var_cpts, vcov.tcfs = vcov.tcfs))
}

NC_est_fun <- function(par_model, x.val, n_p, boxcox, par_boxcox, start = NULL, method.optim,
                       maxit, lower, upper){
  beta_d <- matrix(par_model[1:(2*n_p)], ncol = n_p, nrow = 2, byrow = TRUE)
  Z <- c(1, x.val)
  mu <- c(Z %*% beta_d[1,], Z %*% beta_d[2,])
  sigma_d <- sqrt(par_model[(2*n_p + 2):length(par_model)]^2 + par_model[(2*n_p + 1)]^2)
  sigma2_d <- sigma_d^2
  if(is.null(start)) start <- mu[1] + qnorm(0.3)*sigma_d[1]
  #
  NC_optim <- optim(par = start, fn = NC_normal, par_model = par_model, x.val = x.val, n_p = n_p,
                    method = method.optim, control = list(maxit = maxit), lower = lower, upper = upper)
  cpt_NC_est <- NC_optim$par
  if(boxcox) cpt_NC_est <- boxcox_trans_back(cpt_NC_est, par_boxcox)
  out <- cpt_NC_est
  names(out) <- c("threshold")
  return(out)
}

