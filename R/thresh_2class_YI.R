####==========================================================================####
## This file consists of functions for estimating optimal threshold             ##
## Based on YI method                                                           ##
## Date: 28/05/2021																															##
####==========================================================================####

YI_normal <- function(par, par_model, x.val, n_p){
  beta_d <- matrix(par_model[1:(2*n_p)], ncol = n_p, nrow = 2, byrow = TRUE)
  sigma_d <- sqrt(par_model[(2*n_p + 2):length(par_model)]^2 + par_model[(2*n_p + 1)]^2)
  Z <- c(1, x.val)
  res_tcfs <- numeric(2)
  res_tcfs[1] <- pnorm(par, mean = Z %*% beta_d[1,], sd = sigma_d[1])
  res_tcfs[2] <- 1 - pnorm(par, mean = Z %*% beta_d[2,], sd = sigma_d[2])
  return(sum(res_tcfs) - 1)
}

YI_normal_var <- function(opt_ctps, par_model, vcov_par_model, x.val, n_p){
  term1 <- grad(YI_normal, x = par_model, par = opt_ctps, x.val = x.val, n_p = n_p)
  term2 <- grad(YI_normal, x = opt_ctps, par_model = par_model, x.val = x.val, n_p = n_p)
  grad.tcfs_par <- jacobian(TCF2_normal, x = par_model, x.val = x.val, threshold = opt_ctps, n_p = n_p)
  term3 <- hessian(func = function(par, x.val, n_p){
    YI_normal(par = par[1], par_model = par[-c(1)], x.val = x.val, n_p = n_p)
  }, x = c(opt_ctps, par_model), x.val = x.val, n_p = n_p)
  term3.1 <- term3[1, 1]
  term4 <- jacobian(TCF2_normal, x = opt_ctps, x.val = x.val, par = par_model, n_p = n_p)
  derv.cpt.para <- -solve(term3.1) %*% term3[1, -c(1)]
  derv.GYI.para <- term1 + term2 * derv.cpt.para
  derv.TCF.para <- grad.tcfs_par + term4 %*% derv.cpt.para
  vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  var_J <- as.numeric(derv.GYI.para %*% vcov_par_model %*% t(derv.GYI.para))
  var_cpts <- as.numeric(derv.cpt.para %*% vcov_par_model %*% t(derv.cpt.para))
  return(list(var_J = var_J, var_cpts = var_cpts, vcov.tcfs = vcov.tcfs))
}

YI_est_fun <- function(par_model, x.val, n_p, boxcox, par_boxcox, crit.var = 0.01){
  beta_d <- matrix(par_model[1:(2*n_p)], ncol = n_p, nrow = 2, byrow = TRUE)
  Z <- c(1, x.val)
  mu <- c(Z %*% beta_d[1,], Z %*% beta_d[2,])
  sigma_d <- sqrt(par_model[(2*n_p + 2):length(par_model)]^2 + par_model[(2*n_p + 1)]^2)
  sigma2_d <- sigma_d^2
  #
  var.check <- abs(sigma2_d[1]/sigma2_d[2] - 1) < crit.var
  if(var.check){
    cpt <- (mu[1] + mu[2])/2
  } else{
    cpt <- ((mu[2]*sigma2_d[1] - mu[1]*sigma2_d[2]) -
              sigma_d[1]*sigma_d[2]*sqrt((mu[1] - mu[2])^2 +
                                           (sigma2_d[1] - sigma2_d[2])*log(sigma2_d[1]/sigma2_d[2])))/
      (sigma2_d[1] - sigma2_d[2])
  }
  out <- cpt
  if(boxcox){
    cpt_org <- boxcox_trans_back(cpt, par_boxcox)
    out <- cpt
  }
  names(out) <- c("threshold")
  return(out)
}

