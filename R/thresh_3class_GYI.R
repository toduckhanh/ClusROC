####==========================================================================####
## This file consists of functions for estimating optimal pair of thresholds    ##
## Based on GYI method                                                          ##
## Date: 19/05/2021																															##
####==========================================================================####

GYI_normal <- function(par, par_model, z, n_p){
  if(par[1] > par[2]){
    x <- par[1]
    par[1] <- par[2]
    par[2] <- x
  }
  beta_d <- par_model[1:(3*n_p)]
  sigma_d <- sqrt(par_model[(3*n_p + 2):length(par_model)]^2 + par_model[(3*n_p + 1)]^2)
  mu_est <- z %*% beta_d
  res_tcfs <- numeric(3)
  res_tcfs[1] <- pnorm(par[1], mean = mu_est[1], sd = sigma_d[1])
  res_tcfs[2] <- pnorm(par[2], mean = mu_est[2], sd = sigma_d[2]) -
    pnorm(par[1], mean = mu_est[2], sd = sigma_d[2])
  res_tcfs[3] <- 1 - pnorm(par[2], mean = mu_est[3], sd = sigma_d[3])
  return(sum(res_tcfs))
}

GYI_normal_var <- function(opt_ctps, par_model, vcov_par_model, z, n_p){
  term1 <- grad(GYI_normal, x = par_model, par = opt_ctps, z = z, n_p = n_p)
  term2 <- grad(GYI_normal, x = opt_ctps, par_model = par_model, z = z, n_p = n_p)
  grad.tcfs_par <- jacobian(TCF_normal, x = par_model, z = z, thresholds = opt_ctps, n_p = n_p)
  term3 <- hessian(func = function(par, z, n_p){
    GYI_normal(par = par[1:2], par_model = par[-c(1,2)], z = z, n_p = n_p)
  }, x = c(opt_ctps, par_model), z = z, n_p = n_p)
  term3.1 <- term3[1:2, 1:2]
  term4 <- jacobian(TCF_normal, x = opt_ctps, z = z, par = par_model, n_p = n_p)
  derv.cpt.para <- -solve(term3.1) %*% term3[1:2, -c(1,2)]
  derv.GYI.para <- term1 + term2 %*% derv.cpt.para
  derv.TCF.para <- grad.tcfs_par + term4 %*% derv.cpt.para
  vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  var_J3 <- as.numeric(derv.GYI.para %*% vcov_par_model %*% t(derv.GYI.para))
  vcov_cpts <- derv.cpt.para %*% vcov_par_model %*% t(derv.cpt.para)
  return(list(var_J3 = var_J3, vcov_cpts = vcov_cpts, vcov.tcfs = vcov.tcfs))
}

GYI_est_fun <- function(par_model, z, n_p, boxcox, par_boxcox, crit.var = 0.01){
  beta_d <- par_model[1:(3*n_p)]
  mu <- z %*% beta_d
  sigma_d <- sqrt(par_model[(3*n_p + 2):length(par_model)]^2 + par_model[(3*n_p + 1)]^2)
  sigma2_d <- sigma_d^2
  #
  var1.check <- abs(sigma2_d[1]/sigma2_d[2] - 1) < crit.var
  var2.check <- abs(sigma2_d[2]/sigma2_d[3] - 1) < crit.var
  if(var1.check){
    cpt_1 <- (mu[1] + mu[2])/2
  } else{
    cpt_1 <- ((mu[2]*sigma2_d[1] - mu[1]*sigma2_d[2]) -
                sigma_d[1]*sigma_d[2]*sqrt((mu[1] - mu[2])^2 +
                                             (sigma2_d[1] - sigma2_d[2])*log(sigma2_d[1]/sigma2_d[2])))/
      (sigma2_d[1] - sigma2_d[2])
  }
  if(var2.check){
    cpt_2 <- (mu[2] + mu[3])/2
  } else{
    cpt_2 <- ((mu[3]*sigma2_d[2] - mu[2]*sigma2_d[3]) -
                sigma_d[2]*sigma_d[3]*sqrt((mu[2] - mu[3])^2 +
                                             (sigma2_d[2] - sigma2_d[3])*log(sigma2_d[2]/sigma2_d[3])))/
      (sigma2_d[2] - sigma2_d[3])
  }
  out <- c(cpt_1, cpt_2)
  if(boxcox){
    cpt_1_org <- boxcox_trans_back(cpt_1, par_boxcox)
    cpt_2_org <- boxcox_trans_back(cpt_2, par_boxcox)
    out <- c(cpt_1_org, cpt_2_org)
  }
  names(out) <- c("threshold_1", "threshold_2")
  return(out)
}


