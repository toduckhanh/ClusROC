####==========================================================================####
## This file consists of functions for estimating optimal pair of thresholds    ##
## Based on MV method                                                           ##
## Date: 19/05/2021																															##
####==========================================================================####

MV_normal <- function(par, par_model, z, n_p){
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
  return(prod(res_tcfs))
}

MV_normal_var <- function(opt_ctps, par_model, vcov_par_model, z, n_p){
  term1 <- grad(MV_normal, x = par_model, par = opt_ctps, z = z, n_p = n_p)
  term2 <- grad(MV_normal, x = opt_ctps, par_model = par_model, z = z, n_p = n_p)
  grad.tcfs_par <- jacobian(TCF_normal, x = par_model, z = z, thresholds = opt_ctps, n_p = n_p)
  term3 <- hessian(func = function(par, z, n_p){
    MV_normal(par = par[1:2], par_model = par[-c(1,2)], z = z, n_p = n_p)
  }, x = c(opt_ctps, par_model), z = z, n_p = n_p)
  term3.1 <- term3[1:2, 1:2]
  term4 <- jacobian(TCF_normal, x = opt_ctps, z = z, par = par_model, n_p = n_p)
  derv.cpt.para <- -solve(term3.1) %*% term3[1:2, -c(1,2)]
  derv.MV.para <- term1 + term2 %*% derv.cpt.para
  derv.TCF.para <- grad.tcfs_par + term4 %*% derv.cpt.para
  vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  var_V3 <- as.numeric(derv.MV.para %*% vcov_par_model %*% t(derv.MV.para))
  vcov_cpts <- derv.cpt.para %*% vcov_par_model %*% t(derv.cpt.para)
  return(list(var_V3 = var_V3, vcov_cpts = vcov_cpts, vcov.tcfs = vcov.tcfs))
}


MV_est_fun <- function(par_model, z, n_p, boxcox, par_boxcox, start = NULL, method.optim,
                       maxit, lower, upper){
  beta_d <- par_model[1:(3*n_p)]
  mu <- z %*% beta_d
  sigma_d <- sqrt(par_model[(3*n_p + 2):length(par_model)]^2 + par_model[(3*n_p + 1)]^2)
  sigma2_d <- sigma_d^2
  if(is.null(start)){
    t1_0 <- mu[1] + qnorm(0.3)*sigma_d[1]
    t2_0 <- mu[3] + qnorm(1 - 0.4)*sigma_d[3]
    start <- c(t1_0, t2_0)
  }
  #
  MV_optim <- optim(par = start, fn = MV_normal, par_model = par_model, z = z, n_p = n_p,
                    method = method.optim, control = list(maxit = maxit, fnscale = -1),
                    lower = lower, upper = upper)
  if(MV_optim$par[1] > MV_optim$par[2]){
    cpt_MV_est <- rev(MV_optim$par)
  } else{
    cpt_MV_est <- MV_optim$par
  }
  if(boxcox) cpt_MV_est <- boxcox_trans_back(cpt_MV_est, par_boxcox)
  out <- cpt_MV_est
  names(out) <- c("threshold_1", "threshold_2")
  return(out)
}


