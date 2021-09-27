####==========================================================================####
## This file consists of functions for estimating optimal pair of thresholds    ##
## Based on CtP method                                                          ##
## Date: 19/05/2021																															##
####==========================================================================####

CtP_normal <- function(par, par_model, x.val, n_p){
  if(par[1] > par[2]){
    x <- par[1]
    par[1] <- par[2]
    par[2] <- x
  }
  beta_d <- matrix(par_model[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_d <- sqrt(par_model[(3*n_p + 2):length(par_model)]^2 + par_model[(3*n_p + 1)]^2)
  Z <- c(1, x.val)
  res_tcfs <- numeric(3)
  res_tcfs[1] <- pnorm(par[1], mean = Z %*% beta_d[1,], sd = sigma_d[1])
  res_tcfs[2] <- pnorm(par[2], mean = Z %*% beta_d[2,], sd = sigma_d[2]) -
    pnorm(par[1], mean = Z %*% beta_d[2,], sd = sigma_d[2])
  res_tcfs[3] <- 1 - pnorm(par[2], mean = Z %*% beta_d[3,], sd = sigma_d[3])
  return(sqrt(sum((1 - res_tcfs)^2)))
}

CtP_normal_var <- function(opt_ctps, par_model, vcov_par_model, x.val, n_p){
  term1 <- grad(CtP_normal, x = par_model, par = opt_ctps, x.val = x.val, n_p = n_p)
  term2 <- grad(CtP_normal, x = opt_ctps, par_model = par_model, x.val = x.val, n_p = n_p)
  grad.tcfs_par <- jacobian(TCF_normal, x = par_model, x.val = x.val, thresholds = opt_ctps, n_p = n_p)
  term3 <- hessian(func = function(par, x.val, n_p){
    CtP_normal(par = par[1:2], par_model = par[-c(1,2)], x.val = x.val, n_p = n_p)
  }, x = c(opt_ctps, par_model), x.val = x.val, n_p = n_p)
  term3.1 <- term3[1:2, 1:2]
  term4 <- jacobian(TCF_normal, x = opt_ctps, x.val = x.val, par = par_model, n_p = n_p)
  derv.cpt.para <- -solve(term3.1) %*% term3[1:2, -c(1,2)]
  derv.CtP.para <- term1 + term2 %*% derv.cpt.para
  derv.TCF.para <- grad.tcfs_par + term4 %*% derv.cpt.para
  vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  var_D3 <- as.numeric(derv.CtP.para %*% vcov_par_model %*% t(derv.CtP.para))
  vcov_cpts <- derv.cpt.para %*% vcov_par_model %*% t(derv.cpt.para)
  return(list(var_D3 = var_D3, vcov_cpts = vcov_cpts, vcov.tcfs = vcov.tcfs))
}

CtP_est_fun <- function(par_model, x.val, n_p, boxcox, par_boxcox, start = NULL, method.optim,
                        maxit, lower, upper){
  beta_d <- matrix(par_model[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  Z <- c(1, x.val)
  mu <- c(Z %*% beta_d[1,], Z %*% beta_d[2,], Z %*% beta_d[3,])
  sigma_d <- sqrt(par_model[(3*n_p + 2):length(par_model)]^2 + par_model[(3*n_p + 1)]^2)
  sigma2_d <- sigma_d^2
  if(is.null(start)){
    t1_0 <- mu[1] + qnorm(0.3)*sigma_d[1]
    t2_0 <- mu[3] + qnorm(1 - 0.4)*sigma_d[3]
    start <- c(t1_0, t2_0)
  }
  #
  CtP_optim <- optim(par = start, fn = CtP_normal, par_model = par_model, x.val = x.val,
                     n_p = n_p, method = method.optim, control = list(maxit = maxit),
                     lower = lower, upper = upper)
  if(CtP_optim$par[1] > CtP_optim$par[2]){
    cpt_CtP_est <- rev(CtP_optim$par)
  } else{
    cpt_CtP_est <- CtP_optim$par
  }
  if(boxcox) cpt_CtP_est <- boxcox_trans_back(cpt_CtP_est, par_boxcox)
  out <- cpt_CtP_est
  names(out) <- c("threshold_1", "threshold_2")
  return(out)
}


