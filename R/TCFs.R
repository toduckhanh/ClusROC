TCF_normal <- function(par, x.val, thresholds, n_p, boxcox = FALSE){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_d <- sqrt(par[(3*n_p + 2):(3*n_p + 4)]^2 + par[(3*n_p + 1)]^2)
  if(boxcox) thresholds <- boxcox_trans(thresholds, par[length(par)])
  Z <- c(1, x.val)
  res_tcfs <- numeric(3)
  res_tcfs[1] <- pnorm(thresholds[1], mean = Z %*% beta_d[1,], sd = sigma_d[1])
  res_tcfs[2] <- pnorm(thresholds[2], mean = Z %*% beta_d[2,], sd = sigma_d[2]) -
    pnorm(thresholds[1], mean = Z %*% beta_d[2,], sd = sigma_d[2])
  res_tcfs[3] <- 1 - pnorm(thresholds[2], mean = Z %*% beta_d[3,], sd = sigma_d[3])
  return(res_tcfs)
}

TCF_normal_vcov <- function(par_model, x.val, thresholds, vcov_par_model, n_p, fixed = FALSE, boxcox = FALSE,
                            type_thresholds = c("GYI", "CtP", "MV")){
  type_thresholds <- match.arg(type_thresholds)
  # if(boxcox) thresholds <- boxcox_trans(thresholds, par_model[length(par_model)])
  if(fixed){
    grad.tcfs <- jacobian(TCF_normal, x = par_model, x.val = x.val, thresholds = thresholds, n_p = n_p,
                          boxcox = boxcox)
    vcov.tcfs <- grad.tcfs %*% vcov_par_model %*% t(grad.tcfs)
  } else{
    grad.tcfs_par <- jacobian(TCF_normal, x = par_model, x.val = x.val, thresholds = thresholds, n_p = n_p,
                              boxcox = boxcox)
    grad.tcfs_thresholds <- jacobian(TCF_normal, x = thresholds, x.val = x.val, par = par_model, n_p = n_p)
    term3 <- switch (type_thresholds,
                     GYI = hessian(func = function(par, x.val, n_p){
                       GYI_normal(par = par[1:2], par_model = par[-c(1,2)], x.val = x.val, n_p = n_p)
                     }, x = c(thresholds, par_model), x.val = x.val, n_p = n_p),
                     CtP = hessian(func = function(par, x.val, n_p){
                       CtP_normal(par = par[1:2], par_model = par[-c(1,2)], x.val = x.val, n_p = n_p)
                     }, x = c(thresholds, par_model), x.val = x.val, n_p = n_p),
                     MV = hessian(func = function(par, x.val, n_p){
                       MV_normal(par = par[1:2], par_model = par[-c(1,2)], x.val = x.val, n_p = n_p)
                     }, x = c(thresholds, par_model), x.val = x.val, n_p = n_p)
    )
    term3.1 <- term3[1:2, 1:2]
    derv.cpt.para <- -solve(term3.1) %*% term3[1:2, -c(1,2)]
    derv.TCF.para <- grad.tcfs_par + grad.tcfs_thresholds %*% derv.cpt.para
    vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  }
  return(vcov.tcfs)
}

