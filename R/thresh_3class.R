####==========================================================================####
## This file consists of functions for estimating optimal pair of thresholds    ##
## Date: 26/03/2021																															##
####==========================================================================####

GYI_normal <- function(par, par_model, x.val, n_p){
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
  return(sum(res_tcfs))
}

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

MV_normal <- function(par, par_model, x.val, n_p){
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
  return(prod(res_tcfs))
}


#'@export
optThres3 <- function(out_lme2, x.val, start = NULL, method = "L-BFGS-B", maxit = 200,
                      lower = -Inf, upper = Inf){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  n_p <- out_lme2$n_p
  par_est <- out_lme2$est_para[1:(out_lme2$n_coef + 4)]
  if(is.null(start)){
    t1_0 <- c(1, x.val) %*% par_est[1:n_p] + qnorm(0.3)*sqrt(par_est[3*n_p + 1]^2 + par_est[3*n_p + 2]^2)
    t2_0 <- c(1, x.val) %*% par_est[(2*n_p + 1):(3*n_p)] +
      qnorm(1 - 0.4)*sqrt(par_est[3*n_p + 1]^2 + par_est[3*n_p + 4]^2)
    start = c(t1_0, t2_0)
  }
  ## GYI method
  GYI_optim <- optim(par = start, fn = GYI_normal, par_model = par_est, x.val = x.val, n_p = n_p,
                     method = method, control = list(maxit = maxit, fnscale = -1),
                     lower = lower, upper = upper)
  J3_est <- GYI_optim$value
  if(GYI_optim$par[1] > GYI_optim$par[2]){
    cpt_GYI_est <- rev(GYI_optim$par)
  } else{
    cpt_GYI_est <- GYI_optim$par
  }
  if(out_lme2$boxcox) cpt_GYI_est <- boxcox_trans_back(cpt_GYI_est,
                                                       out_lme2$est_para[length(out_lme2$est_para)])
  TCF_GYI_est <- TCF_normal(par = out_lme2$est_para, x.val = x.val, thresholds = cpt_GYI_est, n_p = n_p,
                            boxcox = out_lme2$boxcox)
  out_GYI <- c(J3_est, cpt_GYI_est, TCF_GYI_est)
  names(out_GYI) <- c("J3", "threshold_1", "threshold_2", "TCF_1", "TCF_2", "TCF_3")
  ## CtP method
  CtP_optim <- optim(par = start, fn = CtP_normal, par_model = par_est, x.val = x.val,
                     n_p = n_p, method = method, control = list(maxit = maxit), lower = lower, upper = upper)
  D3_est <- CtP_optim$value
  if(CtP_optim$par[1] > CtP_optim$par[2]){
    cpt_CtP_est <- rev(CtP_optim$par)
  } else{
    cpt_CtP_est <- CtP_optim$par
  }
  if(out_lme2$boxcox) cpt_CtP_est <- boxcox_trans_back(cpt_CtP_est,
                                                       out_lme2$est_para[length(out_lme2$est_para)])
  TCF_CtP_est <- TCF_normal(par = out_lme2$est_para, x.val = x.val, thresholds = cpt_CtP_est, n_p = n_p,
                            boxcox = out_lme2$boxcox)
  out_CtP <- c(D3_est, cpt_CtP_est, TCF_CtP_est)
  names(out_CtP) <- c("D3", "threshold_1", "threshold_2", "TCF_1", "TCF_2", "TCF_3")
  ## MV method
  MV_optim <- optim(par = start, fn = MV_normal, par_model = par_est, x.val = x.val, n_p = n_p,
                    method = method, control = list(maxit = maxit, fnscale = -1),
                    lower = lower, upper = upper)
  V3_est <- MV_optim$value
  if(MV_optim$par[1] > MV_optim$par[2]){
    cpt_MV_est <- rev(MV_optim$par)
  } else{
    cpt_MV_est <- MV_optim$par
  }
  if(out_lme2$boxcox) cpt_MV_est <- boxcox_trans_back(cpt_MV_est,
                                                       out_lme2$est_para[length(out_lme2$est_para)])
  TCF_MV_est <- TCF_normal(par = out_lme2$est_para, x.val = x.val, thresholds = cpt_MV_est, n_p = n_p,
                            boxcox = out_lme2$boxcox)
  out_MV <- c(V3_est, cpt_MV_est, TCF_MV_est)
  names(out_MV) <- c("V3", "threshold_1", "threshold_2", "TCF_1", "TCF_2", "TCF_3")
  #
  out <- list()
  out$GYI <- out_GYI
  out$CtP <- out_CtP
  out$MV <- out_MV
  return(out)
}




