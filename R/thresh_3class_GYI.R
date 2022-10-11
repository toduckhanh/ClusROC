####========================================================================####
## This file consists of functions for estimating optimal pair of thresholds  ##
## Based on gyi method                                                        ##
## Date: 10/10/2022																														##
####========================================================================####

gyi_normal <- function(par, par_model, z, n_p) {
  if (par[1] > par[2]) {
    x <- par[1]
    par[1] <- par[2]
    par[2] <- x
  }
  beta_d <- par_model[1:(3 * n_p)]
  sigma_d <- sqrt(par_model[(3 * n_p + 2):length(par_model)]^2 +
                    par_model[(3 * n_p + 1)]^2)
  mu_est <- z %*% beta_d
  res_tcfs <- numeric(3)
  res_tcfs[1] <- pnorm(par[1], mean = mu_est[1], sd = sigma_d[1])
  res_tcfs[2] <- pnorm(par[2], mean = mu_est[2], sd = sigma_d[2]) -
    pnorm(par[1], mean = mu_est[2], sd = sigma_d[2])
  res_tcfs[3] <- 1 - pnorm(par[2], mean = mu_est[3], sd = sigma_d[3])
  return(sum(res_tcfs))
}

gyi_normal_var <- function(opt_ctps, par_model, vcov_par_model, z, n_p) {
  term1 <- grad(gyi_normal, x = par_model, par = opt_ctps, z = z, n_p = n_p)
  term2 <- grad(gyi_normal, x = opt_ctps, par_model = par_model, z = z,
                n_p = n_p)
  grad_tcfs_par <- jacobian(tcf_normal, x = par_model, z = z,
                            thresholds = opt_ctps, n_p = n_p)
  term3 <- hessian(func = function(par, z, n_p) {
    gyi_normal(par = par[1:2], par_model = par[-c(1, 2)], z = z, n_p = n_p)
  }, x = c(opt_ctps, par_model), z = z, n_p = n_p)
  term3_1 <- term3[1:2, 1:2]
  term4 <- jacobian(tcf_normal, x = opt_ctps, z = z, par = par_model, n_p = n_p)
  derv_cpt_para <- -solve(term3_1) %*% term3[1:2, -c(1, 2)]
  derv_gyi_para <- term1 + term2 %*% derv_cpt_para
  derv_tcf_para <- grad_tcfs_par + term4 %*% derv_cpt_para
  vcov_tcfs <- derv_tcf_para %*% vcov_par_model %*% t(derv_tcf_para)
  var_j3 <- as.numeric(derv_gyi_para %*% vcov_par_model %*% t(derv_gyi_para))
  vcov_cpts <- derv_cpt_para %*% vcov_par_model %*% t(derv_cpt_para)
  return(list(var_j3 = var_j3, vcov_cpts = vcov_cpts, vcov_tcfs = vcov_tcfs))
}

gyi_est_fun <- function(par_model, z, n_p, boxcox, par_boxcox,
                        crit_var = 0.01) {
  beta_d <- par_model[1:(3 * n_p)]
  mu <- z %*% beta_d
  sigma_d <- sqrt(par_model[(3 * n_p + 2):length(par_model)]^2 +
                    par_model[(3 * n_p + 1)]^2)
  sigma2_d <- sigma_d^2
  #
  var1_check <- abs(sigma2_d[1] / sigma2_d[2] - 1) < crit_var
  var2_check <- abs(sigma2_d[2] / sigma2_d[3] - 1) < crit_var
  if (var1_check) {
    cpt_1 <- (mu[1] + mu[2]) / 2
  } else {
    a1 <- mu[2] * sigma2_d[1] - mu[1] * sigma2_d[2]
    b1 <- sigma_d[1] * sigma_d[2]
    c1 <- (mu[1] - mu[2])^2
    d1 <- (sigma2_d[1] - sigma2_d[2])
    cpt_1 <- (a1 - b1 * sqrt(c1 + d1 * log(sigma2_d[1] / sigma2_d[2]))) / d1
  }
  if (var2_check) {
    cpt_2 <- (mu[2] + mu[3]) / 2
  } else {
    a2 <- (mu[3] * sigma2_d[2] - mu[2] * sigma2_d[3])
    b2 <- sigma_d[2] * sigma_d[3]
    c2 <- (mu[2] - mu[3])^2
    d2 <- (sigma2_d[2] - sigma2_d[3])
    cpt_2 <- (a2 - b2 * sqrt(c2 + d2 * log(sigma2_d[2] / sigma2_d[3]))) / d2
  }
  out <- c(cpt_1, cpt_2)
  if (boxcox) {
    cpt_1_org <- boxcox_trans_back(cpt_1, par_boxcox)
    cpt_2_org <- boxcox_trans_back(cpt_2, par_boxcox)
    out <- c(cpt_1_org, cpt_2_org)
  }
  names(out) <- c("threshold_1", "threshold_2")
  return(out)
}
