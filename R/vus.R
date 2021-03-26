####==========================================================================####
## This file consists of functions for estimating covariate-specific VUS        ##
## Date: 26/03/2021																															##
####==========================================================================####

#' @import utils
#' @import numDeriv
#' @import stats

theta_1 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f01 <- function(s){
    dnorm((s - mu_2)/sigma_e[2])*pnorm((s - mu_1)/sigma_e[1])*
      pnorm((s - mu_3)/sigma_e[3], lower.tail = FALSE)/sigma_e[2]
  }
  return(integrate(f01, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

theta_2 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f02 <- function(s){
    dnorm((s - mu_2)/sigma_e[2])*pnorm((s - mu_1)/sigma_e[1])*pnorm((s - mu_3)/sqrt(2*sigma_c^2 + sigma_e[3]^2),
                                                                    lower.tail = FALSE)/sigma_e[2]
  }
  return(integrate(f02, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

theta_3 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f03 <- function(s){
    dnorm((s - mu_2)/sqrt(2*sigma_c^2 + sigma_e[2]^2))*pnorm((s - mu_1)/sigma_e[1])*
      pnorm((s - mu_3)/sigma_e[3], lower.tail = FALSE)/sqrt(2*sigma_c^2 + sigma_e[2]^2)
  }
  return(integrate(f03, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

theta_4 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f04 <- function(s){
    dnorm((s - mu_2)/sigma_e[2])*pnorm((s - mu_1)/sqrt(2*sigma_c^2 + sigma_e[1]^2))*
      pnorm((s - mu_3)/sigma_e[3], lower.tail = FALSE)/sigma_e[2]
  }
  return(integrate(f04, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

theta_5 <- function(par, x.val, n_p, n_c, n_k, n, subdivisions = 1000, ...){
  beta_d <- matrix(par[1:(3*n_p)], ncol = n_p, nrow = 3, byrow = TRUE)
  sigma_e <- par[(3*n_p + 2):length(par)]
  sigma_c <- par[(3*n_p + 1)]
  Z <- c(1, x.val)
  mu_1 <- as.numeric(Z %*% beta_d[1,])
  mu_2 <- as.numeric(Z %*% beta_d[2,])
  mu_3 <- as.numeric(Z %*% beta_d[3,])
  f05 <- function(s){
    dnorm((s - mu_2)/sqrt(sigma_c^2 + sigma_e[2]^2))*pnorm((s - mu_1)/sqrt(sigma_c^2 + sigma_e[1]^2))*
      pnorm((s - mu_3)/sqrt(sigma_c^2 + sigma_e[3]^2), lower.tail = FALSE)/sqrt(sigma_c^2 + sigma_e[2]^2)
  }
  return(integrate(f05, lower = -Inf, upper = Inf, subdivisions = subdivisions, ...)$value)
}

#' @export
VUS <- function(out_lme2, x.val, apVar = FALSE, ci.level = 0.95, subdivisions = 1000, ...){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  n_p <- out_lme2$n_p
  n_c <- out_lme2$cls
  n_k <- out_lme2$n_c
  n <- out_lme2$n
  par <- out_lme2$est_para[1:(out_lme2$n_coef + 4)]
  vcov_par_model <- out_lme2$vcov_sand[1:(out_lme2$n_coef + 4), 1:(out_lme2$n_coef + 4)]
  theta_est_1 <- theta_1(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  theta_est_2 <- theta_2(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  theta_est_3 <- theta_3(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  theta_est_4 <- theta_4(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  theta_est_5 <- theta_5(par = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n,
                         subdivisions = subdivisions, ...)
  ## check if all clusters/families have the same cluster size, if yes, the weight calculation_k is simplified
  unique.n_k <- unique(n_k)
  equal.n_k.flag <- ifelse(length(unique.n_k) == 1, TRUE, FALSE)
  if(equal.n_k.flag){
    n_k.div.N <- (unique.n_k/n)^3;
    p.sss <- n_c*n_k.div.N;
    p.ssk <- n_c*(n_c - 1)*n_k.div.N
    p.ijk <- n_c*(n_c - 1)*(n_c - 2)*n_k.div.N
  } else{
    n_k3 <- n_k^3
    p.sss <- sum(n_k3)/n^3
    n_k.sqr <- n_k^2
    outer.n_kSQR.n_k <- outer(n_k.sqr, n_k)
    diag(outer.n_kSQR.n_k) <- NA
    p.ssk <- sum(outer.n_kSQR.n_k, na.rm = TRUE)/(n^3)
    combo <- combn(x = 1:n_c, m = 3)
    p.ijk <- sum(apply(combo, 2, function(idx) prod(n_k[idx])))/(n^3)
    p.ijk <- 6*p.ijk
  }
  p.skk <- p.sks <- p.ssk
  vus_est <- theta_est_1*p.sss + theta_est_2*p.ssk + theta_est_3*p.sks + theta_est_4*p.skk + theta_est_5*p.ijk
  if(isTRUE(apVar) & !all(is.na(vcov_par_model))){
    jac_vus <- rbind(jacobian(theta_1, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n),
                     jacobian(theta_2, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n),
                     jacobian(theta_3, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n),
                     jacobian(theta_4, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n),
                     jacobian(theta_5, x = par, x.val = x.val, n_p = n_p, n_c = n_c, n_k = n_k, n = n))
    Sig_1 <- jac_vus %*% vcov_par_model %*% t(jac_vus)
    pp <- c(p.sss, p.ssk, p.sks, p.skk, p.ijk)
    vus_sd <- as.numeric(sqrt(pp %*% Sig_1 %*% pp))
    vus_ci <- vus_est + c(-1,1)*qnorm((1 + ci.level)/2)*vus_sd
    out <- c(vus_est, vus_sd, vus_ci)
    names(out) <- c("Est", "SE", "Low.ci", "Up.ci")
    return(out)
  } else return(vus_est)
}

