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

GYI_normal_var <- function(opt_ctps, par_model, vcov_par_model, x.val, n_p){
  term1 <- grad(GYI_normal, x = par_model, par = opt_ctps, x.val = x.val, n_p = n_p)
  term2 <- grad(GYI_normal, x = opt_ctps, par_model = par_model, x.val = x.val, n_p = n_p)
  grad.tcfs_par <- jacobian(TCF_normal, x = par_model, x.val = x.val, thresholds = opt_ctps, n_p = n_p)
  term3 <- hessian(func = function(par, x.val, n_p){
    GYI_normal(par = par[1:2], par_model = par[-c(1,2)], x.val = x.val, n_p = n_p)
  }, x = c(opt_ctps, par_model), x.val = x.val, n_p = n_p)
  term3.1 <- term3[1:2, 1:2]
  term4 <- jacobian(TCF_normal, x = opt_ctps, x.val = x.val, par = par_model, n_p = n_p)
  derv.cpt.para <- -solve(term3.1) %*% term3[1:2, -c(1,2)]
  derv.GYI.para <- term1 + term2 %*% derv.cpt.para
  derv.TCF.para <- grad.tcfs_par + term4 %*% derv.cpt.para
  vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  var_J3 <- as.numeric(derv.GYI.para %*% vcov_par_model %*% t(derv.GYI.para))
  vcov_cpts <- derv.cpt.para %*% vcov_par_model %*% t(derv.cpt.para)
  return(list(var_J3 = var_J3, vcov_cpts = vcov_cpts, vcov.tcfs = vcov.tcfs))
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

MV_normal_var <- function(opt_ctps, par_model, vcov_par_model, x.val, n_p){
  term1 <- grad(MV_normal, x = par_model, par = opt_ctps, x.val = x.val, n_p = n_p)
  term2 <- grad(MV_normal, x = opt_ctps, par_model = par_model, x.val = x.val, n_p = n_p)
  grad.tcfs_par <- jacobian(TCF_normal, x = par_model, x.val = x.val, thresholds = opt_ctps, n_p = n_p)
  term3 <- hessian(func = function(par, x.val, n_p){
    MV_normal(par = par[1:2], par_model = par[-c(1,2)], x.val = x.val, n_p = n_p)
  }, x = c(opt_ctps, par_model), x.val = x.val, n_p = n_p)
  term3.1 <- term3[1:2, 1:2]
  term4 <- jacobian(TCF_normal, x = opt_ctps, x.val = x.val, par = par_model, n_p = n_p)
  derv.cpt.para <- -solve(term3.1) %*% term3[1:2, -c(1,2)]
  derv.MV.para <- term1 + term2 %*% derv.cpt.para
  derv.TCF.para <- grad.tcfs_par + term4 %*% derv.cpt.para
  vcov.tcfs <- derv.TCF.para %*% vcov_par_model %*% t(derv.TCF.para)
  var_V3 <- as.numeric(derv.MV.para %*% vcov_par_model %*% t(derv.MV.para))
  vcov_cpts <- derv.cpt.para %*% vcov_par_model %*% t(derv.cpt.para)
  return(list(var_V3 = var_V3, vcov_cpts = vcov_cpts, vcov.tcfs = vcov.tcfs))
}

optThres3_core <- function(para, x.val, n_p, n_coef, boxcox, start = NULL, method = "L-BFGS-B", maxit = 200,
                           lower = -Inf, upper = Inf){
  par_est <- para[1:(n_coef + 4)]
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
  if(boxcox) cpt_GYI_est <- boxcox_trans_back(cpt_GYI_est, para[length(para)])
  TCF_GYI_est <- TCF_normal(par = para, x.val = x.val, thresholds = cpt_GYI_est, n_p = n_p, boxcox = boxcox)
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
  if(boxcox) cpt_CtP_est <- boxcox_trans_back(cpt_CtP_est, para[length(para)])
  TCF_CtP_est <- TCF_normal(par = para, x.val = x.val, thresholds = cpt_CtP_est, n_p = n_p,
                            boxcox = boxcox)
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
  if(boxcox) cpt_MV_est <- boxcox_trans_back(cpt_MV_est, para[length(para)])
  TCF_MV_est <- TCF_normal(par = para, x.val = x.val, thresholds = cpt_MV_est, n_p = n_p, boxcox = boxcox)
  out_MV <- c(V3_est, cpt_MV_est, TCF_MV_est)
  names(out_MV) <- c("V3", "threshold_1", "threshold_2", "TCF_1", "TCF_2", "TCF_3")
  #
  out <- list()
  out$GYI <- out_GYI
  out$CtP <- out_CtP
  out$MV <- out_MV
  return(out)
}


#'@export
optThres3 <- function(out_lme2, x.val, apVar = TRUE, bootstrap = FALSE,
                      nR = 1000, type.boot = c("cluster", "stratified"), data,
                      parallel = FALSE, ncpus = ifelse(parallel, 2, NULL),
                      start = NULL, method = "L-BFGS-B", maxit = 200, lower = -Inf, upper = Inf){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  n_p <- out_lme2$n_p
  n_coef <- out_lme2$n_coef
  fit <- list()
  if(inherits(x.val, "numeric") & n_p == 2){ # 1 covariate
    fit$x.val <- x.val
    fit$n_p <- n_p
    n_x <- length(x.val)
    temp_thres <- list()
    for(i in 1:n_x){
      out <- optThres3_core(para = out_lme2$est_para, x.val = x.val[i], n_p = n_p, n_coef = n_coef,
                            boxcox = out_lme2$boxcox, start = start, method = method, maxit = maxit,
                            lower = lower, upper = upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) x[c("threshold_1", "threshold_2")])),
                                    Method = factor(c("Generalized Youden Index", "Closest to Perfection",
                                                      "Max Volume"),
                                                    levels = c("Generalized Youden Index",
                                                               "Closest to Perfection", "Max Volume")),
                                    x = x.val[i], row.names = NULL)
    }
    fit$thres3 <- do.call(rbind, temp_thres)
    if(apVar){
      fit$vcov.thres3 <- list()
      if(!out_lme2$boxcox & !bootstrap){ ## asymptotic variance under normal distribution
        for(i in 1:n_x){
          fit$vcov.thres3$GYI[[i]] <- GYI_normal_var(opt_ctps = as.numeric(temp_thres[[i]][1, 1:2]),
                                                     par_model = out_lme2$est_para,
                                                     vcov_par_model = out_lme2$vcov_sand, x.val = x.val[i],
                                                     n_p = n_p)$vcov_cpts
          fit$vcov.thres3$CtP[[i]] <- CtP_normal_var(opt_ctps = as.numeric(temp_thres[[i]][2, 1:2]),
                                                     par_model = out_lme2$est_para,
                                                     vcov_par_model = out_lme2$vcov_sand, x.val = x.val[i],
                                                     n_p = n_p)$vcov_cpts
          fit$vcov.thres3$MV[[i]] <- MV_normal_var(opt_ctps = as.numeric(temp_thres[[i]][3, 1:2]),
                                                   par_model = out_lme2$est_para,
                                                   vcov_par_model = out_lme2$vcov_sand, x.val = x.val[i],
                                                   n_p = n_p)$vcov_cpts
        }
      } else{
        for(i in 1:n_x){
          out_bts <- boot_lme2(B = nR, name.test = out_lme2$name.test, name.class = out_lme2$name.class,
                               name.covars = out_lme2$name.covars, name.clust = out_lme2$name.clust,
                               data = data, x.val = x.val[i], type = type.boot, boxcox = out_lme2$boxcox,
                               parallel = parallel, ncpus = ncpus)
          thres_GYI_bts <- thres_CtP_bts <- thres_MV_bts <- matrix(0, nrow = 2, ncol = nR)
          for(k in 1:nR){
            temp <- optThres3_core(para = out_bts[,k], x.val = x.val[i], n_p = n_p, n_coef = n_coef,
                                   boxcox = out_lme2$boxcox)
            thres_GYI_bts[,k] <- temp$GYI[c("threshold_1", "threshold_2")]
            thres_CtP_bts[,k] <- temp$CtP[c("threshold_1", "threshold_2")]
            thres_MV_bts[,k] <- temp$MV[c("threshold_1", "threshold_2")]
          }
          fit$vcov.thres3$GYI[[i]] <- cov(t(thres_GYI_bts))
          fit$vcov.thres3$CtP[[i]] <- cov(t(thres_CtP_bts))
          fit$vcov.thres3$MV[[i]] <- cov(t(thres_MV_bts))
        }
      }
    }
  }
  class(fit) <- "optThres3"
  return(fit)
}

#'@export
plot.optThres3 <- function(x, plot.ellip = FALSE, ci.level = ifelse(plot.ellip, 0.95, NULL),
                           colors = NULL, xlims, ylims, file.name = NULL){
  if(isFALSE(inherits(x, "optThres3"))) stop("x was not from optThres3()!")
  if(inherits(x$x.val, "numeric")& x$n_p == 2){
    n_x <- length(x$x.val)
    dt_thres <- x$thres3[, 1:3]
    colnames(dt_thres) <- c("x", "y", "method")
    dt_thres$pts <- as.factor(rep(1:n_x, each = 3))
    dt_thres_list <- split(dt_thres, dt_thres$method)
    dt_ell_thres_GYI <- dt_ell_thres_CtP <- dt_ell_thres_MV <- data.frame()
    for(i in 1:n_x){
      uu <- ellipse(center = as.numeric(dt_thres_list$`Generalized Youden Index`[i, 1:2]),
                    shape = x$vcov.thres3$GYI[[i]],
                    radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_GYI <- rbind(dt_ell_thres_GYI, data.frame(uu, pts = as.factor(i)))
      uu1 <- ellipse(center = as.numeric(dt_thres_list$`Closest to Perfection`[i, 1:2]),
                     shape = x$vcov.thres3$CtP[[i]],
                     radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_CtP <- rbind(dt_ell_thres_CtP, data.frame(uu1, pts = as.factor(i)))
      uu2 <- ellipse(center = as.numeric(dt_thres_list$`Max Volume`[i, 1:2]),
                     shape = x$vcov.thres3$MV[[i]],
                     radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_MV <- rbind(dt_ell_thres_MV, data.frame(uu2, pts = as.factor(i)))
    }
    dt_ell_thres <- rbind(dt_ell_thres_GYI, dt_ell_thres_CtP, dt_ell_thres_MV)
    dt_ell_thres$method <- factor(rep(c("Generalized Youden Index", "Closest to Perfection", "Max Volume"),
                                      each = 52*n_x),
                                  levels = c("Generalized Youden Index", "Closest to Perfection", "Max Volume"))
    pp <- ggplot(dt_thres, aes(x = x, y = y, colour = pts)) +
      facet_grid( ~ method) +
      geom_point(size = 0.5) +
      geom_path(data = dt_ell_thres, aes(x = x, y = y, group = pts), size = 0.5) +
      xlab("Optimal threshold 1") + ylab("Optimal threshold 2") +
      theme_bw() +
      theme(legend.position = "bottom")
    if(is.null(colors)){
      pp <- pp + scale_color_manual(name = "Value(s) of covariate:", breaks = as.character(1:n_x),
                                    values = topo.colors(n_x),
                                    labels = x$x.val)
    } else{
      if(length(colors) != n_x) stop(paste("Number of colors need to be equal:", n_x))
      else{
        pp <- pp + scale_color_manual(name = "Value(s) of covariate:", breaks = as.character(1:n_x),
                                      values = colors, labels = x$x.val)
      }
    }
    if(!missing(xlims) & !missing(ylims)){
      pp <- pp + xlim(xlims) + ylim(ylims)
    }
  }
  if(!is.null(file.name)){
    ggexport(pp, filename = file.name)
  }
  pp
}

