####==========================================================================####
## This file consists of functions for estimating optimal pair of thresholds    ##
## Date: 25/05/2021																															##
####==========================================================================####

optThres3_core <- function(method, para, x.val, n_p, n_coef, boxcox, start, method.optim,
                           maxit, lower, upper){
  par_model <- para[1:(n_coef + 4)]
  if(boxcox) par_boxcox <- para[length(para)]
  out <- list()
  if("GYI" %in% method){
    out$GYI <- GYI_est_fun(par_model = par_model, x.val = x.val, n_p = n_p, boxcox = boxcox,
                           par_boxcox = par_boxcox, crit.var = 0.01)
  }
  if("CtP" %in% method){
    out$CtP <- CtP_est_fun(par_model = par_model, x.val = x.val, n_p = n_p, boxcox = boxcox,
                           par_boxcox = par_boxcox, start = start, method.optim = method.optim,
                           maxit = maxit, lower = lower, upper = upper)
  }
  if("MV" %in% method){
    out$MV <- MV_est_fun(par_model = par_model, x.val = x.val, n_p = n_p, boxcox = boxcox,
                         par_boxcox = par_boxcox, start = start, method.optim = method.optim,
                         maxit = maxit, lower = lower, upper = upper)
  }
  return(out)
}

optThres3_se <- function(method, thres_est, out_lme2, x.val, n_p, n_coef, bootstrap, nR, type.boot, data,
                         parallel, ncpus, start, method.optim, maxit, lower, upper){
  out <- list()
  ## asymptotic variance under normal distribution
  if(!out_lme2$boxcox & !bootstrap){
    if("GYI" %in% method){
      out$cov_GYI <- GYI_normal_var(
        opt_ctps = as.numeric(thres_est[thres_est$Method == "Generalized Youden Index Approach", 1:2]),
        par_model = out_lme2$est_para,
        vcov_par_model = out_lme2$vcov_sand, x.val = x.val, n_p = n_p)$vcov_cpts
    }
    if("CtP" %in% method){
      out$cov_CtP <- CtP_normal_var(
        opt_ctps = as.numeric(thres_est[thres_est$Method == "Closest to Perfection Approach", 1:2]),
        par_model = out_lme2$est_para,
        vcov_par_model = out_lme2$vcov_sand, x.val = x.val, n_p = n_p)$vcov_cpts
    }
    if("MV" %in% method){
      out$cov_MV <- MV_normal_var(
        opt_ctps = as.numeric(thres_est[thres_est$Method == "Max Volume Approach", 1:2]),
        par_model = out_lme2$est_para,
        vcov_par_model = out_lme2$vcov_sand, x.val = x.val, n_p = n_p)$vcov_cpts
    }
  } else{ ## bootstrap
    out_bts <- boot_lme2(B = nR, name.test = out_lme2$name.test, name.class = out_lme2$name.class,
                         name.covars = out_lme2$name.covars, name.clust = out_lme2$name.clust,
                         data = data, x.val = x.val, type = type.boot, boxcox = out_lme2$boxcox,
                         parallel = parallel, ncpus = ncpus)
    if("GYI" %in% method){
      thres_GYI_bts <- matrix(0, nrow = 2, ncol = nR)
      for(k in 1:nR){
        thres_GYI_bts[,k] <- optThres3_core(method = "GYI", para = out_bts[,k], x.val = x.val, n_p = n_p,
                                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                            method.optim = method.optim, maxit = maxit, lower = lower,
                                            upper = upper)$GYI
      }
      out$cov_GYI <- cov(na.omit(t(thres_GYI_bts)))
    }
    if("CtP" %in% method){
      thres_CtP_bts <- matrix(0, nrow = 2, ncol = nR)
      for(k in 1:nR){
        thres_CtP_bts[,k] <- optThres3_core(method = "CtP", para = out_bts[,k], x.val = x.val, n_p = n_p,
                                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                            method.optim = method.optim, maxit = maxit, lower = lower,
                                            upper = upper)$CtP
      }
      out$cov_CtP <- cov(na.omit(t(thres_CtP_bts)))
    }
    if("MV" %in% method){
      thres_MV_bts <- matrix(0, nrow = 2, ncol = nR)
      for(k in 1:nR){
        thres_MV_bts[,k] <- optThres3_core(method = "MV", para = out_bts[,k], x.val = x.val, n_p = n_p,
                                           n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                           method.optim = method.optim, maxit = maxit, lower = lower,
                                           upper = upper)$MV
      }
      out$cov_MV <- cov(na.omit(t(thres_MV_bts)))
    }
  }
  return(out)
}

#'@export
optThres3control <- function(method.optim = c("L-BFGS-B", "BFGS", "Nelder-Mead", "CG", "SANN", "Brent"),
                             start = NULL, maxit = 200, lower = -Inf, upper = Inf, nR = 250,
                             type.boot = c("cluster", "stratified"),
                             parallel = FALSE, ncpus = NULL){
  if(parallel & is.null(ncpus)) ncpus <- 2
  list(method.optim = match.arg(method.optim), start = start, maxit = maxit, lower = lower, upper = upper,
       nR = nR, type.boot = match.arg(type.boot), parallel = parallel, ncpus = ncpus)
}

#'@export
optThres3 <- function(method = c("GYI", "CtP", "MV"), out_lme2, x.val, apVar = TRUE,
                      data, control = list()){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  if(out_lme2$n_coef/out_lme2$n_p != 3) stop("There is not a case of three-class setting!")
  if(missing(method)) stop("Please, input the selection method(s)!")
  controlvals <- optThres3control()
  if (!missing(control)) {
    controlvals[names(control)] <- control
  }
  ## checking the method
  methodtemp <- substitute(me, list(me = method))
  okMethod <- c("GYI", "CtP", "MV")
  if(length(methodtemp) > 3) stop(gettextf("The maximum number of selection methods are 3!"), domain = NA)
  if(any(duplicated(methodtemp))) stop(gettextf("The selection methods need to be unique!"), domain = NA)
  if(!is.character(methodtemp)) methodtemp <- deparse(methodtemp)
  if(any(sapply(methodtemp, function(x) !is.element(x, okMethod)))){
    stop(gettextf("the selection method(s) should be: %s", paste(sQuote(okMethod), collapse = ", ")),
         domain = NA)
  }
  methodtemp <- methodtemp[c(which(methodtemp == "GYI"), which(methodtemp == "CtP"), which(methodtemp == "MV"))]
  Method <- character(length(methodtemp))
  Method[methodtemp == "GYI"] <- "Generalized Youden Index Approach"
  Method[methodtemp == "CtP"] <- "Closest to Perfection Approach"
  Method[methodtemp == "MV"] <- "Max Volume Approach"
  Method <- factor(Method, levels = c("Generalized Youden Index Approach",
                                      "Closest to Perfection Approach", "Max Volume Approach"))
  n_p <- out_lme2$n_p
  n_coef <- out_lme2$n_coef
  fit <- list()
  fit$method <- methodtemp
  if(n_p == 1){# no covariate
    if(!missing(x.val)) {
      if(!is.null(x.val)) warning("Sepecified value(s) of covariate(s) are not used!", call. = FALSE)
    }
    x.val <- NULL
    fit$x.val <- x.val
    fit$n_p <- n_p
    n_x <- 1
    temp_thres <- list()
    for(i in 1:n_x){
      out <- optThres3_core(method = methodtemp, para = out_lme2$est_para, x.val = x.val[i], n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) x[c("threshold_1", "threshold_2")])),
                                    Method = Method, row.names = NULL)
    }
    if(apVar){
      if(!out_lme2$boxcox){ ## asymptotic variance under normal distribution
        bootstrap <- FALSE
        if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing!")
      } else {
        bootstrap <- TRUE
        if(missing(data)) stop("The original data is required to process bootstrap procedure!")
      }
      fit$vcov.thres3 <- list()
      for(i in 1:n_x){
        fit$vcov.thres3[[i]] <- optThres3_se(method = methodtemp, thres_est = temp_thres[[i]],
                                             out_lme2 = out_lme2, x.val = x.val[i], n_p = n_p,
                                             n_coef = n_coef, bootstrap = bootstrap, nR = controlvals$nR,
                                             type.boot = controlvals$type.boot, data = data,
                                             parallel = controlvals$parallel, ncpus = controlvals$ncpus,
                                             start = controlvals$start, method.optim = controlvals$method.optim,
                                             maxit = controlvals$maxit, lower = controlvals$lower,
                                             upper = controlvals$upper)
      }
    }
  }
  if(n_p == 2){
    if(missing(x.val)) stop("Please input specific value(s) of covariate.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariate.")
    if(!inherits(x.val, "numeric")) stop("For case of 1 covariate, please input a number or a vector.")
    if(any(is.na(x.val))) stop("NA value(s) not allowed!")
    fit$x.val <- x.val
    fit$n_p <- n_p
    n_x <- length(x.val)
    temp_thres <- list()
    for(i in 1:n_x){
      out <- optThres3_core(method = methodtemp, para = out_lme2$est_para, x.val = x.val[i], n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) x[c("threshold_1", "threshold_2")])),
                                    Method = Method, x = x.val[i], row.names = NULL)
    }
    if(apVar){
      if(!out_lme2$boxcox){ ## asymptotic variance under normal distribution
        bootstrap <- FALSE
        if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing!")
      } else {
        bootstrap <- TRUE
        if(missing(data)) stop("The original data is required to process bootstrap procedure!")
      }
      fit$vcov.thres3 <- list()
      for(i in 1:n_x){
        fit$vcov.thres3[[i]] <- optThres3_se(method = methodtemp, thres_est = temp_thres[[i]],
                                             out_lme2 = out_lme2, x.val = x.val[i], n_p = n_p,
                                             n_coef = n_coef, bootstrap = bootstrap, nR = controlvals$nR,
                                             type.boot = controlvals$type.boot, data = data,
                                             parallel = controlvals$parallel, ncpus = controlvals$ncpus,
                                             start = controlvals$start, method.optim = controlvals$method.optim,
                                             maxit = controlvals$maxit, lower = controlvals$lower,
                                             upper = controlvals$upper)
      }
    }
  }
  if(n_p > 2){ # multiple covariates
    if(missing(x.val)) stop("Please input specific value(s) of covariates.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariates.")
    if(inherits(x.val, "numeric")){
      if(length(x.val) != n_p - 1) stop(paste("For case of", n_p - 1, "covariates, please input a vector of", n_p - 1, "values of covariates."))
    }
    if(inherits(x.val, "matrix")) {
      if(ncol(x.val) != n_p - 1) stop(paste("For case of m points of", n_p - 1, "covariates, please input a matrix with", n_p - 1, "columns and m rows containing values of covariates."))
    }
    if(any(is.na(x.val))) stop("NA value(s) not allowed!")
    x.val <- matrix(x.val, ncol = n_p - 1, byrow = FALSE)
    fit$x.val <- x.val
    fit$n_p <- n_p
    n_x <- nrow(x.val)
    temp_thres <- list()
    for(i in 1:n_x){
      out <- optThres3_core(method = methodtemp, para = out_lme2$est_para, x.val = x.val[i,], n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) x[c("threshold_1", "threshold_2")])),
                                    Method = Method, x = t(x.val[i,]), row.names = NULL)
    }
    if(apVar){
      if(!out_lme2$boxcox){ ## asymptotic variance under normal distribution
        bootstrap <- FALSE
        if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing!")
      } else {
        bootstrap <- TRUE
        if(missing(data)) stop("The original data is required to process bootstrap procedure!")
      }
      fit$vcov.thres3 <- list()
      for(i in 1:n_x){
        fit$vcov.thres3[[i]] <- optThres3_se(method = methodtemp, thres_est = temp_thres[[i]],
                                             out_lme2 = out_lme2, x.val = x.val[i,], n_p = n_p,
                                             n_coef = n_coef, bootstrap = bootstrap, nR = controlvals$nR,
                                             type.boot = controlvals$type.boot, data = data,
                                             parallel = controlvals$parallel, ncpus = controlvals$ncpus,
                                             start = controlvals$start, method.optim = controlvals$method.optim,
                                             maxit = controlvals$maxit, lower = controlvals$lower,
                                             upper = controlvals$upper)
      }
    }
  }
  fit$thres3 <- do.call(rbind, temp_thres)
  class(fit) <- "optThres3"
  return(fit)
}

#'@export
plot.optThres3 <- function(x, ci.level = 0.95, colors = NULL, xlims, ylims, size.point = 0.5,
                           size.path = 0.5, names.labels, file.name = NULL, ...){
  if(isFALSE(inherits(x, "optThres3"))) stop("x was not from optThres3()!")
  if(x$n_p == 1){
    n_x <- 1
    labels <- "Intercept"
    if(missing(names.labels)) names.labels <- " "
  }
  if(x$n_p == 2) {
    n_x <- length(x$x.val)
    labels <- x$x.val
    if(missing(names.labels)) names.labels <- "Value(s) of covariate:"
  }
  if(x$n_p == 3) {
    n_x <- nrow(x$x.val)
    labels <- apply(x$x.val, 1, function(y) paste0("(", paste(y, collapse = ", "), ")"))
    if(missing(names.labels)) names.labels <- "Value(s) of covariates:"
  }
  dt_thres <- x$thres3[, 1:3]
  colnames(dt_thres) <- c("x", "y", "method")
  dt_thres$pts <- as.factor(rep(1:n_x, each = length(x$method)))
  dt_thres_list <- split(dt_thres, dt_thres$method)
  dt_ell_thres <- list()
  if("GYI" %in% x$method){
    dt_ell_thres_GYI <- data.frame()
    for(i in 1:n_x){
      uu <- ellipse(center = as.numeric(dt_thres_list$`Generalized Youden Index Approach`[i, 1:2]),
                    shape = x$vcov.thres3[[i]]$cov_GYI,
                    radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_GYI <- rbind(dt_ell_thres_GYI, data.frame(uu, pts = as.factor(i)))
    }
    dt_ell_thres_GYI$method <- factor(rep(c("Generalized Youden Index Approach"), 52*n_x),
                                      levels = c("Generalized Youden Index Approach"))
    dt_ell_thres$GYI <- dt_ell_thres_GYI
  }
  if("CtP" %in% x$method){
    dt_ell_thres_CtP <- data.frame()
    for(i in 1:n_x){
      uu <- ellipse(center = as.numeric(dt_thres_list$`Closest to Perfection Approach`[i, 1:2]),
                    shape = x$vcov.thres3[[i]]$cov_CtP,
                    radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_CtP <- rbind(dt_ell_thres_CtP, data.frame(uu, pts = as.factor(i)))
    }
    dt_ell_thres_CtP$method <- factor(rep(c("Closest to Perfection Approach"), 52*n_x),
                                      levels = c("Closest to Perfection Approach"))
    dt_ell_thres$CtP <- dt_ell_thres_CtP
  }
  if("MV" %in% x$method){
    dt_ell_thres_MV <- data.frame()
    for(i in 1:n_x){
      uu <- ellipse(center = as.numeric(dt_thres_list$`Max Volume Approach`[i, 1:2]),
                    shape = x$vcov.thres3[[i]]$cov_MV,
                    radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_MV <- rbind(dt_ell_thres_MV, data.frame(uu, pts = as.factor(i)))
    }
    dt_ell_thres_MV$method <- factor(rep(c("Max Volume Approach"), 52*n_x), levels = c("Max Volume Approach"))
    dt_ell_thres$MV <- dt_ell_thres_MV
  }
  dt_ell_thres <- do.call(rbind, dt_ell_thres)
  pp <- ggplot(dt_thres, aes(x = x, y = y, colour = pts)) +
    facet_grid( ~ method) +
    geom_point(size = size.point) +
    geom_path(data = dt_ell_thres, aes(x = x, y = y, group = pts), size = size.path) +
    xlab("Optimal threshold 1") + ylab("Optimal threshold 2") +
    theme_bw() +
    theme(legend.position = "bottom", strip.text.x = element_text(size = 7))
  if(is.null(colors)){
    colors <- topo.colors(n_x)
  } else{
    if(length(colors) != n_x) stop(paste("Number of colors need to be equal:", n_x))
  }
  pp <- pp + scale_color_manual(name = names.labels, breaks = as.character(1:n_x),
                                values = colors, labels = labels)
  if(!missing(xlims) & !missing(ylims)){
    pp <- pp + xlim(xlims) + ylim(ylims)
  }
  if(!is.null(file.name)){
    ggexport(pp, filename = file.name, ...)
  }
  pp
}
#
