####==========================================================================####
## This file consists of functions for estimating optimal threshold             ##
## Date: 28/05/2021																															##
####==========================================================================####

optThres2_core <- function(method, para, x.val, n_p, n_coef, boxcox, start, method.optim,
                           maxit, lower, upper){
  par_model <- para[1:(n_coef + 3)]
  if(boxcox) par_boxcox <- para[length(para)]
  out <- list()
  if("YI" %in% method){
    out$YI <- YI_est_fun(par_model = par_model, x.val = x.val, n_p = n_p, boxcox = boxcox,
                         par_boxcox = par_boxcox, crit.var = 0.01)
  }
  if("NC" %in% method){
    out$NC <- NC_est_fun(par_model = par_model, x.val = x.val, n_p = n_p, boxcox = boxcox,
                         par_boxcox = par_boxcox, start = start, method.optim = method.optim,
                         maxit = maxit, lower = lower, upper = upper)
  }
  if("MA" %in% method){
    out$MA <- MA_est_fun(par_model = par_model, x.val = x.val, n_p = n_p, boxcox = boxcox,
                         par_boxcox = par_boxcox, start = start, method.optim = method.optim,
                         maxit = maxit, lower = lower, upper = upper)
  }
  return(out)
}

optThres2_se <- function(method, thres_est, out_lme2, x.val, n_p, n_coef, bootstrap, nR, type.boot, data,
                         parallel, ncpus, start, method.optim, maxit, lower, upper){
  out <- list()
  ## asymptotic variance under normal distribution
  if(!out_lme2$boxcox & !bootstrap){
    if("YI" %in% method){
      out$var_YI <- YI_normal_var(opt_ctps = as.numeric(thres_est[thres_est$Method == "Youden Index Approach", 1]),
                                  par_model = out_lme2$est_para,
                                  vcov_par_model = out_lme2$vcov_sand, x.val = x.val, n_p = n_p)$var_cpts
    }
    if("NC" %in% method){
      out$var_NC <- NC_normal_var(opt_ctps = as.numeric(thres_est[thres_est$Method == "Northwest Corner Approach", 1]),
                                  par_model = out_lme2$est_para,
                                  vcov_par_model = out_lme2$vcov_sand, x.val = x.val, n_p = n_p)$var_cpts
    }
    if("MA" %in% method){
      out$var_MA <- MA_normal_var(opt_ctps = as.numeric(thres_est[thres_est$Method == "Max Area Approach", 1]),
                                  par_model = out_lme2$est_para,
                                  vcov_par_model = out_lme2$vcov_sand, x.val = x.val, n_p = n_p)$var_cpts
    }
  } else{## bootstrap
    out_bts <- boot_lme2(B = nR, name.test = out_lme2$name.test, name.class = out_lme2$name.class,
                         name.covars = out_lme2$name.covars, name.clust = out_lme2$name.clust,
                         data = data, x.val = x.val, type = type.boot, boxcox = out_lme2$boxcox,
                         parallel = parallel, ncpus = ncpus)
    if("YI" %in% method){
      thres_YI_bts <- numeric(nR)
      for(k in 1:nR){
        thres_YI_bts[k] <- optThres2_core(method = "YI", para = out_bts[,k], x.val = x.val, n_p = n_p,
                                          n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                          method.optim = method.optim, maxit = maxit, lower = lower,
                                          upper = upper)$YI
      }
      out$var_YI <- var(thres_YI_bts, na.rm = TRUE)
    }
    if("NC" %in% method){
      thres_NC_bts <- numeric(nR)
      for(k in 1:nR){
        thres_NC_bts[k] <- optThres2_core(method = "NC", para = out_bts[,k], x.val = x.val, n_p = n_p,
                                          n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                          method.optim = method.optim, maxit = maxit, lower = lower,
                                          upper = upper)$NC
      }
      out$var_NC <- var(thres_NC_bts, na.rm = TRUE)
    }
    if("MA" %in% method){
      thres_MA_bts <- numeric(nR)
      for(k in 1:nR){
        thres_MA_bts[k] <- optThres2_core(method = "MA", para = out_bts[,k], x.val = x.val, n_p = n_p,
                                          n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                          method.optim = method.optim, maxit = maxit, lower = lower,
                                          upper = upper)$MA
      }
      out$var_MA <- var(thres_MA_bts, na.rm = TRUE)
    }
  }
  return(out)
}

#'@export
optThres2 <- function(method = c("YI", "NC", "MA"), out_lme2, x.val, apVar = TRUE,
                      data, ci.level = 0.95, control = list()){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  if(out_lme2$n_coef/out_lme2$n_p != 2) stop("There is not a case of two-class setting!")
  if(missing(method)) stop("Please, input the selection method(s)!")
  controlvals <- optThres3control()
  if (!missing(control)) {
    controlvals[names(control)] <- control
  }
  ## checking the method
  methodtemp <- substitute(me, list(me = method))
  okMethod <- c("YI", "NC", "MA")
  if(length(methodtemp) > 3) stop(gettextf("The maximum number of selection methods are 3!"), domain = NA)
  if(any(duplicated(methodtemp))) stop(gettextf("The selection methods need to be unique!"), domain = NA)
  if(!is.character(methodtemp)) methodtemp <- deparse(methodtemp)
  if(any(sapply(methodtemp, function(x) !is.element(x, okMethod)))){
    stop(gettextf("the selection method(s) should be: %s", paste(sQuote(okMethod), collapse = ", ")),
         domain = NA)
  }
  methodtemp <- methodtemp[c(which(methodtemp == "YI"), which(methodtemp == "NC"), which(methodtemp == "MA"))]
  Method <- character(length(methodtemp))
  Method[methodtemp == "YI"] <- "Youden Index Approach"
  Method[methodtemp == "NC"] <- "Northwest Corner Approach"
  Method[methodtemp == "MA"] <- "Max Area Approach"
  Method <- factor(Method,
                   levels = c("Youden Index Approach", "Northwest Corner Approach", "Max Area Approach"))
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
      out <- optThres2_core(method = methodtemp, para = out_lme2$est_para, x.val = x.val[i], n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(threshold = sapply(out, function(x) x[c("threshold")]),
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
      fit$vcov.thres2 <- list()
      for(i in 1:n_x){
        fit$vcov.thres2[[i]] <- optThres2_se(method = methodtemp, thres_est = temp_thres[[i]],
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
      out <- optThres2_core(method = methodtemp, para = out_lme2$est_para, x.val = x.val[i], n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(threshold = sapply(out, function(x) x[c("threshold")]),
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
      fit$vcov.thres2 <- list()
      for(i in 1:n_x){
        fit$vcov.thres2[[i]] <- optThres2_se(method = methodtemp, thres_est = temp_thres[[i]],
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
      out <- optThres2_core(method = methodtemp, para = out_lme2$est_para, x.val = x.val[i,], n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(threshold = sapply(out, function(x) x[c("threshold")]),
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
      fit$vcov.thres2 <- list()
      for(i in 1:n_x){
        fit$vcov.thres2[[i]] <- optThres2_se(method = methodtemp, thres_est = temp_thres[[i]],
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
  fit$thres2 <- do.call(rbind, temp_thres)
  fit$ci.level <- ci.level
  class(fit) <- "optThres2"
  return(fit)
}



#' @export
print.optThres2 <- function(x, digits = max(3, getOption("digits") - 2),
                            dig.tst = max(1, min(5, digits - 1)), ...){
  if(isFALSE(inherits(x, "optThres2"))) stop("The object is not optThres2!")
  cat("\n")
  cat("CALL: ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  n_method <- length(x$method)
  if(x$n_p == 1){
    labels <- "Intercept"
  }
  if(x$n_p == 2) {
    labels <- rep(as.character(x$x.val), each = n_method)
  }
  if(x$n_p == 3) {
    labels <- rep(apply(x$x.val, 1, function(y) paste0("(", paste(y, collapse = ", "), ")")), each = n_method)
  }
  if(!is.null(x$vcov.thres2)){
    thres_se <- sqrt(as.numeric(sapply(x$vcov.thres2, unlist)))
    ci.tab <- cbind(x$thres2$threshold - thres_se*qnorm((1 + x$ci.level)/2),
                    x$thres2$threshold + thres_se*qnorm((1 + x$ci.level)/2))
    ci.tab <- format(round(ci.tab, digits = digits))
    infer_tab <- data.frame(labels, x$thres2$Method, x$thres2$threshold, thres_se,
                            apply(matrix(ci.tab[,1:2], ncol = 2, byrow = FALSE), 1,
                                  function(y) paste0("(", paste(y, collapse = ", "), ")")))
    infer_tab[,3:4] <- signif(infer_tab[,3:4], digits = digits)
    colnames(infer_tab) <- c("Covariates Values", "Method", "Est.", "Std.Error",
                             paste0(x$ci.level*100, "% CI"))
    cat("Covariate-specific optimal thresholds: \n")
    print(infer_tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, ...)
  }
  else{
    infer_tab <- data.frame(labels, x$thres2$Method, x$thres2$threshold)
    infer_tab[,3] <- signif(infer_tab[,3], digits = digits)
    colnames(infer_tab) <- c("Covariates Values", "Method", "Est.")
    cat("Covariate-specific optimal thresholds: \n")
    print(infer_tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, ...)
  }
  cat("\n")
  invisible(x)
}
