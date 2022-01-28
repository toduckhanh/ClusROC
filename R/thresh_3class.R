####==========================================================================####
## This file consists of functions for estimating optimal pair of thresholds    ##
## Date: 25/05/2021																															##
####==========================================================================####

#' @import grDevices
#' @import car
#' @import ggplot2
#'
optThres3_core <- function(method, para, z, n_p, n_coef, boxcox, start, method.optim,
                           maxit, lower, upper){
  par_model <- para[1:(n_coef + 4)]
  if(boxcox) par_boxcox <- para[length(para)]
  out <- list()
  if("GYI" %in% method){
    out$GYI <- GYI_est_fun(par_model = par_model, z = z, n_p = n_p, boxcox = boxcox,
                           par_boxcox = par_boxcox, crit.var = 0.01)
  }
  if("CtP" %in% method){
    out$CtP <- CtP_est_fun(par_model = par_model, z = z, n_p = n_p, boxcox = boxcox,
                           par_boxcox = par_boxcox, start = start, method.optim = method.optim,
                           maxit = maxit, lower = lower, upper = upper)
  }
  if("MV" %in% method){
    out$MV <- MV_est_fun(par_model = par_model, z = z, n_p = n_p, boxcox = boxcox,
                         par_boxcox = par_boxcox, start = start, method.optim = method.optim,
                         maxit = maxit, lower = lower, upper = upper)
  }
  return(out)
}

optThres3_se <- function(method, thres_est, out_lme2, z, n_p, n_coef, bootstrap, nR, data,
                         parallel, ncpus, start, method.optim, maxit, lower, upper){
  out <- list()
  ## asymptotic variance under normal distribution
  if(!out_lme2$boxcox & !bootstrap){
    if("GYI" %in% method){
      out$cov_GYI <- GYI_normal_var(
        opt_ctps = as.numeric(thres_est[thres_est$Method == "Generalized Youden Index", 1:2]),
        par_model = out_lme2$est_para,
        vcov_par_model = out_lme2$vcov_sand, z = z, n_p = n_p)$vcov_cpts
    }
    if("CtP" %in% method){
      out$cov_CtP <- CtP_normal_var(
        opt_ctps = as.numeric(thres_est[thres_est$Method == "Closest to Perfection", 1:2]),
        par_model = out_lme2$est_para,
        vcov_par_model = out_lme2$vcov_sand, z = z, n_p = n_p)$vcov_cpts
    }
    if("MV" %in% method){
      out$cov_MV <- MV_normal_var(
        opt_ctps = as.numeric(thres_est[thres_est$Method == "Max Volume", 1:2]),
        par_model = out_lme2$est_para,
        vcov_par_model = out_lme2$vcov_sand, z = z, n_p = n_p)$vcov_cpts
    }
  } else{ ## cluster bootstrap
    out_bts <- boot_lme2(B = nR, name.test = out_lme2$name.test, name.class = out_lme2$name.class,
                         name.covars = out_lme2$name.covars, name.clust = out_lme2$name.clust,
                         data = data, z = z, type = "cluster", boxcox = out_lme2$boxcox,
                         parallel = parallel, ncpus = ncpus)
    if("GYI" %in% method){
      thres_GYI_bts <- matrix(0, nrow = 2, ncol = nR)
      for(k in 1:nR){
        thres_GYI_bts[,k] <- optThres3_core(method = "GYI", para = out_bts[,k], z = z, n_p = n_p,
                                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                            method.optim = method.optim, maxit = maxit, lower = lower,
                                            upper = upper)$GYI
      }
      out$cov_GYI <- cov(na.omit(t(thres_GYI_bts)))
    }
    if("CtP" %in% method){
      thres_CtP_bts <- matrix(0, nrow = 2, ncol = nR)
      for(k in 1:nR){
        thres_CtP_bts[,k] <- optThres3_core(method = "CtP", para = out_bts[,k], z = z, n_p = n_p,
                                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                            method.optim = method.optim, maxit = maxit, lower = lower,
                                            upper = upper)$CtP
      }
      out$cov_CtP <- cov(na.omit(t(thres_CtP_bts)))
    }
    if("MV" %in% method){
      thres_MV_bts <- matrix(0, nrow = 2, ncol = nR)
      for(k in 1:nR){
        thres_MV_bts[,k] <- optThres3_core(method = "MV", para = out_bts[,k], z = z, n_p = n_p,
                                           n_coef = n_coef, boxcox = out_lme2$boxcox, start = start,
                                           method.optim = method.optim, maxit = maxit, lower = lower,
                                           upper = upper)$MV
      }
      out$cov_MV <- cov(na.omit(t(thres_MV_bts)))
    }
  }
  return(out)
}

optThres3control <- function(method.optim = c("L-BFGS-B", "BFGS", "Nelder-Mead"),
                             start = NULL, maxit = 200, lower = -Inf, upper = Inf, nR = 250,
                             parallel = FALSE, ncpus = NULL){
  if(parallel & is.null(ncpus)) ncpus <- 2
  list(method.optim = match.arg(method.optim), start = start, maxit = maxit, lower = lower, upper = upper,
       nR = nR, parallel = parallel, ncpus = ncpus)
}


### ---- Estimate the optimal pair of thresholds under normal assumption ----
#' @title Estimate the covariate-specific optimal pair of thresholds for clustered data.
#'
#' @description \code{optThres3} estimates covariate-specific optimal pair of thresholds of a continuous diagnostic test in a clustered design when subjects can be diagnosed in three ordinal groups.
#'
#' @param method  The method to be used. See 'Details'.
#' @param out_lme2  an object of class "lme2", a result of a call to \code{\link{lme2}}.
#' @param x.val  specific value(s) of covariate(s) where the ROC surface are computed. In case non-covariate, no value is needed to specify. In case of one covariate, \code{x.val} should be a number. In case of \eqn{p} covariates (\eqn{p > 1}), \code{x.val} should be a vector containing \eqn{p} values of the covariates.
#' @param apVar  a logical value. If set to TRUE, the variance-covariance matrix of covariate-specific optimal thresholds is computed.
#' @param data  a data frame containing the variables in the model, which used for performing bootstrap procedure to estimate the variance-covariance matrix of covariate-specific optimal thresholds.
#' @param control  a list of control parameters. See 'Details'
#'
#' @details
#' This function implements three estimation methods discussed in To et al. (2021) for estimating covariate-specific optimal pair of thresholds in a clustered design when subjects diagnosed in three ordinal groups. The estimators are based on the results of fitting the linear mixed-effect model on the diagnostic tests, which is done by using \code{\link{lme2}} with REML approach.
#'
#' Before applying the estimation, a quick check for the monotone ordering assumption will be performed. That is, for given values of covariates, three predicted means of three diagnostic groups will be compared. If the assumption does not meet, the covariate-specific optimal pair of thresholds at the values of covariates will be not estimated.
#'
#' Method \code{"GYI"} is Generalized Youden Index. This maximizes the sum of three covariate-specific True Class Fractions - TCFs (or total of correct classification rates) to obtain the optimal pair of thresholds. Method \code{"CtP"} is Closest to Pefection. By using this method, the optimal pair of thresholds is obtained by minimizing the distance, in the unit cube, between the point of three covariate-specific TCFs and the top corner (1, 1, 1). Method \code{"MV"} is Max Volume, which searches for thresholds that maximize the volume of a box under the covariate-specific ROC surface, and the volume is defined as the product of the three covariate-specific TCFs. The user can select more than one method.
#'
#' The asymptotic variance-covariance matrix of the estimated covariate-specific optimal thresholds is approximated through the Delta method under the normal assumptions. If the Box-Cox transformation was applied for the linear mixed-effect on the diagnostic test (\code{\link{lme2}}), a nonparametric bootstrap procedure for clustered data will be applied to obtain the asymptotic covariance matrix the estimated covariate-specific optimal thresholds in original scale (see To et al. 2021, for more details).
#'
#' The \code{control} argument is a list that can supply any of the following components:
#' \tabular{ll}{
#' \code{method.optim} \tab \cr
#' \tab Name of optimization method to be used. There are three options: \code{"L-BFGS-B"}, \code{"BFGS"} and \code{"Nelder-Mead"}. \cr
#' \code{start} \tab \cr
#' \tab Initial values for the thresholds to be optimized over. If it is \code{NULL}, a starting point will be automatically obtained. \cr
#' \code{maxit} \tab \cr
#' \tab The maximum number of iterations. Defaults to 200. \cr
#' \code{lower, upper} \tab \cr
#' \tab Bounds on the variables for the "L-BFGS-B" method. Defaults are \code{-Inf} and \code{Inf}. \cr
#' \code{nR} \tab \cr
#' \tab The number of bootstrap replicates for estimating the covariance matrix under Box-Cox transformation. Defaults to 250. \cr
#' \code{parallel} \tab \cr
#' \tab A logical value. If TRUE, a parallel computing is employed in the bootstrap resampling process. \cr
#' \code{ncpus} \tab \cr
#' \tab Number of processes to be used in parallel computing. Default is 2.
#' }
#'
#' @return \code{optThres3} returns an object of class inheriting from "optThres3" class, which is a list containing at least the following components:
#'
#' \item{call}{the matched call.}
#' \item{method}{the methods were used to obtain the optima threholds.}
#' \item{thres3}{a vector or matrix containing the estimated optimal thresholds.}
#' \item{vcov.thres3}{a matrix or list of matrices containing the estimated variance-covariance matrices.}
#' \item{tcfs}{a vector or matrix containing the estimated TCFs at the optimal thresholds.}
#' \item{mess_order}{a diagnostic message for monontone ordering of means at given covariates' values.}
#' \item{x.val}{value(s) of covariate(s).}
#' \item{n_p}{total numbers of the regressors in the model.}
#'
#' Generic functions such as \code{print} and \code{plot} have methods to show the results.
#'
#' @references
#' To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2021)
#' ``ROC estimation and threshold selection criteria in three-class classification problems for clustered data''. \emph{Submitted}.
#'
#' @examples
#'
#' data(data_3class)
#' ## One covariate
#' out1 <- lme2(fixed.formula = Y ~ X1, name.class = "D", name.clust = "id_Clus",
#'              data = data_3class)
#'
#' ### Estimate covariate-specific optimal thresholds at multiple values of one covariate,
#' ### with 3 methods
#' out_thres_1 <- optThres3(method = c("GYI", "MV", "CtP"), out_lme2 = out1, x.val = 1, apVar = TRUE)
#' print(out_thres_1)
#' plot(out_thres_1)
#'
#'## Two covariates
#' out2 <- lme2(fixed.formula = Y ~ X1 + X2, name.class = "D", name.clust = "id_Clus",
#'              data = data_3class)
#'
#' ### Estimate covariate-specific optimal thresholds at one point, with 3 methods
#' out_thres_2 <- optThres3(method = c("GYI", "MV", "CtP"), out_lme2 = out2, x.val = c(1, 0),
#'                          apVar = TRUE)
#' print(out_thres_2)
#' plot(out_thres_2)
#'
#' ### Estimate covariate-specific optimal thresholds at three points, with 3 methods
#' out_thres_3 <- optThres3(method = c("GYI", "MV", "CtP"), out_lme2 = out2,
#'                          x.val = rbind(c(-0.5, 0), c(0.5, 0), c(0.5, 1)), apVar = TRUE)
#' print(out_thres_3)
#' plot(out_thres_3, colors = c("forestgreen", "blue"))
#'
#'@export
optThres3 <- function(method = c("GYI", "CtP", "MV"), out_lme2, x.val, apVar = TRUE,
                      data, control = list()){
  ## Check all conditions
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  n_p <- out_lme2$n_p
  n_coef <- out_lme2$n_coef
  if(n_coef/n_p != 3) stop("There is not a case of three-class setting!")
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
  Method[methodtemp == "GYI"] <- "Generalized Youden Index"
  Method[methodtemp == "CtP"] <- "Closest to Perfection"
  Method[methodtemp == "MV"] <- "Max Volume"
  Method <- factor(Method, levels = c("Generalized Youden Index", "Closest to Perfection", "Max Volume"))
  ##
  if(n_p == 1){
    if(!missing(x.val)) {
      if(!is.null(x.val)) warning("Sepecified value(s) of covariate(s) are not used!", call. = FALSE)
    }
    x.val <- NULL
  }
  if(n_p == 2){
    if(missing(x.val)) stop("Please input specific value(s) of covariate.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariate.")
    if(!inherits(x.val, "numeric")) stop("For the case of 1 covariate, please input a number or a vector.")
    if(any(is.na(x.val))) stop("NA value(s) are not allowed!")
  }
  if(n_p > 2){
    if(missing(x.val)) stop("Please input specific value(s) of covariates.")
    if(is.null(x.val)) stop("Please input specific value(s) of covariates.")
    n_vb <- attr(out_lme2$terms, "n_vb")
    if(inherits(x.val, "numeric")){
      if(length(x.val) != n_vb) stop(paste("For case of", n_vb, "covariates, please input a vector of", n_vb, "values of covariates."))
    }
    if(inherits(x.val, "matrix")) {
      if(ncol(x.val) != n_vb) stop(paste("For case of m points of", n_vb, "covariates, please input a matrix with", n_vb, "columns and m rows containing values of covariates."))
    }
    if(any(is.na(x.val))) stop("NA value(s) not allowed!")
    x.val <- matrix(x.val, ncol = n_vb, byrow = FALSE)
  }
  ##
  if(apVar){
    if(!out_lme2$boxcox){ ## asymptotic variance under normal distribution
      bootstrap <- FALSE
      if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing!")
      if(any(is.na(out_lme2$vcov_sand))) stop("There are NA values in the estimated covariance matrix of parameters. Unable to estimate standard error.")
    } else {
      bootstrap <- TRUE
      if(missing(data)) stop("The original data is required to process bootstrap procedure!")
    }
  }
  ##
  call <- match.call()
  fit <- list()
  fit$call <- call
  fit$method <- methodtemp
  ## Check the ordering of means: mu_1 < mu_2 < mu_3
  par_model <- out_lme2$est_para
  Z <- make_data(out_lme2, x.val, n_p)
  res_check <- check_mu_order(Z, par_model, n_p)
  if(all(res_check$status == 0))
    stop("The assumption of montone ordering DOES NOT hold for all the value(s) of the covariate(s)")
  if(any(res_check$status == 0)){
    mess_order <- paste("The assumption of montone ordering DOES NOT hold for some points. The points number:",
                        paste(which(res_check$status == 0), collapse = ", "), "are deleted from analysis!")
    fit$mess_order <- mess_order
    message(mess_order)
  }
  Z <- res_check$Z_new
  ##
  if(n_p == 1){# no covariate
    n_x <- 1
    temp_thres <- temp_tcfs <- list()
    for(i in 1:n_x){
      out <- optThres3_core(method = methodtemp, para = par_model, z = Z, n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) x[c("threshold_1", "threshold_2")])),
                                    Method = Method, row.names = NULL)
      temp_tcfs[[i]] <- t(mapply(function(x, y){
        TCF_normal(par = par_model, z = Z, thresholds = c(x, y), n_p = n_p, boxcox = out_lme2$boxcox)
      }, x = temp_thres[[i]]$threshold_1, y = temp_thres[[i]]$threshold_2))
    }
    if(apVar){
      fit$vcov.thres3 <- list()
      for(i in 1:n_x){
        fit$vcov.thres3[[i]] <- optThres3_se(method = methodtemp, thres_est = temp_thres[[i]],
                                             out_lme2 = out_lme2, z = Z, n_p = n_p,
                                             n_coef = n_coef, bootstrap = bootstrap, nR = controlvals$nR,
                                             data = data, parallel = controlvals$parallel,
                                             ncpus = controlvals$ncpus, start = controlvals$start,
                                             method.optim = controlvals$method.optim,
                                             maxit = controlvals$maxit, lower = controlvals$lower,
                                             upper = controlvals$upper)
      }
    }
  }
  if(n_p == 2){
    x.val <- x.val[res_check$status != 0]
    n_x <- length(x.val)
    temp_thres <- temp_tcfs <- list()
    for(i in 1:n_x){
      out <- optThres3_core(method = methodtemp, para = par_model, z = Z[[i]], n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) x[c("threshold_1", "threshold_2")])),
                                    Method = Method, x = x.val[i], row.names = NULL)
      temp_tcfs[[i]] <- t(mapply(function(x, y){
        TCF_normal(par = par_model, z = Z[[i]], thresholds = c(x, y), n_p = n_p, boxcox = out_lme2$boxcox)
      }, x = temp_thres[[i]]$threshold_1, y = temp_thres[[i]]$threshold_2))
    }
    if(apVar){
      fit$vcov.thres3 <- list()
      for(i in 1:n_x){
        fit$vcov.thres3[[i]] <- optThres3_se(method = methodtemp, thres_est = temp_thres[[i]],
                                             out_lme2 = out_lme2, z = Z[[i]], n_p = n_p,
                                             n_coef = n_coef, bootstrap = bootstrap, nR = controlvals$nR,
                                             data = data, parallel = controlvals$parallel,
                                             ncpus = controlvals$ncpus, start = controlvals$start,
                                             method.optim = controlvals$method.optim,
                                             maxit = controlvals$maxit, lower = controlvals$lower,
                                             upper = controlvals$upper)
      }
    }
  }
  if(n_p > 2){ # multiple covariates
    x.val <- matrix(x.val[res_check$status != 0,], ncol = n_vb, byrow = FALSE)
    n_x <- nrow(x.val)
    temp_thres <- temp_tcfs <- list()
    for(i in 1:n_x){
      out <- optThres3_core(method = methodtemp, para = par_model, z = Z[[i]], n_p = n_p,
                            n_coef = n_coef, boxcox = out_lme2$boxcox, start = controlvals$start,
                            method.optim = controlvals$method.optim, maxit = controlvals$maxit,
                            lower = controlvals$lower, upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) x[c("threshold_1", "threshold_2")])),
                                    Method = Method, x = t(x.val[i,]), row.names = NULL)
      temp_tcfs[[i]] <- t(mapply(function(x, y){
        TCF_normal(par = par_model, z = Z[[i]], thresholds = c(x, y), n_p = n_p, boxcox = out_lme2$boxcox)
      }, x = temp_thres[[i]]$threshold_1, y = temp_thres[[i]]$threshold_2))
    }
    if(apVar){
      fit$vcov.thres3 <- list()
      for(i in 1:n_x){
        fit$vcov.thres3[[i]] <- optThres3_se(method = methodtemp, thres_est = temp_thres[[i]],
                                             out_lme2 = out_lme2, z = Z[[i]], n_p = n_p,
                                             n_coef = n_coef, bootstrap = bootstrap, nR = controlvals$nR,
                                             data = data, parallel = controlvals$parallel,
                                             ncpus = controlvals$ncpus, start = controlvals$start,
                                             method.optim = controlvals$method.optim,
                                             maxit = controlvals$maxit, lower = controlvals$lower,
                                             upper = controlvals$upper)
      }
    }
  }
  fit$x.val <- x.val
  fit$n_p <- n_p
  fit$thres3 <- do.call(rbind, temp_thres)
  fit$tcfs <- do.call(rbind, temp_tcfs)
  class(fit) <- "optThres3"
  return(fit)
}


## ---- The function plot.optThres3 ----
#' @title Plot confidence regions for covariate-specific optimal thresholds.
#'
#' @description Plot covariate-specific optimal thresholds and the corresponding confidence regions.
#'
#' @method plot optThres3
#' @param x an object of class "optThres3", a result of a call to \code{\link{optThres3}}.
#' @param ci.level  a confidence level to be used for constructing the confidence interval; default is 0.95.
#' @param colors  specification color(s) for confidence regions. If specified, the number of colors needs to be equal the number of covariates
#' @param xlims,ylims numeric vectors of length 2, giving the x and y coordinates ranges (the possible range of the first and the second threshold).
#' @param size.point,size.path  numeric value of length 1, giving the amount by which plotting points and lines.
#' @param names.labels a optional character vector giving the label name of covariates.
#' @param file.name File name to create on disk.
#' @param ... further arguments passed to \code{\link{ggexport}} method, for example, \code{width}, \code{height}.
#'
#' @details \code{plot.optThres3} provides plots for displaying the confidence regions and point estimates of covariate-specific optimal thresholds. The plots are based on \code{ggplot()}.
#'
#' @seealso \code{\link{optThres3}}
#'
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
  if(x$n_p > 2) {
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
      uu <- ellipse(center = as.numeric(dt_thres_list$`Generalized Youden Index`[i, 1:2]),
                    shape = x$vcov.thres3[[i]]$cov_GYI,
                    radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_GYI <- rbind(dt_ell_thres_GYI, data.frame(uu, pts = as.factor(i)))
    }
    dt_ell_thres_GYI$method <- factor(rep(c("Generalized Youden Index"), 52*n_x),
                                      levels = c("Generalized Youden Index"))
    dt_ell_thres$GYI <- dt_ell_thres_GYI
  }
  if("CtP" %in% x$method){
    dt_ell_thres_CtP <- data.frame()
    for(i in 1:n_x){
      uu <- ellipse(center = as.numeric(dt_thres_list$`Closest to Perfection`[i, 1:2]),
                    shape = x$vcov.thres3[[i]]$cov_CtP,
                    radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_CtP <- rbind(dt_ell_thres_CtP, data.frame(uu, pts = as.factor(i)))
    }
    dt_ell_thres_CtP$method <- factor(rep(c("Closest to Perfection"), 52*n_x),
                                      levels = c("Closest to Perfection"))
    dt_ell_thres$CtP <- dt_ell_thres_CtP
  }
  if("MV" %in% x$method){
    dt_ell_thres_MV <- data.frame()
    for(i in 1:n_x){
      uu <- ellipse(center = as.numeric(dt_thres_list$`Max Volume`[i, 1:2]),
                    shape = x$vcov.thres3[[i]]$cov_MV,
                    radius = sqrt(qchisq(ci.level, 2)), draw = FALSE)
      dt_ell_thres_MV <- rbind(dt_ell_thres_MV, data.frame(uu, pts = as.factor(i)))
    }
    dt_ell_thres_MV$method <- factor(rep(c("Max Volume"), 52*n_x), levels = c("Max Volume"))
    dt_ell_thres$MV <- dt_ell_thres_MV
  }
  dt_ell_thres <- do.call(rbind, dt_ell_thres)
  pp <- ggplot(dt_thres, aes_string(x = "x", y = "y", colour = "pts")) +
    facet_grid( ~ method) +
    geom_point(size = size.point) +
    geom_path(data = dt_ell_thres, aes_string(x = "x", y = "y", group = "pts"), size = size.path) +
    xlab("Optimal threshold 1") + ylab("Optimal threshold 2") +
    theme_bw() +
    theme(legend.position = "bottom", strip.text.x = element_text(size = 9))
  if(is.null(colors)){
    colors <- topo.colors(n_x)
  } else{
    if(length(colors) != n_x) stop(paste("Number of colors needs to be equal:", n_x))
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

## ---- The function print.optThres3 ----
#' @title Print summary results of \code{optThres3}
#'
#' @description \code{print.optThres3} prints the results for the output of function \code{\link{optThres3}}.
#'
#' @method print optThres3
#' @param x an object of class "optThres3", a result of a call to \code{\link{optThres3}}.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param call logical. If \code{TRUE}, the matched call will be printed.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.optThres3} shows a nice format of the summary table for covariate-specific optimal thresholds estimates.
#'
#' @seealso \code{\link{optThres3}}
#'
#' @export
print.optThres3 <- function(x, digits = 3, call = TRUE, ...){
  if(isFALSE(inherits(x, "optThres3"))) stop("The object is not optThres3!")
  cat("\n")
  if(call){
    cat("CALL: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  }
  if(!is.null(x$mess_order)){
    cat("NOTE: ", x$mess_order, "\n \n", sep = "")
  }
  n_method <- length(x$method)
  if(x$n_p == 1){
    labels <- "Intercept"
  }
  if(x$n_p == 2) {
    labels <- rep(format(x$x.val), each = n_method)
  }
  if(x$n_p > 2) {
    labels <- rep(apply(x$x.val, 1, function(y) paste0("(", paste(y, collapse = ", "), ")")), each = n_method)
  }
  infer_tab <- data.frame(labels, x$thres3$Method, x$thres3$threshold_1, x$thres3$threshold_2, x$tcfs[,1],
                          x$tcfs[,2], x$tcfs[,3])
  infer_tab[,3:7] <- signif(infer_tab[,3:7], digits = digits)
  colnames(infer_tab) <- c("Covariate(s) Values", "Method", "Threshold 1", "Threshold 2", "TCF 1", "TCF 2",
                           "TCF 3")
  cat("Covariate-specific optimal pair of thresholds: \n")
  print(infer_tab, quote = FALSE, right = TRUE, na.print = "--", row.names = FALSE, ...)
  cat("\n")
  invisible(x)
}


