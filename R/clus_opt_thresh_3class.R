####========================================================================####
## This file consists of functions for estimating optimal pair of thresholds  ##
## Date: 10/10/2022																														##
####========================================================================####

#' @import grDevices
#' @import car
#' @import ggplot2
#'
clus_opt_thres3_core <- function(method, para, z, n_p, n_coef, boxcox, start,
                                 method_optim, maxit, lower, upper) {
  par_model <- para[1:(n_coef + 4)]
  if (boxcox) {
    par_boxcox <- para[length(para)]
  }
  out <- list()
  if ("GYI" %in% method) {
    out$gyi <- gyi_est_fun(par_model = par_model, z = z, n_p = n_p,
                           boxcox = boxcox, par_boxcox = par_boxcox,
                           crit_var = 0.01)
  }
  if ("CtP" %in% method) {
    out$ctp <- ctp_est_fun(par_model = par_model, z = z, n_p = n_p,
                           boxcox = boxcox, par_boxcox = par_boxcox,
                           start = start, method_optim = method_optim,
                           maxit = maxit, lower = lower, upper = upper)
  }
  if ("MV" %in% method) {
    out$mv <- mv_est_fun(par_model = par_model, z = z, n_p = n_p,
                         boxcox = boxcox, par_boxcox = par_boxcox,
                         start = start, method_optim = method_optim,
                         maxit = maxit, lower = lower, upper = upper)
  }
  return(out)
}

clus_opt_thres3_se <- function(method, thres_est, out_clus_lme, z, n_p, n_coef,
                               bootstrap, n_boot, data, parallel, ncpus, start,
                               method_optim, maxit, lower, upper) {
  out <- list()
  ## asymptotic variance under normal distribution
  if (!out_clus_lme$boxcox && !bootstrap) {
    if ("GYI" %in% method) {
      out$cov_gyi <- gyi_normal_var(
        opt_ctps = as.numeric(
          thres_est[thres_est$Method == "Generalized Youden Index", 1:2]),
        par_model = out_clus_lme$est_para,
        vcov_par_model = out_clus_lme$vcov_sand, z = z, n_p = n_p)$vcov_cpts
    }
    if ("CtP" %in% method) {
      out$cov_ctp <- ctp_normal_var(
        opt_ctps = as.numeric(
          thres_est[thres_est$Method == "Closest to Perfection", 1:2]),
        par_model = out_clus_lme$est_para,
        vcov_par_model = out_clus_lme$vcov_sand, z = z, n_p = n_p)$vcov_cpts
    }
    if ("MV" %in% method) {
      out$cov_mv <- mv_normal_var(
        opt_ctps = as.numeric(thres_est[thres_est$Method == "Max Volume", 1:2]),
        par_model = out_clus_lme$est_para,
        vcov_par_model = out_clus_lme$vcov_sand, z = z, n_p = n_p)$vcov_cpts
    }
  } else { ## cluster bootstrap
    out_bts <- boot_clus_lme(out_clus_lme = out_clus_lme, data = data, z = z,
                             n_boot = n_boot, type = "cluster",
                             boxcox = out_clus_lme$boxcox, parallel = parallel,
                             ncpus = ncpus)
    if ("GYI" %in% method) {
      thres_gyi_bts <- matrix(0, nrow = 2, ncol = n_boot)
      for (k in 1:n_boot) {
        thres_gyi_bts[, k] <- clus_opt_thres3_core(
          method = "GYI", para = out_bts[, k], z = z, n_p = n_p,
          n_coef = n_coef, boxcox = out_clus_lme$boxcox, start = start,
          method_optim = method_optim, maxit = maxit, lower = lower,
          upper = upper)$gyi
      }
      out$cov_gyi <- cov(na.omit(t(thres_gyi_bts)))
    }
    if ("CtP" %in% method) {
      thres_ctp_bts <- matrix(0, nrow = 2, ncol = n_boot)
      for (k in 1:n_boot) {
        thres_ctp_bts[, k] <- clus_opt_thres3_core(
          method = "CtP", para = out_bts[, k], z = z, n_p = n_p,
          n_coef = n_coef, boxcox = out_clus_lme$boxcox, start = start,
          method_optim = method_optim, maxit = maxit, lower = lower,
          upper = upper)$ctp
      }
      out$cov_ctp <- cov(na.omit(t(thres_ctp_bts)))
    }
    if ("MV" %in% method) {
      thres_mv_bts <- matrix(0, nrow = 2, ncol = n_boot)
      for (k in 1:n_boot){
        thres_mv_bts[, k] <- clus_opt_thres3_core(
          method = "MV", para = out_bts[, k], z = z, n_p = n_p,
          n_coef = n_coef, boxcox = out_clus_lme$boxcox, start = start,
          method_optim = method_optim, maxit = maxit, lower = lower,
          upper = upper)$mv
      }
      out$cov_mv <- cov(na.omit(t(thres_mv_bts)))
    }
  }
  return(out)
}

clus_opt_thres3_control <- function(method_optim = c("L-BFGS-B", "BFGS",
                                                     "Nelder-Mead"),
                                    start = NULL, maxit = 200, lower = -Inf,
                                    upper = Inf, n_boot = 250, parallel = FALSE,
                                    ncpus = NULL) {
  if (parallel && is.null(ncpus)) {
    ncpus <- 2
  }
  list(method_optim = match.arg(method_optim), start = start, maxit = maxit,
       lower = lower, upper = upper, n_boot = n_boot, parallel = parallel,
       ncpus = ncpus)
}


### ---- Estimate the optimal pair of thresholds under normal assumption ----
#' @title Estimation of the covariate-specific optimal pair of thresholds for clustered data.
#'
#' @description \code{clus_opt_thres3} estimates covariate-specific optimal pair of thresholds of a continuous diagnostic test in a clustered design, with three classes of diseases.
#'
#' @param method  the method to be used. See 'Details'.
#' @param out_clus_lme  an object of class "clus_lme", i.e., a result of \code{\link{clus_lme}} call.
#' @param newdata  a data frame (containing specific value(s) of covariate(s)) in which to look for variables with which to estimate covariate-specific optimal pair of thresholds. In absence of covariate, no values have to be specified.
#' @param ap_var  logical value. If set to \code{TRUE}, the variance-covariance matrix of (estimated) covariate-specific optimal thresholds is estimated.
#' @param data  a data frame containing the variables to be used when performing a bootstrap procedure to estimate the variance-covariance matrix, in case of Box-Cox transformation.
#' @param control  a list of control parameters. See 'Details'.
#'
#' @details
#' This function implements estimation methods discussed in To et al. (2022) for covariate-specific optimal pair of thresholds in a clustered design with three ordinal groups. The estimators are based on the results from \code{\link{clus_lme}} function, which fits the linear mixed-effect model by using REML approach.
#'
#' Before performing estimation, a check for the monotone ordering assumption is performed. This means that, for the fixed values of covariates, three predicted mean values for test results in three diagnostic groups are compared. If the assumption is not meet, the covariate-specific optimal pair of thresholds at the values of covariates are not estimated.
#'
#' The estimation procedure uses three criteria. Method \code{"GYI"} is Generalized Youden Index, which maximizes the sum of three covariate-specific True Class Fractions - TCFs. Method \code{"CtP"} is based on Closest to Pefection approach. By using this method, the optimal pair of thresholds is obtained by minimizing the distance, in the unit cube, between a generic point on the covariate-specific ROC surface and the top corner (1, 1, 1). Method \code{"MV"} is based on Maximum Volume approach, which searches for thresholds that maximize the volume of a box under the covariate-specific ROC surface. The user can select more than one method. This function allows to estimate covariate-specific optimal pair of thresholds at multiple points for covariates.
#'
#' The asymptotic variance-covariance matrix of the (estimated) covariate-specific optimal thresholds is estimated by using the Delta method under the normal assumption. If the Box-Cox transformation is applied to the linear mixed-effect model, a nonparametric bootstrap procedure for clustered data will be used to obtain the estimated asymptotic covariance matrix (see To et al. 2022, for more details).
#'
#' The \code{control} argument is a list that can supply any of the following components:
#' \describe{
#'   \item{\code{method_optim}}{Optimization method to be used. There are three options: \code{"L-BFGS-B"}, \code{"BFGS"} and \code{"Nelder-Mead"}. Default is \code{"L-BFGS-B"}.}
#'   \item{\code{start}}{Starting values in the optimization procedure. If it is \code{NULL}, a starting point will be automatically obtained.}
#'   \item{\code{maxit}}{The maximum number of iterations. Default is 200.}
#'   \item{\code{lower, upper}}{Possible bounds on the threshold range, for the optimization based on "L-BFGS-B" method. Defaults are \code{-Inf} and \code{Inf}.}
#'   \item{\code{n_boot}}{Number of bootstrap replicates for estimating the covariance matrix (when Box-Cox transformation is applied). Default is 250.}
#'   \item{\code{parallel}}{A logical value. If set to \code{TRUE}, a parallel computing is employed in the bootstrap resampling process.}
#'   \item{\code{ncpus}}{Number of processes to be used in parallel computing. Default is 2.}
#'}
#'
#' @return \code{clus_opt_thres3} returns an object of "clus_opt_thres3" class, which is a list containing at least the following components:
#'
#' \item{call}{the matched call.}
#' \item{method}{the methods used to obtain the estimated optimal pair of threholds.}
#' \item{thres3}{a vector or matrix containing the estimated optimal thresholds.}
#' \item{thres3_se}{a vector or matrix containing the estimated standard errors.}
#' \item{vcov_thres3}{a matrix or list of matrices containing the estimated variance-covariance matrices.}
#' \item{tcfs}{a vector or matrix containing the estimated TCFs at the optimal thresholds.}
#' \item{mess_order}{a diagnostic message from checking the monontone ordering.}
#' \item{newdata}{value(s) of covariate(s).}
#' \item{n_p}{total number of regressors in the model.}
#'
#' Generic functions such as \code{print} and \code{plot} are also used to show the results.
#'
#' @references
#' To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2022)
#' ``Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data''. \emph{Statistical Methods in Medical Research}, \bold{7}, 31, 1325-1341.
#'
#' @examples
#' data(data_3class)
#' ## One covariate
#' out1 <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
#'                  name_clust = "id_Clus", data = data_3class)
#'
#' ### Estimate covariate-specific optimal thresholds at multiple values of one covariate,
#' ### with 3 methods
#' out_thres_1 <- clus_opt_thres3(method = c("GYI", "MV", "CtP"),
#'                                out_clus_lme = out1,
#'                                newdata = data.frame(X1 = 1), ap_var = TRUE)
#' print(out_thres_1)
#' plot(out_thres_1)
#'
#'## Two covariates
#' out2 <- clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
#'                  name_clust = "id_Clus", data = data_3class)
#'
#' ### Estimate covariate-specific optimal thresholds at one point, with 3 methods
#' out_thres_2 <- clus_opt_thres3(method = c("GYI", "MV", "CtP"),
#'                                out_clus_lme = out2,
#'                                newdata = data.frame(X1 = 1, X2 = 0),
#'                                ap_var = TRUE)
#' print(out_thres_2)
#' plot(out_thres_2)
#'
#' ### Estimate covariate-specific optimal thresholds at three points, with 3 methods
#' out_thres_3 <- clus_opt_thres3(method = c("GYI", "MV", "CtP"),
#'                                out_clus_lme = out2,
#'                                newdata = data.frame(X1 = c(-0.5, 0.5, 0.5),
#'                                                     X2 = c(0, 0, 1)),
#'                                ap_var = TRUE)
#' print(out_thres_3)
#' plot(out_thres_3, colors = c("forestgreen", "blue"))
#'
#'@export
clus_opt_thres3 <- function(method = c("GYI", "CtP", "MV"), out_clus_lme,
                            newdata, ap_var = TRUE, data, control = list()) {
  ## Check all conditions
  if (isFALSE(inherits(out_clus_lme, "clus_lme"))) {
    stop("out_clus_lme was not from clus_lme()!")
  }
  n_p <- out_clus_lme$n_p
  n_coef <- out_clus_lme$n_coef
  if (n_coef / n_p != 3) {
    stop("There is not a case of three-class setting!")
  }
  if (missing(method)) {
    stop("Please, input the selection method(s)!")
  }
  controlvals <- clus_opt_thres3_control()
  if (!missing(control)) {
    controlvals[names(control)] <- control
  }
  ## checking the method
  methodtemp <- substitute(me, list(me = method))
  ok_method <- c("GYI", "CtP", "MV")
  if (length(methodtemp) > 3) {
    stop(gettextf("The maximum number of selection methods are 3!"),
         domain = NA)
  }
  if (any(duplicated(methodtemp))) {
    stop(gettextf("The selection methods need to be unique!"), domain = NA)
  }
  if (!is.character(methodtemp)) {
    methodtemp <- deparse(methodtemp)
  }
  if (any(sapply(methodtemp, function(x) !is.element(x, ok_method)))) {
    stop(gettextf("the selection method(s) should be: %s",
                  paste(sQuote(ok_method), collapse = ", ")),
         domain = NA)
  }
  methodtemp <- methodtemp[c(which(methodtemp == "GYI"),
                             which(methodtemp == "CtP"),
                             which(methodtemp == "MV"))]
  mmethod <- character(length(methodtemp))
  mmethod[methodtemp == "GYI"] <- "Generalized Youden Index"
  mmethod[methodtemp == "CtP"] <- "Closest to Perfection"
  mmethod[methodtemp == "MV"] <- "Max Volume"
  mmethod <- factor(mmethod, levels = c("Generalized Youden Index",
                                        "Closest to Perfection", "Max Volume"))
  ##
  if (n_p == 1) {
    if (!missing(newdata)) {
      if (!is.null(newdata)) {
        warning("Sepecified value(s) of covariate(s) are not used!",
                call. = FALSE)
      }
    }
    newdata <- NULL
  } else {
    if (missing(newdata)) {
      stop("Please input a data frame including specific value(s) of covariate(s).")
    }
    if (is.null(newdata)) {
      stop("Please input a data frame including specific value(s) of covariate(s).")
    }
    if (!inherits(newdata, "data.frame")) {
      stop("Please input a data frame including specific value(s) of covariate(s).")
    }
    if (any(is.na(newdata))) {
      stop("NA value(s) are not allowed!")
    }
  }
  ##
  if (ap_var) {
    if (!out_clus_lme$boxcox) { ## asymptotic variance under normal distribution
      bootstrap <- FALSE
      if (is.null(out_clus_lme$vcov_sand)) {
        stop("The estimated covariance matrix of parameters was missing!")
      }
      if (any(is.na(out_clus_lme$vcov_sand))) {
        stop("There are NA values in the estimated covariance matrix of parameters. Unable to estimate standard error.")
      }
    } else {
      bootstrap <- TRUE
      if (missing(data)) {
        stop("The original data is required to process bootstrap procedure!")
      }
    }
  }
  ##
  call <- match.call()
  fit <- list()
  fit$call <- call
  fit$method <- methodtemp
  ## Check the ordering of means: mu_1 < mu_2 < mu_3
  par_model <- out_clus_lme$est_para
  z <- make_data(out_clus_lme, newdata, n_p)
  res_check <- check_mu_order(z, par_model, n_p)
  if (all(res_check$status == 0)) {
    stop("The assumption of montone ordering DOES NOT hold for all the value(s) of the covariate(s)")
  }
  if (any(res_check$status == 0)) {
    mess_order <- paste("The assumption of montone ordering DOES NOT hold for some points. The points number:",
                        paste(which(res_check$status == 0), collapse = ", "),
                        "are deleted from analysis!")
    fit$mess_order <- mess_order
    message(mess_order)
  }
  z <- res_check$z_new
  ##
  if (n_p == 1) {# without covariate
    fit$newdata <- newdata
    n_x <- 1
    temp_thres <- temp_tcfs <- list()
    for (i in 1:n_x) {
      out <- clus_opt_thres3_core(
        method = methodtemp, para = par_model, z = z[[i]], n_p = n_p,
        n_coef = n_coef, boxcox = out_clus_lme$boxcox,
        start = controlvals$start, method_optim = controlvals$method_optim,
        maxit = controlvals$maxit, lower = controlvals$lower,
        upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) {
        x[c("threshold_1", "threshold_2")]
        })), Method = mmethod, row.names = NULL)
      temp_tcfs[[i]] <- t(mapply(function(x, y) {
        tcf_normal(par = par_model, z = z[[i]], thresholds = c(x, y), n_p = n_p,
                   boxcox = out_clus_lme$boxcox)
      }, x = temp_thres[[i]]$threshold_1, y = temp_thres[[i]]$threshold_2))
    }
    if (ap_var) {
      temp_se_thres <- list()
      fit$vcov_thres3 <- list()
      for (i in 1:n_x) {
        out_var <- clus_opt_thres3_se(
          method = methodtemp, thres_est = temp_thres[[i]],
          out_clus_lme = out_clus_lme, z = z[[i]], n_p = n_p, n_coef = n_coef,
          bootstrap = bootstrap, n_boot = controlvals$n_boot, data = data,
          parallel = controlvals$parallel, ncpus = controlvals$ncpus,
          start = controlvals$start, method_optim = controlvals$method_optim,
          maxit = controlvals$maxit, lower = controlvals$lower,
          upper = controlvals$upper)
        fit$vcov_thres3[[i]] <- out_var
        se_thres3 <- t(sapply(out_var, function(x) sqrt(diag(x))))
        temp_se_thres[[i]] <- data.frame(threshold_1 = se_thres3[, 1],
                                         threshold_2 = se_thres3[, 2],
                                         Method = mmethod, row.names = NULL)
      }
    }
  } else {
    fit$newdata <- as.data.frame(newdata[res_check$status != 0, ])
    names(fit$newdata) <- names(newdata)
    n_x <- nrow(fit$newdata)
    temp_thres <- temp_tcfs <- list()
    for (i in 1:n_x) {
      out <- clus_opt_thres3_core(
        method = methodtemp, para = par_model, z = z[[i]], n_p = n_p,
        n_coef = n_coef, boxcox = out_clus_lme$boxcox,
        start = controlvals$start, method_optim = controlvals$method_optim,
        maxit = controlvals$maxit, lower = controlvals$lower,
        upper = controlvals$upper)
      temp_thres[[i]] <- data.frame(t(sapply(out, function(x) {
        x[c("threshold_1", "threshold_2")]
        })), Method = mmethod, x = fit$newdata[i, ], row.names = NULL)
      temp_tcfs[[i]] <- t(mapply(function(x, y) {
        tcf_normal(par = par_model, z = z[[i]], thresholds = c(x, y), n_p = n_p,
                   boxcox = out_clus_lme$boxcox)
      }, x = temp_thres[[i]]$threshold_1, y = temp_thres[[i]]$threshold_2))
    }
    if (ap_var) {
      temp_se_thres <- list()
      fit$vcov_thres3 <- list()
      for (i in 1:n_x) {
        out_var <- clus_opt_thres3_se(
          method = methodtemp, thres_est = temp_thres[[i]],
          out_clus_lme = out_clus_lme, z = z[[i]], n_p = n_p, n_coef = n_coef,
          bootstrap = bootstrap, n_boot = controlvals$n_boot, data = data,
          parallel = controlvals$parallel, ncpus = controlvals$ncpus,
          start = controlvals$start, method_optim = controlvals$method_optim,
          maxit = controlvals$maxit, lower = controlvals$lower,
          upper = controlvals$upper)
        fit$vcov_thres3[[i]] <- out_var
        se_thres3 <- t(sapply(out_var, function(x) sqrt(diag(x))))
        temp_se_thres[[i]] <- data.frame(threshold_1 = se_thres3[, 1],
                                         threshold_2 = se_thres3[, 2],
                                         Method = mmethod, x = fit$newdata[i, ],
                                         row.names = NULL)
      }
    }
  }
  fit$n_p <- n_p
  fit$thres3 <- do.call(rbind, temp_thres)
  fit$tcfs <- do.call(rbind, temp_tcfs)
  if (ap_var) {
    fit$thres3_se <- do.call(rbind, temp_se_thres)
  }
  class(fit) <- "clus_opt_thres3"
  return(fit)
}


## ---- The function plot.clus_opt_thres3 ----
#' @title Plot of confidence regions for covariate-specific optimal pair of thresholds.
#'
#' @description This function plots confidence regions for covariate-specific optimal pair of thresholds.
#'
#' @method plot clus_opt_thres3
#' @param x  an object of class "clus_opt_thres3", i.e., a result of \code{\link{clus_opt_thres3}}.
#' @param ci_level  confidence level to be used for constructing the confidence regions; default is 0.95.
#' @param colors  a string vector for the name(s) specifying color(s) to be used for drawing confidence regions. If specified, the dimension of the vector needs to be equal the number of considered points (each point corresponds to a set of values for the covariates).
#' @param xlims,ylims numeric vectors of dimension 2, giving the limits for x and y axes in the plot.
#' @param size_point,size_path  numeric values, indicating sizes for point(s) and line(s) in the plot.
#' @param names_labels a optional character vector giving the label name for covariates.
#' @param file_name File name to create on disk.
#' @param ... further arguments used with \code{\link{ggexport}} function, for example, \code{width}, \code{height}.
#'
#' @details \code{plot.clus_opt_thres3} provides plots of confidence regions (and point estimates) of covariate-specific optimal pair of thresholds. The plots are based on \code{ggplot()}.
#'
#' @return \code{plot.clus_opt_thres3} returns plots of confidence regions of covariate-specific optimal pair of thresholds.
#'
#' @seealso \code{\link{clus_opt_thres3}}
#'
#'@export
plot.clus_opt_thres3 <- function(x, ci_level = 0.95, colors = NULL, xlims,
                                 ylims, size_point = 0.5, size_path = 0.5,
                                 names_labels, file_name = NULL, ...) {
  if (isFALSE(inherits(x, "clus_opt_thres3"))) {
    stop("x was not from clus_opt_thres3()!")
  }
  if (x$n_p == 1) {
    n_x <- 1
    labels <- "Intercept"
    if (missing(names_labels)) {
      names_labels <- " "
    }
  }
  if (x$n_p == 2) {
    n_x <- nrow(x$newdata)
    labels <- apply(x$newdata, 1, function(y) paste0(y))
    if (missing(names_labels)) {
      names_labels <- "Value(s) of covariate:"
    }
  }
  if (x$n_p > 2) {
    n_x <- nrow(x$newdata)
    labels <- apply(x$newdata, 1, function(y) {
      paste0("(", paste(y, collapse = ", "), ")")
    })
    if (missing(names_labels)) {
      names_labels <- "Value(s) of covariates:"
    }
  }
  dt_thres <- x$thres3[, 1:3]
  colnames(dt_thres) <- c("x", "y", "method")
  dt_thres$pts <- as.factor(rep(1:n_x, each = length(x$method)))
  pp <- ggplot(dt_thres, aes_string(x = "x", y = "y", colour = "pts")) +
    facet_grid(~ method) +
    geom_point(size = size_point) +
    xlab("Optimal threshold 1") +
    ylab("Optimal threshold 2") +
    theme_bw() +
    theme(legend.position = "bottom", strip.text.x = element_text(size = 9))
  if (!is.null(x$vcov_thres3)) {
    c_alp <- sqrt(qchisq(ci_level, 2))
    dt_thres_list <- split(dt_thres, dt_thres$method)
    dt_ell_thres <- list()
    if ("GYI" %in% x$method) {
      dt_ell_thres_gyi <- data.frame()
      for (i in 1:n_x){
        uu <- ellipse(
          center = as.numeric(dt_thres_list$`Generalized Youden Index`[i, 1:2]),
          shape = x$vcov_thres3[[i]]$cov_gyi, radius = c_alp, draw = FALSE)
        dt_ell_thres_gyi <- rbind(dt_ell_thres_gyi,
                                  data.frame(uu, pts = as.factor(i)))
      }
      dt_ell_thres_gyi$method <- factor(rep(c("Generalized Youden Index"),
                                            52 * n_x),
                                        levels = c("Generalized Youden Index"))
      dt_ell_thres$gyi <- dt_ell_thres_gyi
    }
    if ("CtP" %in% x$method) {
      dt_ell_thres_ctp <- data.frame()
      for (i in 1:n_x) {
        uu <- ellipse(
          center = as.numeric(dt_thres_list$`Closest to Perfection`[i, 1:2]),
          shape = x$vcov_thres3[[i]]$cov_ctp, radius = c_alp, draw = FALSE)
        dt_ell_thres_ctp <- rbind(dt_ell_thres_ctp,
                                  data.frame(uu, pts = as.factor(i)))
      }
      dt_ell_thres_ctp$method <- factor(rep(c("Closest to Perfection"),
                                            52 * n_x),
                                        levels = c("Closest to Perfection"))
      dt_ell_thres$ctp <- dt_ell_thres_ctp
    }
    if ("MV" %in% x$method) {
      dt_ell_thres_mv <- data.frame()
      for (i in 1:n_x) {
        uu <- ellipse(
          center = as.numeric(dt_thres_list$`Max Volume`[i, 1:2]),
          shape = x$vcov_thres3[[i]]$cov_mv,
          radius = c_alp, draw = FALSE)
        dt_ell_thres_mv <- rbind(dt_ell_thres_mv,
                                 data.frame(uu, pts = as.factor(i)))
      }
      dt_ell_thres_mv$method <- factor(rep(c("Max Volume"), 52 * n_x),
                                       levels = c("Max Volume"))
      dt_ell_thres$mv <- dt_ell_thres_mv
    }
    dt_ell_thres <- do.call(rbind, dt_ell_thres)
    pp <- pp +
      geom_path(data = dt_ell_thres,
                aes_string(x = "x", y = "y", group = "pts"), size = size_path)
  }
  if (is.null(colors)) {
    colors <- topo.colors(n_x)
  } else {
    if (length(colors) != n_x) {
      stop(paste("Number of colors needs to be equal:", n_x))
    }
  }
  pp <- pp +
    scale_color_manual(name = names_labels, breaks = as.character(1:n_x),
                       values = colors, labels = labels)
  if (!missing(xlims) && !missing(ylims)) {
    pp <- pp + xlim(xlims) + ylim(ylims)
  }
  if (!is.null(file_name)) {
    ggexport(pp, filename = file_name, ...)
  }
  pp
}
#

## ---- The function print.clus_opt_thres3 ----
#' @title Print summary results from \code{clus_opt_thres3}
#'
#' @description \code{print.clus_opt_thres3} displays the results of the output from \code{\link{clus_opt_thres3}}.
#'
#' @method print clus_opt_thres3
#' @param x an object of class "clus_opt_thres3", a result of \code{\link{clus_opt_thres3}} call.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param call logical. If set to \code{TRUE}, the matched call will be printed.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.clus_opt_thres3} shows a summary table for covariate-specific optimal pair of thresholds estimates.
#'
#' @return \code{print.clus_opt_thres3} returns a summary table for results of covariate-specific optimal pair of thresholds estimation.
#'
#' @seealso \code{\link{clus_opt_thres3}}
#'
#' @export
print.clus_opt_thres3 <- function(x, digits = 3, call = TRUE, ...) {
  if (isFALSE(inherits(x, "clus_opt_thres3"))) {
    stop("The object is not clus_opt_thres3!")
  }
  cat("\n")
  if (call) {
    cat("CALL: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n \n", sep = "")
  }
  if (!is.null(x$mess_order)) {
    cat("NOTE: ", x$mess_order, "\n \n", sep = "")
  }
  n_method <- length(x$method)
  if (x$n_p == 1) {
    labels <- "Intercept"
  }
  if (x$n_p == 2) {
    labels <- rep(apply(x$newdata, 1, function(y) paste0(y)), each = n_method)
  }
  if (x$n_p > 2) {
    labels <- rep(apply(x$newdata, 1, function(y) {
      paste0("(", paste(y, collapse = ", "), ")")
    }), each = n_method)
  }
  infer_tab <- data.frame(labels, x$thres3$Method, x$thres3$threshold_1,
                          x$thres3$threshold_2, x$tcfs[, 1], x$tcfs[, 2],
                          x$tcfs[, 3])
  infer_tab[, 3:7] <- signif(infer_tab[, 3:7], digits = digits)
  colnames(infer_tab) <- c("Covariate(s) Values", "Method", "Threshold 1",
                           "Threshold 2", "TCF 1", "TCF 2", "TCF 3")
  cat("Covariate-specific optimal pair of thresholds: \n")
  print(infer_tab, quote = FALSE, right = TRUE, na.print = "--",
        row.names = FALSE, ...)
  cat("\n")
  if (!is.null(x$thres3_se)) {
    se_tab <- data.frame(labels, x$thres3_se$Method, x$thres3_se$threshold_1,
                         x$thres3_se$threshold_2)
    se_tab[, 3:4] <- signif(se_tab[, 3:4], digits = digits)
    colnames(se_tab) <- c("Covariate(s) Values", "Method", "SE. Threshold 1",
                          "SE. Threshold 2")
    cat("Standard errors of Covariate-specific optimal pair of thresholds: \n")
    print(se_tab, quote = FALSE, right = TRUE, na.print = "--",
          row.names = FALSE, ...)
    cat("\n")
  }
  invisible(x)
}
