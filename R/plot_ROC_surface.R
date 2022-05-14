####==========================================================================####
## This file consists of functions for estimating and ploting the ROC surface   ##
## Date: 26/03/2021																															##
####==========================================================================####

#' @import rgl
#' @import utils
#' @import graphics

shade.ellips <- function(orgi, sig, lev){
  t1 <- sig[2, 2]
  sig[2, 2] <- sig[3, 3]
  sig[3, 3] <- t1
  t1 <- sig[1, 2]
  sig[1, 2] <- sig[1, 3]
  sig[1, 3] <- t1
  sig[lower.tri(sig)] <- sig[upper.tri(sig)]
  ellips <- ellipse3d(sig, centre = orgi[c(1,3,2)], t = sqrt(qchisq(lev, 3)))
  return(ellips)
}

### ---- Plot the ROC surface under normal assumption ----
#' @title Plot an estimated covariate-specific ROC surface for clustered data.
#'
#' @description \code{ROCsurface} estimates and makes a 3D plot of a covariate-specific ROC surface for a continuous diagnostic test, in a clustered design, with three ordinal groups.
#'
#' @param out_lme2  an object of class "lme2", a result of \code{\link{lme2}} call.
#' @param newdata  a data frame with 1 row (containing specific value(s) of covariate(s)) in which to look for variables with which to estimate covariate-specific ROC. In absence of covariate, no values have to be specified.
#' @param step.tcf  number: increment to be used in the grid for \eqn{p1 = tcf1} and \eqn{p3 = tcf3}.
#' @param main  the main title for plot.
#' @param file.name  	File name to create on disk.
#' @param ellips  a logical value. If set to \code{TRUE}, the function adds an ellipsoidal confidence region for TCFs (True Class Fractions), at a specified pair of values for the thresholds, to the plot of estimated covariate-specific ROC surface.
#' @param thresholds  a specified pair of thresholds, used to construct the ellipsoidal confidence region for TCFs.
#' @param ci.level  a confidence level to be used for constructing the ellipsoidal confidence region; default is 0.95.
#'
#' @details
#' This function implements a method in To et al. (2022) for estimating covariate-specific ROC surface of a continuous diagnostic test in a clustered design, with three ordinal groups. The estimator is based on the results from \code{\link{lme2}} with REML approach.
#'
#' Before performing estimation, a check for the monotone ordering assumption is performed. This means that, for the fixed values of covariates, three predicted mean values for test results in three diagnostic groups are compared. If the assumption is not meet, the covariate-specific ROC surface at the values of covariates is not estimated.
#'
#' The ellipsoidal confidence region for TCFs at a given pair of thresholds, if required, is constructed by using normal approximation and is plotted in the ROC surface space. The confidence level (default) is 0.95. Note that, if the Box-Cox transformation is applied for the linear mixed-effect model, the pair of thresholds must be input in the original scale. If the constructed confidence region for TCFs is outside the unit cube, a probit transformation will be automatically applied to obtain an appropriate confidence region, which is inside the unit cube (see Bantis et. al., 2017).
#'
#' @return \code{ROCsurface} returns a 3D \code{rgl} plot of the estimated covariate-specific ROC surface.
#'
#' @references
#' Bantis, L. E., Nakas, C. T., Reiser, B., Myall, D., and Dalrymple-Alford, J. C. (2017).
#' ``Construction of joint confidence regions for the optimal true class fractions of Receiver Operating Characteristic (ROC) surfaces and manifolds''. \emph{Statistical methods in medical research}, \bold{26}, 3, 1429-1442.
#'
#' To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2022)
#' ``Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data''. \emph{Statistical Methods in Medical Research}, DOI: 10.1177/09622802221089029.
#'
#' @examples
#' \dontrun{
#' data(data_3class)
#' ## One covariate
#' out1 <- lme2(fixed.formula = Y ~ X1, name.class = "D", name.clust = "id_Clus",
#'              data = data_3class)
#'
#' ### plot only covariate-specific ROC surface
#' ROCsurface(out_lme2 = out1, newdata = data.frame(X1 = 1))
#'
#' ### plot covariate-specific ROC surface and a 95% ellipsoidal confidence region for TCFs
#' ROCsurface(out_lme2 = out1, newdata = data.frame(X1 = 1), ellips = TRUE,
#'            thresholds = c(0.9, 3.95))
#'
#' ## Two covariates
#' out2 <- lme2(fixed.formula = Y ~ X1 + X2, name.class = "D", name.clust = "id_Clus",
#'              data = data_3class)
#'
#' ### plot only covariate-specific ROC surface
#' ROCsurface(out_lme2 = out2, newdata = data.frame(X1 = 1, X2 = 1))
#'
#' ### plot covariate-specific ROC surface and a 95% ellipsoidal confidence region for TCFs
#' ROCsurface(out_lme2 = out2, newdata = data.frame(X1 = 1, X2 = 1), ellips = TRUE,
#'            thresholds = c(0.9, 3.95))
#' }
#'
#' @export
ROCsurface <- function(out_lme2, newdata, step.tcf = 0.01, main = NULL, file.name = NULL,
                       ellips = FALSE, thresholds = NULL, ci.level = ifelse(ellips, 0.95, NULL)){
  # define parameters
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  if(out_lme2$n_coef/out_lme2$n_p != 3) stop("There is not a case of three-class setting!")
  n_p <- out_lme2$n_p
  if(n_p == 1) {
    if(!missing(newdata)) {
      if(!is.null(newdata)) warning("Sepecified value(s) of covariate(s) are not used!", call. = FALSE)
    }
    newdata <- NULL
  }else {
    if(missing(newdata)) stop("Please input a data frame including specific value(s) of covariate(s).")
    if(is.null(newdata)) stop("Please input a data frame including specific value(s) of covariate(s).")
    if(!inherits(newdata, "data.frame") | nrow(newdata) != 1)
      stop("The number of rows in newdata must be equal 1.")
    if(any(is.na(newdata))) stop("NA value(s) not allowed!")
  }
  if(ellips){
    if(is.null(thresholds)) stop("Need to assign the pair of thresholds!")
    else {
      if(!inherits(thresholds, "numeric") | length(thresholds) != 2) stop("Please input a the pair of thresholds!")
      else {
        if(thresholds[1] > thresholds[2]) stop("The 1st threshold needs to less than 2nd threshold!")
        else {
          if(is.null(out_lme2$vcov_sand)) stop("The estimated covariance matrix of parameters was missing!")
          if(any(is.na(out_lme2$vcov_sand))) stop("There are NA values in the estimated covariance matrix of parameters. Unable to estimate standard error of TCFs.")
        }
      }
    }
  }
  ## main
  par_model <- out_lme2$est_para
  beta_d <- par_model[1:(3*n_p)]
  sigma_d <- sqrt(par_model[(3*n_p + 2):(3*n_p + 4)]^2 + par_model[(3*n_p + 1)]^2)
  Z <- make_data(out_lme2, newdata, n_p)
  mu_d <- Z[[1]] %*% beta_d
  if((mu_d[1] < mu_d[2])*(mu_d[2] < mu_d[3]) == 0) stop("The monotone ordering assumption DOES NOT hold for the value(s) of the covariate(s)")
  a12 <- sigma_d[2]/sigma_d[1]
  a32 <- sigma_d[2]/sigma_d[3]
  b12 <- (mu_d[1] - mu_d[2])/sigma_d[1]
  b32 <- (mu_d[3] - mu_d[2])/sigma_d[3]
  # estimate the ROC surface by given p1 = TCF1 and p3 = TCF3
  p1 <- p3 <- seq(0, 1, by = step.tcf)
  dt_p1_p3 <- expand.grid(x = p1, y = p3)
  rocs <- apply(dt_p1_p3, 1, function(x){
    if(qnorm(x[1], mean = mu_d[1], sd = sigma_d[1]) <= qnorm(1 - x[2], mean = mu_d[3], sd = sigma_d[3])){
      res <- pnorm((qnorm(1 - x[2]) + b32)/a32) - pnorm((qnorm(x[1]) + b12)/a12)
    } else res <- 0
    return(res)
  })
  fit <- cbind(dt_p1_p3, rocs)
  colnames(fit) <- c("TCF1", "TCF3", "TCF2")
  out <- list()
  out$fit <- fit
  tcf1 <- matrix(p1, length(p1), length(p1), byrow = FALSE)
  tcf3 <- matrix(p3, length(p1), length(p1), byrow = TRUE)
  tcf2 <- matrix(rocs, length(p1), length(p1), byrow = FALSE)
  open3d()
  my.userMatrix <- rbind(c(-0.8370321, -0.5446390, -0.0523976, 0),
                         c(0.1272045, -0.2868422, 0.9494949, 0),
                         c(-0.5321618, 0.7880925, 0.3093767, 0),
                         c(0, 0, 0, 1))
  par3d(windowRect = 50 + c(0, 0, 640, 640), userMatrix = my.userMatrix)
  if(is.null(main)) main <- "Covariate-specific ROC surface"
  # plot3d(0, 0, 0, type = "n", box = FALSE, xlab = "TCF 1", ylab = "TCF 3", zlab = "TCF 2",
  #        xlim = c(0,1), ylim = c(0,1), zlim = c(0,1))
  plot3d(0, 0, 0, type = "n", box = FALSE, xlab = "", ylab = "", zlab = "",
         xlim = c(0,1), ylim = c(0,1), zlim = c(0,1), axes = FALSE)
  axes3d(edges = c("x--", "y--", "z--"))
  mtext3d("TCF 1", "x--", line = 2, at = 0.35)
  mtext3d("TCF 2", "z--", line = 4, at = 0.55)
  mtext3d("TCF 3", "y--", line = 4, at = 0.15, level = 2)
  bgplot3d({
    plot.new()
    title(main = main, line = 1)
  })
  surface3d(tcf1, tcf3, tcf2, col = "gray40", alpha = 0.5)
  if(ellips){
    if(inherits(thresholds, "numeric")){
      tcfs_ellips <- TCF_normal(par = par_model, z = Z[[1]], thresholds = thresholds, n_p = n_p,
                                boxcox = out_lme2$boxcox)
      vcov_tcfs <- TCF_normal_vcov(par_model = par_model, z = Z[[1]], thresholds = thresholds,
                                   vcov_par_model = out_lme2$vcov_sand, n_p = n_p, fixed = TRUE,
                                   boxcox = out_lme2$boxcox)
    }
    ellip.tcf <- shade.ellips(orgi = tcfs_ellips, sig = vcov_tcfs, lev = ci.level)
    if(max(c(ellip.tcf$vb[1,], ellip.tcf$vb[2,], ellip.tcf$vb[3,])) > 1 |
       min(c(ellip.tcf$vb[1,], ellip.tcf$vb[2,], ellip.tcf$vb[3,])) < 0){
      tcfs_ellips_prob <- probit(tcfs_ellips)
      jac_prob <- suppressWarnings(jacobian(func = probit, x = tcfs_ellips))
      if(!any(is.na(jac_prob))){
        vcov_tcfs_prob <- jac_prob %*% vcov_tcfs %*% jac_prob
        ellip.tcf_prob <- shade.ellips(orgi = tcfs_ellips_prob, sig = vcov_tcfs_prob, lev = ci.level)
        ellip.tcf <- ellip.tcf_prob
        ellip.tcf$vb[1:3,] <- pnorm(ellip.tcf_prob$vb[1:3,])
        out$massge_probit <- "A probit transformation was applied to guarantee the confidence region inside the unit cube.\n"
        cat(out$massge_probit)
      } else{
        out$massge_probit0 <- "Could not use probit transformation, due to point estimates of TCF 1/TCF 2/TCF 3 are very close to 0 or 1.\n"
        cat(out$massge_probit0)
        out$massge_probit1 <- "The confidence region of (TCF 1, TCF 2, TCF 3) is truncated to (0, 1).\n"
        cat(out$massge_probit1)
        ellip.tcf$vb[1, ellip.tcf$vb[1,] > 1] <- 1
        ellip.tcf$vb[2, ellip.tcf$vb[2,] > 1] <- 1
        ellip.tcf$vb[3, ellip.tcf$vb[3,] > 1] <- 1
        ellip.tcf$vb[1, ellip.tcf$vb[1,] < 0] <- 0
        ellip.tcf$vb[2, ellip.tcf$vb[2,] < 0] <- 0
        ellip.tcf$vb[3, ellip.tcf$vb[3,] < 0] <- 0
      }
    }
    plot3d(ellip.tcf, box = FALSE, col = "green", alpha = 0.5, xlim = c(0,1),
           ylim = c(0, 1), zlim = c(0, 1), xlab = " ", ylab = " ", zlab = " ",
           add = TRUE)
    plot3d(tcfs_ellips[1], tcfs_ellips[3], tcfs_ellips[2], type = "s", col = "red",
           radius = 0.01, add = TRUE)
  }
  # play3d(spin3d(axis = c(0, 0, 1), rpm = 12.25), duration = 2)
  # play3d(spin3d(axis = c(0, 1, 0), rpm = 0.3), duration = 2)
  if(!is.null(file.name)){
    if(!grepl(".png", file.name)) file.name <- paste0(file.name, ".png")
    rgl.snapshot(file.name)
  }
  invisible(out)
}

