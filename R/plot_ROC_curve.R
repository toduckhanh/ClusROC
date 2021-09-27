####==========================================================================####
## This file consists of functions for estimating and plotting the ROC curve     ##
## Date: 27/05/2021																															##
####==========================================================================####


ROCfunc <- function(par, Z, pp, n_p){
  beta_d <- matrix(par[1:(2*n_p)], ncol = 2, nrow = n_p, byrow = FALSE)
  sigma_d <- sqrt(par[(2*n_p + 2):(2*n_p + 3)]^2 + par[(2*n_p + 1)]^2)
  mu_d <- Z %*% beta_d
  return(pnorm((mu_d[2] - mu_d[1] + sigma_d[1]*qnorm(pp))/sigma_d[2]))
}

#' @export
ROCcurve <- function(out_lme2, x.val, step.fpf = 0.01, plot = TRUE, ci.roc = FALSE,
                     ci.level = ifelse(ci.roc, 0.95, NULL), trans.roc = c("probit", "logit"),
                     main = NULL, file.name = NULL){
  # ellips = FALSE, threshold = NULL, ci.level = ifelse(ellips, 0.95, NULL)
  # define parameters
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  if(out_lme2$n_coef/out_lme2$n_p != 2) stop("There is not a case of two-class setting!")
  n_p <- out_lme2$n_p
  par_model <- out_lme2$est_para
  Z <- c(1, x.val)
  p1 <- seq(1e-9, 1 - 1e-9, by = step.fpf)
  roc <- ROCfunc(par_model, Z, p1, n_p)
  fit <- data.frame(FPF = p1, TPF = roc, name = rep("Covariate-specific ROC curve", length(p1)))
  if(ci.roc){
    roc_jac <- jacobian(func = ROCfunc, x = par_model, Z = Z, pp = p1, n_p = n_p)
    roc_se <- sqrt(apply(roc_jac, 1, function(x, y) as.numeric(x %*% y %*% x), y = out_lme2$vcov_sand))
    trans.roc <- match.arg(trans.roc)
    delta <- switch(trans.roc,
                    probit = qnorm(roc),
                    logit = qlogis(roc))
    delta_grad <- switch(trans.roc,
                         probit = 1/(dnorm(delta)),
                         logit = 1/(roc*(1 - roc)))
    delta_sd <- delta_grad*roc_se
    delta_ci <- t(mapply(FUN = function(x, y) x + c(-1, 1)*qnorm((1 + ci.level)/2)*y,
                         x = delta, y = delta_sd))
    roc.ci <- switch(trans.roc,
                     probit = pnorm(delta_ci),
                     logit = plogis(delta_ci))
    fit$lci <- roc.ci[,1]
    fit$uci <- roc.ci[,2]
  }
  out <- list()
  out$fit <- fit
  if(plot){
    if(ci.roc){
      pp <- ggplot(data = fit, aes(x = FPF, y = TPF, ymin = lci, ymax = uci)) +
        facet_grid( ~ name) +
        geom_line(colour = "blue") +
        geom_ribbon(alpha = 0.3) +
        geom_abline(intercept = 0, slope = 1, linetype = 2) +
        xlab("1 - Specificity") + ylab("Sensitivity") +
        xlim(0, 1) + ylim(0, 1) +
        theme_bw()
    } else{
      pp <- ggplot(out$fit, aes(x = FPF, y = TPF)) +
        facet_grid( ~ name) +
        geom_line(colour = "blue") +
        geom_abline(intercept = 0, slope = 1, linetype = 2) +
        xlab("1 - Specificity") + ylab("Sensitivity") +
        xlim(0, 1) + ylim(0, 1) +
        theme_bw()
    }
  }
  if(!is.null(file.name)){
    ggexport(pp, filename = file.name, ... = ...)
  }
  pp
}


