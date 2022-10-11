boxcox_trans <- function(x, lambda) {
  if (lambda != 0) {
    y <- (x^lambda - 1) / lambda
  } else {
    y <- log(x)
  }
  return(y)
}

boxcox_trans_back <- function(x, lambda) {
  if (lambda != 0) {
    y <- (x * lambda + 1)^(1 / lambda)
  } else {
    y <- exp(x)
  }
  return(y)
}

probit <- function(x) qnorm(x)

make_data <- function(out_clus_lme, newdata, n_p) {
  if (isFALSE(inherits(out_clus_lme, "clus_lme"))) {
    stop("out_clus_lme was not from clus_lme()!")
  }
  terms <- out_clus_lme$terms
  fm <- delete.response(terms)
  class_disease <- factor(attr(terms, "levl_class"),
                          levels = attr(terms, "levl_class"))
  id_coef <- c(seq(1, out_clus_lme$n_coef - 2, by = 3),
               seq(2, out_clus_lme$n_coef - 1, by = 3),
               seq(3, out_clus_lme$n_coef, by = 3))
  if (n_p == 1) {
    newdata <- data.frame(class_disease)
    names(newdata) <- as.character(attr(fm, "variable"))[2]
    dt <- model.matrix(fm, data = newdata)[, id_coef]
    res <- list()
    res[[1]] <- dt
  } else {
    nt <- nrow(newdata)
    newdata_list <- split(newdata, 1:nt)
    fm_list <- lapply(newdata_list, function(x) {
      y <- data.frame(class_disease,
                      do.call(rbind, replicate(3, x, simplify = FALSE)))
      names(y)[1] <- as.character(attr(fm, "variable"))[2]
      return(y)
    })
    res <- lapply(fm_list, function(x) {
      m <- model.frame(fm, x, xlev = attr(terms, "xlevels"))
      x <- model.matrix(fm, m,
                        contrasts.arg = attr(terms, "contrasts"))[, id_coef]
      return(x)
    })
  }
  return(res)
}

# Check ordering of multiple points
check_mu_order <- function(z, par, n_p) {
  # z: design matrix
  beta_fit <- par[1:(3 * n_p)]
  mu_est <- sapply(z, function(x) x %*% beta_fit)
  status <- apply(mu_est, 2, function(x) (x[1] < x[2]) * (x[2] < x[3]))
  if (n_p == 1) {
    z_new <- z_del <- z
    mu_del <- mu_est
  } else {
    z_new <- z[status != 0]
    z_del <- z[status == 0]
    mu_del <- mu_est[, status == 0]
  }
  res <- list(mu_est = mu_est, status = status, z_new = z_new, z_del = z_del,
              mu_del = mu_del)
  return(res)
}

### The following functions are used for bootstrap process
check_order <- function(par, z, n_class, n_p) { # Check ordering of one point
  beta_fit <- par[1:(n_class * n_p)]
  mu_est <- z %*% beta_fit
  flag1 <- mu_est[1] - mu_est[2] < 0
  flag <- flag1
  if (n_class == 3) {
    flag2 <- mu_est[2] - mu_est[3] < 0
    flag <- flag1 * flag2
  }
  return(flag)
}

## only for Box-Cox transformation
check_sign <- function(par, z, n_class, n_p) {
  beta_fit <- par[1:(n_class * n_p)]
  mu_est <- z %*% beta_fit
  flag1 <- mu_est[1] > 0
  flag2 <- mu_est[2] > 0
  flag <- flag1 * flag2
  if (n_class == 3) {
    flag3 <- mu_est[3] > 0
    flag <- flag * flag3
  }
  return(flag)
}
