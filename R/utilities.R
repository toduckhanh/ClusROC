boxcox_trans <- function(x, lambda){
  if(lambda != 0) y <- (x^lambda - 1)/lambda
  else y <- log(x)
  return(y)
}

boxcox_trans_back <- function(x, lambda){
  if(lambda != 0) y <- (x*lambda + 1)^(1/lambda)
  else y <- exp(x)
  return(y)
}

probit <- function(x) qnorm(x)

make_data <- function(out_lme2, newdata, n_p){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  Terms <- out_lme2$terms
  fm <- delete.response(Terms)
  class.disease <- factor(attr(Terms, "levl.class"), levels = attr(Terms, "levl.class"))
  id_coef <- c(seq(1, out_lme2$n_coef - 2, by = 3), seq(2, out_lme2$n_coef - 1, by = 3),
               seq(3, out_lme2$n_coef, by = 3))
  if(n_p == 1){
    newdata <- data.frame(class.disease)
    names(newdata) <- as.character(attr(fm, "variable"))[2]
    dt <- model.matrix(fm, data = newdata)[,id_coef]
    res <- list()
    res[[1]] <- dt
  }
  else {
    n_dt <- nrow(data)
    newdata_list <- split(newdata, 1:nrow(newdata))
    fm_list <- lapply(newdata_list, function(x) {
      y <- data.frame(class.disease, do.call(rbind, replicate(3, x, simplify = FALSE)))
      names(y)[1] <- as.character(attr(fm, "variable"))[2]
      return(y)
    })
    res <- lapply(fm_list, function(x){
      m <- model.frame(fm, x, xlev = attr(Terms, "xlevels"))
      X <- model.matrix(fm, m, contrasts.arg = attr(Terms, "contrasts"))[,id_coef]
      return(X)
    })
  }
  return(res)
}

# Check ordering of multiple points
check_mu_order <- function(Z, par, n_p){
  # Z: design matrix
  beta_fit <- par[1:(3*n_p)]
  mu_est <- sapply(Z, function(x) x %*% beta_fit)
  status <- apply(mu_est, 2, function(x) (x[1] < x[2])*(x[2] < x[3]))
  if(n_p == 1) {
    Z_new <- Z_del <- Z
    mu_del <- mu_est
  }
  else {
    Z_new <- Z[status != 0]
    Z_del <- Z[status == 0]
    mu_del <- mu_est[, status == 0]
  }
  res <- list(mu_est = mu_est, status = status, Z_new = Z_new, Z_del = Z_del, mu_del = mu_del)
  return(res)
}

### The following functions are used for bootstrap process
check_order <- function(par, z, n_class, n_p){ # Check ordering of one point
  beta_fit <- par[1:(n_class*n_p)]
  mu_est <- z %*% beta_fit
  flag1 <- mu_est[1] - mu_est[2] < 0
  flag <- flag1
  if(n_class == 3){
    flag2 <- mu_est[2] - mu_est[3] < 0
    flag <- flag1*flag2
  }
  return(flag)
}

check_sign <- function(par, z, n_class, n_p){ ## only for Box-Cox transformation
  beta_fit <- par[1:(n_class*n_p)]
  mu_est <- z %*% beta_fit
  flag1 <- mu_est[1] > 0
  flag2 <- mu_est[2] > 0
  flag <- flag1*flag2
  if(n_class == 3){
    flag3 <- mu_est[3] > 0
    flag <- flag*flag3
  }
  return(flag)
}






