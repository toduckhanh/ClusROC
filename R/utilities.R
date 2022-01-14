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

make_data <- function(out_lme2, x.val, n_p){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  Terms <- out_lme2$terms
  fm <- delete.response(Terms)
  names.data <- as.character(attr(fm, "variable"))[-1]
  class.disease <- factor(attr(Terms, "levl.class"), levels = attr(Terms, "levl.class"))
  id_coef <- c(seq(1, out_lme2$n_coef - 2, by = 3), seq(2, out_lme2$n_coef - 1, by = 3),
               seq(3, out_lme2$n_coef, by = 3))
  if(n_p == 1){
    new_data <- data.frame(class.disease)
    names(new_data) <- names.data
    dt <- model.matrix(fm, data = new_data)[,id_coef]
    res <- dt
  }
  if(n_p == 2){
    n_x <- length(x.val)
    class.disease.2 <- rep(class.disease, n_x)
    x.val.2 <- rep(x.val, each = 3)
    new_data <- data.frame(class.disease.2, x.val.2)
    names(new_data) <- names.data
    dt <- model.matrix(fm, data = new_data)[,id_coef]
    res <- lapply(split(dt, rep(1:n_x, each = 3)), matrix, ncol = ncol(dt))
  }
  if(n_p > 2){
    n_vb <- attr(out_lme2$terms, "n_vb")
    if(!inherits(x.val, "matrix")) x.val <- matrix(x.val, ncol = n_vb, byrow = FALSE)
    n_x <- nrow(x.val)
    class.disease.2 <- rep(class.disease, n_x)
    x.val.2 <- apply(x.val, 2, function(x) rep(x, each = 3))
    new_data <- data.frame(class.disease.2, x.val.2)
    names(new_data) <- names.data
    dt <- model.matrix(fm, data = new_data)[,id_coef]
    res <- lapply(split(dt, rep(1:n_x, each = 3)), matrix, ncol = ncol(dt))
  }
  return(res)
}

check_mu_order <- function(Z, par, n_p){ # Check ordering of multiple points
  # Z: design matrix
  beta_fit <- par[1:(3*n_p)]
  if(n_p == 1) {
    mu_est <- beta_fit
    status <- (mu_est[1] < mu_est[2])*(mu_est[2] < mu_est[3])
    Z_new <- Z_del <- Z
    mu_del <- mu_est
  }
  if(n_p == 2){
    mu_est <- sapply(Z, function(x) x %*% beta_fit)
    status <- apply(mu_est, 2, function(x) (x[1] < x[2])*(x[2] < x[3]))
    Z_new <- Z[status != 0]
    Z_del <- Z[status == 0]
    mu_del <- mu_est[, status == 0]
  }
  if(n_p > 2){
    mu_est <- sapply(Z, function(x) x %*% beta_fit)
    status <- apply(mu_est, 2, function(x) (x[1] < x[2])*(x[2] < x[3]))
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






