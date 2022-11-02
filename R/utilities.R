##
check_fixed_formula <- function(fixed_formula, fun_call, names_vars) {
  res <- NULL
  if (missing(fixed_formula) || !inherits(fixed_formula, "formula") ||
      length(fixed_formula) != 3) {
    stop("agrument \"fixed_formula\" must be a formula of the form \"resp ~ pred\"")
  } else {
    all_names_mf <- all.vars(fixed_formula)
    name_test <- all_names_mf[1]
    names_covars <- all_names_mf[-1]
    if (!is.element(name_test, names_vars)) {
      stop(paste(name_test, " is not included in ", deparse(fun_call$data)))
    }
    if (!all(is.element(names_covars, names_vars))) {
      stop(paste("One or some covariates are not found in ",
                 deparse(fun_call$data)))
    }
    res <- list(name_test = name_test, names_covars = names_covars)
  }
  return(res)
}

##
check_class <- function(name_class, names_vars, data) {
  res <- NULL
  if (missing(name_class) || !inherits(name_class, "character") ||
      length(name_class) != 1 || !is.element(name_class, names_vars)) {
    stop("agrument \"name_class\" was either missing or wrong name!")
  } else {
    n_class <- length(unique(data[, name_class]))
    if (n_class != 3) {
      stop("agrument \"name_class\" must have 3 levels or classes!")
    }
    res <- list(n_class = n_class)
  }
  return(res)
}

##
check_clust <- function(name_clust, names_vars, name_test,
                        names_covars, name_class) {
  res <- NULL
  if (missing(name_clust) || !inherits(name_clust, "character") ||
      length(name_clust) != 1 || !is.element(name_clust, names_vars)) {
    stop("agrument \"name_clust\" was either missing or wrong name!")
  } else {
    if (is.element(name_clust, c(name_test, names_covars, name_class))) {
      stop("agrument \"name_clust\" cannot be the name of neither test, covariates nor classes!")
    }
  }
  return(res)
}

##
check_levl_class <- function(trace, levl_class, temp_levl) {
  if (is.null(levl_class)) {
    if (trace) {
      cat("The ordered levels of classes are specified by the order of \n averages of the test values for each class:\n")
      cat(paste(temp_levl, collapse = " < "), "\n")
    }
    levl_class <- temp_levl
  } else {
    if (any(is.na(levl_class)) || !inherits(levl_class, "character") ||
        length(levl_class) != 3) {
      stop("agrument levl_class must be a character vector with length 3 without NA.")
    }
    if (all(levl_class == temp_levl)) {
      if (trace) {
        cat("The orders of \"levl_class\" are the same as \n the orders of averages of tests values for each class:\n")
        cat(paste(levl_class, collapse = " < "), "\n")
      }
    } else {
      if (trace) {
        cat("The orders of \"levl_class\" are not the same as \n the orders of averages of tests values for each class\n")
        cat("The correct one should be:\n")
        cat(paste(temp_levl, collapse = " < "), "\n")
      }
      levl_class <- temp_levl
    }
  }
  return(list(levl_class = levl_class))
}

##
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

check_newdata_vus <- function(fixed_formula, newdata, n_p) {
  if (n_p == 1) {
    if (!missing(newdata)) {
      if (!is.null(newdata)) {
        cat("Warning message:\n")
        cat("Sepecified value(s) of covariate(s) are not used!\n")
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
    newdata <- model.frame(formula = delete.response(terms(fixed_formula)),
                           newdata)
    attr(newdata, "terms") <- NULL
    newdata <- as.data.frame(newdata)
  }
  return(list(newdata = newdata))
}

check_newdata_roc <- function(fixed_formula, newdata, n_p) {
  if (n_p == 1) {
    if (!missing(newdata)) {
      if (!is.null(newdata)) {
        cat("Warning message:\n")
        cat("Sepecified value(s) of covariate(s) are not used!\n")
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
    if (!inherits(newdata, "data.frame") || nrow(newdata) != 1) {
      stop("The number of rows in newdata must be equal 1.")
    }
    if (any(is.na(newdata))) {
      stop("NA value(s) not allowed!")
    }
    newdata <- model.frame(formula = delete.response(terms(fixed_formula)),
                           newdata)
    attr(newdata, "terms") <- NULL
    newdata <- as.data.frame(newdata)
  }
  return(list(newdata = newdata))
}

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

new_data_check <- function(res_check) {
  mess_order <- NULL
  if (all(res_check$status == 0)) {
    stop("The assumption of montone ordering DOES NOT hold for all the value(s) of the covariate(s)")
  }
  if (any(res_check$status == 0)) {
    mess_order <- paste("The assumption of montone ordering DOES NOT hold for some points. The points number:",
                        paste(which(res_check$status == 0), collapse = ", "),
                        "are deleted from analysis!")
    message(mess_order)
  }
  return(list(z_new = res_check$z_new, mess_order = mess_order))
}

check_ellip_roc <- function(ellips, thresholds, vcov_sand) {
  res <- NULL
  if (ellips) {
    if (is.null(thresholds)) {
      stop("Need to assign the pair of thresholds!")
    } else {
      if (!inherits(thresholds, "numeric") || length(thresholds) != 2) {
        stop("Please input a the pair of thresholds!")
      } else {
        if (thresholds[1] > thresholds[2]) {
          stop("The 1st threshold needs to less than 2nd threshold!")
        } else {
          if (is.null(vcov_sand)) {
            stop("The estimated covariance matrix of parameters was missing!")
          }
          if (any(is.na(vcov_sand))) {
            stop("There are NA values in the estimated covariance matrix of parameters. Unable to estimate standard error of TCFs.")
          }
        }
      }
    }
  }
  return(res)
}

##
get_method_opt_thres3 <- function(method) {
  if (missing(method)) {
    stop("Please, input the selection method(s)!")
  }
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
  return(list(methodtemp = methodtemp, mmethod = mmethod))
}

get_labels_opt_thres3 <- function(n_p, newdata, names_labels) {
  if (n_p == 1) {
    n_x <- 1
    point_labels <- "Intercept"
    if (missing(names_labels)) {
      names_labels <- " "
    }
  }
  if (n_p == 2) {
    n_x <- nrow(newdata)
    point_labels <- apply(newdata, 1, function(y) paste0(y))
    if (missing(names_labels)) {
      names_labels <- "Value(s) of covariate:"
    }
  }
  if (n_p > 2) {
    n_x <- nrow(newdata)
    point_labels <- apply(newdata, 1, function(y) {
      paste0("(", paste(y, collapse = ", "), ")")
    })
    if (missing(names_labels)) {
      names_labels <- "Value(s) of covariates:"
    }
  }
  return(list(n_x = n_x, point_labels = point_labels,
              names_labels = names_labels))
}

##
check_apvar_opt_thres3 <- function(ap_var, boxcox, vcov_sand) {
  bootstrap <- FALSE
  if (ap_var) {
    if (!boxcox) { ## asymptotic variance under normal distribution
      if (is.null(vcov_sand)) {
        stop("The estimated covariance matrix of parameters was missing!")
      }
      if (any(is.na(vcov_sand))) {
        stop("There are NA values in the estimated covariance matrix of parameters. Unable to estimate standard error.")
      }
    } else {
      bootstrap <- TRUE
    }
  }
  return(bootstrap)
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
