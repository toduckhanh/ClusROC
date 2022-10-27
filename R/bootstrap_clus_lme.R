####========================================================================####
## This file consists of functions for nonparametric bootstrap for            ##
## the cluster-effect models under Box-Cox transformation                     ##
####========================================================================####

#' @import snow
#' @import doSNOW
#' @import foreach

boot_clus_lme <- function(out_clus_lme, z, n_boot = 250,
                          type = c("cluster", "stratified"), boxcox = TRUE,
                          parallel = FALSE, ncpus = ifelse(parallel, 2, NULL)) {
  if (isFALSE(inherits(out_clus_lme, "clus_lme"))) {
    stop("out_clus_lme was not from clus_lme()!")
  }
  levl_class <- as.character(attr(out_clus_lme$terms, "levl_class"))
  data <- out_clus_lme$data
  data[, out_clus_lme$name_class] <- factor(data[, out_clus_lme$name_class],
                                            levels = levl_class)
  data_list <- split(data, data[, out_clus_lme$name_clust])
  fixed_formula <- as.formula(paste(out_clus_lme$name_test, "~",
                                    out_clus_lme$name_covars))
  type <- match.arg(type)
  # re-sampling process
  bts_func <- function(n_boot, data_list, fixed_formula, name_class, name_clust,
                       levl_class, type, boxcox, z) {
    flag <- 0
    while (flag == 0) {
      flag_data <- 0
      while (flag_data == 0) {
        if (type == "cluster") {
          n_cluster <- length(data_list)
          id <- sample(1:n_cluster, n_cluster, replace = TRUE)
          data_list_bts <- data_list[id]
          nc_b <- sapply(data_list_bts, nrow)
          data_bts <- do.call(rbind, data_list_bts)
          data_bts[, name_clust] <- rep(1:n_cluster, nc_b)
        } else {
          n_c <- sapply(data_list, nrow)
          id_n_c <- sapply(n_c, function(y) {
            sample(1:y, size = y, replace = TRUE)
          })
          data_list_bts <- mapply(function(x, y) x[y, ], x = data_list,
                                  y = id_n_c, SIMPLIFY = FALSE)
          data_bts <- do.call(rbind, data_list_bts)
        }
        if (isFALSE(0 %in% table(data_bts[, name_class]))) flag_data <- 1
      }
      out <- try(clus_lme(fixed_formula = fixed_formula,
                          name_class = name_class, name_clust = name_clust,
                          data = data_bts, levl_class = levl_class,
                          boxcox = boxcox, ap_var = FALSE, trace = FALSE),
                 silent = TRUE)
      if (!inherits(out, "try-error")) {
        flag <- check_order(out$est_para, z = z, n_class = length(levl_class),
                            n_p = out$n_p)
        if (boxcox) {
          flag <- flag * check_sign(out$est_para, z = z, length(levl_class),
                                    n_p = out$n_p)
        }
      }
    }
    return(out$est_para)
  }
  if (!parallel) {
    n_par <- length(out_clus_lme$est_para)
    out_boot <- matrix(nrow = n_par, ncol = n_boot)
    for (i in 1:n_boot) {
      out_boot[, i] <- bts_func(n_boot = n_boot, data_list = data_list,
                                fixed_formula = fixed_formula,
                                name_class = out_clus_lme$name_class,
                                name_clust = out_clus_lme$name_clust,
                                levl_class = levl_class, type = type,
                                boxcox = boxcox, z = z)
    }
  } else {
    cl <- makeCluster(rep("localhost", ncpus), type = "SOCK")
    registerDoSNOW(cl)
    clusterEvalQ(cl, list(library(ClusROC)))
    name_class <- out_clus_lme$name_class
    name_clust <- out_clus_lme$name_clust
    out_boot <- parSapply(cl, X = 1:n_boot, FUN = bts_func,
                          data_list = data_list, fixed_formula = fixed_formula,
                          name_class = name_class, name_clust = name_clust,
                          levl_class = levl_class,
                          type = type, boxcox = boxcox, z = z)
    stopCluster(cl)
  }
  return(out_boot)
}
