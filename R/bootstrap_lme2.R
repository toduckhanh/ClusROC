####==========================================================================####
## This file consists of functions for nonparametric bootstrap for              ##
## the cluster-effect models under Box-Cox transformation                       ##
## Date: 27/03/2021																															##
####==========================================================================####

#' @import snow
#' @import doSNOW
#' @import foreach

boot_lme2 <- function(out_lme2, data, z, B = 250, type = c("cluster", "stratified"), boxcox = TRUE,
                      parallel = FALSE, ncpus = ifelse(parallel, 2, NULL)){
  if(isFALSE(inherits(out_lme2, "lme2"))) stop("out_lme2 was not from lme2()!")
  levl.class <- as.character(attr(out_lme2$terms, "levl.class"))
  data[, out_lme2$name.class] <- factor(data[, out_lme2$name.class], levels = levl.class)
  data_list <- split(data, data[, out_lme2$name.clust])
  fixed.formula <- as.formula(paste(out_lme2$name.test, "~", out_lme2$name.covars))
  type <- match.arg(type)
  # re-sampling process
  bts_func <- function(B, data_list, fixed.formula, name.class, name.clust, levl.class, type, boxcox, z){
    flag <- 0
    while(flag == 0){
      flag_data <- 0
      while(flag_data == 0){
        if(type == "cluster"){
          n_cluster <- length(data_list)
          id <- sample(1:n_cluster, n_cluster, replace = TRUE)
          data_list_bts <- data_list[id]
          nc_b <- sapply(data_list_bts, nrow)
          data_bts <- do.call(rbind, data_list_bts)
          data_bts[, name.clust] <- rep(1:length(id), nc_b)
        } else{
          n_c <- sapply(data_list, nrow)
          id_n_c <- sapply(n_c, function(y) sample(1:y, size = y, replace = TRUE))
          data_list_bts <- mapply(function(x, y) x[y,], x = data_list, y = id_n_c, SIMPLIFY = FALSE)
          data_bts <- do.call(rbind, data_list_bts)
        }
        if(isFALSE(0 %in% table(data_bts[, name.class]))) flag_data <- 1
      }
      out <- try(lme2(fixed.formula = fixed.formula, name.class = name.class, name.clust = name.clust,
                      data = data_bts, levl.class = levl.class, boxcox = boxcox, apVar = FALSE, trace = FALSE),
                 silent = TRUE)
      if(!inherits(out, "try-error")) {
        flag <- check_order(out$est_para, z = z, n_class = length(levl.class), n_p = out$n_p)
        if(boxcox){
          flag <- flag*check_sign(out$est_para, z = z, length(levl.class), n_p = out$n_p)
        }
      }
    }
    return(out$est_para)
  }
  if(!parallel){
    n_par <- length(out_lme2$est_para)
    out_boot <- matrix(nrow = n_par, ncol = B)
    for(i in 1:B){
      out_boot[,i] <- bts_func(B = B, data_list = data_list, fixed.formula = fixed.formula,
                               name.class = out_lme2$name.class, name.clust = out_lme2$name.clust,
                               levl.class = levl.class, type = type, boxcox = boxcox, z = z)
    }
  } else{
    cl <- makeCluster(rep("localhost", ncpus), type = "SOCK")
    registerDoSNOW(cl)
    clusterEvalQ(cl, list(library(ClusROC)))
    name.class <- out_lme2$name.class
    name.clust <- out_lme2$name.clust
    out_boot <- parSapply(cl, X = 1:B, FUN = bts_func, data_list = data_list, fixed.formula = fixed.formula,
                          name.class = name.class, name.clust = name.clust, levl.class = levl.class,
                          type = type, boxcox = boxcox, z = z)
    stopCluster(cl)
  }
  return(out_boot)
}
