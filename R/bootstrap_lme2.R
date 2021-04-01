####==========================================================================####
## This file consists of functions for nonparametric bootstrap for              ##
## the cluster-effect models under Box-Cox transformation                       ##
## Date: 27/03/2021																															##
####==========================================================================####

check_order <- function(par, x.val, n_class, n_p){
  beta_fit <- par[1:(n_class*n_p)]
  Z <- c(1, x.val)
  if(n_class == 2){
    flag <- (beta_fit[1:n_p] - beta_fit[(n_p + 1):(2*n_p)]) %*% Z < 0
  }
  if(n_class == 3){
    flag1 <- (beta_fit[1:n_p] - beta_fit[(n_p + 1):(2*n_p)]) %*% Z < 0
    flag2 <- (beta_fit[(n_p + 1):(2*n_p)] - beta_fit[(2*n_p + 1):(3*n_p)]) %*% Z < 0
    flag <- flag1*flag2
  }
  return(flag)
}

#' @import snow
#' @import doSNOW
#' @import foreach

boot_lme2 <- function(name.test, name.class, name.covars, name.clust, data, levl.class = NULL, x.val,
                      B = 250, type = c("cluster", "stratified"), boxcox = TRUE,
                      parallel = FALSE, ncpus = ifelse(parallel, 2, NULL), ...){
  form.mean <- as.formula(paste(name.test, "~", name.class))
  mean.temp <- aggregate(form.mean, FUN = mean, data = data)
  temp.levl <- mean.temp[order(mean.temp[,2]), 1]
  if(is.null(levl.class)){
    levl.class <- temp.levl
  } else{
    if(!all(levl.class == temp.levl)){
      levl.class <- temp.levl
    }
  }
  data[, name.class] <- factor(data[, name.class], levels = levl.class)
  data_list <- split(data, data[,name.clust])
  type <- match.arg(type)
  # re-sampling process
  bts_func <- function(B, data_list, name.test, name.class, name.covars, name.clust, levl.class, type, boxcox,
                       x.val){
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
          data_bts[,name.clust] <- rep(1:length(id), nc_b)
        } else{
          n_c <- sapply(data_list, nrow)
          id_n_c <- sapply(n_c, function(y) sample(1:y, size = y, replace = TRUE))
          data_list_bts <- mapply(function(x, y) x[y,], x = data_list, y = id_n_c, SIMPLIFY = FALSE)
          data_bts <- do.call(rbind, data_list_bts)
        }
        if(isFALSE(0 %in% table(data_bts[,name.class]))) flag_data <- 1
      }
      out <- try(lme2(name.test = name.test, name.class = name.class, name.covars = name.covars,
                      name.clust = name.clust, data = data_bts, levl.class = levl.class,
                      boxcox = boxcox, apVar = FALSE, trace = FALSE),
                 silent = FALSE)
      if(class(out) != "try-error") {
        print(out$est_para)
        print(length(levl.class))
        print(out$n_p)
        print(x.val)
        flag <- check_order(out$est_para, x.val = x.val, n_class = length(levl.class), n_p = out$n_p)
      }
    }
    return(out$est_para)
  }
  if(!parallel){
    if(missing(name.covars)) n_p <- 2*length(levl.class) + 1 + as.numeric(boxcox)
    else n_p <- (length(name.covars) + 2)*length(levl.class) + 1 + as.numeric(boxcox)
    out_boot <- matrix(nrow = n_p, ncol = B)
    for(i in 1:B){
      out_boot[,i] <- bts_func(B = B, data_list = data_list, name.test = name.test, name.class = name.class,
                               name.covars = name.covars, name.clust = name.clust, levl.class = levl.class,
                               type = type, boxcox = boxcox, x.val = x.val)
    }
  } else{
    cl <- makeCluster(rep("localhost", ncpus), type = "SOCK")
    registerDoSNOW(cl)
    clusterEvalQ(cl, list(library(ClusROC)))
    out_boot <- parSapply(cl, X = 1:B, FUN = bts_func, data_list = data_list, name.test = name.test,
                          name.class = name.class, name.covars = name.covars, name.clust = name.clust,
                          levl.class = levl.class, type = type, boxcox = boxcox, x.val = x.val)
    stopCluster(cl)
  }
  return(out_boot)
}
