test_that("clus_opt_thres3 works with no covariate", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ 1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_opt_thres3(method = c("GYI", "MV", "CtP"),
                             out_clus_lme = out_md, ap_var = TRUE),
    class = "clus_opt_thres3")
})

test_that("clus_opt_thres3 works with one contiuous covariate", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_opt_thres3(method = c("GYI", "MV", "CtP"),
                             out_clus_lme = out_md,
                             newdata = data.frame(X1 = 1), ap_var = TRUE),
    class = "clus_opt_thres3")
})

test_that("clus_opt_thres3 works with two covariates", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_opt_thres3(method = c("GYI", "MV", "CtP"),
                             out_clus_lme = out_md,
                             newdata = data.frame(X1 = c(0.5, 1), X2 = c(0, 0)),
                             ap_var = TRUE),
    class = "clus_opt_thres3")
})

test_that("clus_opt_thres3 works with Box-Cox transformation,
          standard bootstrap process with 250 bootstrap replicates", {
  data(data_3class_bcx)
  out_md <- clus_lme(fixed_formula = Y ~ X, name_class = "D",
                     name_clust = "id_Clus", data = data_3class_bcx,
                     boxcox = TRUE)
  expect_s3_class(
    object = clus_opt_thres3(method = c("GYI", "MV", "CtP"),
                             out_clus_lme = out_md,
                             newdata = data.frame(X = c(1, 1.2)),
                             ap_var = TRUE),
    class = "clus_opt_thres3")
})

test_that("clus_opt_thres3 works with Box-Cox transformation,
          parallel bootstrap process with 250 bootstrap replicates", {
  data(data_3class_bcx)
  out_md <- clus_lme(fixed_formula = Y ~ X, name_class = "D",
                     name_clust = "id_Clus", data = data_3class_bcx,
                     boxcox = TRUE)
  expect_s3_class(
    object = clus_opt_thres3(method = c("GYI", "MV", "CtP"),
                             out_clus_lme = out_md,
                             newdata = data.frame(X = c(1, 1.2)),
                             ap_var = TRUE,
                             control = list(parallel = TRUE, ncpus = 4)),
    class = "clus_opt_thres3")
})

test_that("clus_opt_thres3 does not work if not input newdata", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_error(
    object = clus_opt_thres3(method = c("GYI", "MV", "CtP"),
                             out_clus_lme = out_md))
})

