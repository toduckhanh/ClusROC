test_that("clus_tcfs works with no covariate", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ 1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_tcfs(out_md, thresholds = c(1, 4)),
    class = "clus_tcfs")
})

test_that("clus_tcfs works with one contiuous covariate", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_tcfs(out_md, newdata = data.frame(X1 = 0.5),
                       thresholds = c(1, 4)),
    class = "clus_tcfs")
})

test_that("clus_tcfs works with two covariates", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_tcfs(out_md, newdata = data.frame(X1 = 0.5, X2 = 0),
                       thresholds = c(1, 4)),
    class = "clus_tcfs")
})

test_that("clus_tcfs works with Box-Cox transformation", {
  data("data_3class_bcx")
  out_md <- clus_lme(fixed_formula = Y ~ X, name_class = "D",
                     name_clust = "id_Clus", data = data_3class_bcx,
                     boxcox = TRUE)
  expect_s3_class(
    object = clus_tcfs(out_md, newdata = data.frame(X = 1),
                       thresholds = c(1, 4)),
    class = "clus_tcfs")
})

test_that("clus_tcfs does not work in case of one covariate without newdata", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_error(
    object = clus_tcfs(out_md, thresholds = c(1, 4)))
})

test_that("clus_tcfs does not work if not input a pair of thresholds", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_error(
    object = clus_tcfs(out_md, newdata = data.frame(X1 = 0.5)))
})
