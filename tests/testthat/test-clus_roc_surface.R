test_that("clus_roc_surface works without covariate", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ 1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  clus_roc_surface(out_md)
  expect_known_scene("plot_roc_surface_0_cov")
})

test_that("clus_roc_surface works with one covariate", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  clus_roc_surface(out_md, newdata = data.frame(X1 = 0.5))
  expect_known_scene("plot_roc_surface_1_cov")
})

test_that("clus_roc_surface works with two covariates", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  clus_roc_surface(out_md, newdata = data.frame(X1 = 0.5, X2 = 0))
  expect_known_scene("plot_roc_surface_2_cov")
})

test_that("clus_roc_surface works with Box-Cox transformation", {
  data("data_3class_bcx")
  out_md <- clus_lme(fixed_formula = Y ~ X, name_class = "D",
                     name_clust = "id_Clus", data = data_3class_bcx,
                     boxcox = TRUE)
  clus_roc_surface(out_md, newdata = data.frame(X = 1))
  expect_known_scene("plot_roc_surface_bcx")
})

test_that("clus_roc_surface only works with one point of covariate", {
  data(data_3class)
  out_md <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                     name_clust = "id_Clus", data = data_3class)
  expect_error(
    clus_roc_surface(out_md, newdata = data.frame(X1 = c(0.5, 1)))
  )
})
