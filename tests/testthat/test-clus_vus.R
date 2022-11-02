test_that("clus_vus works with one contiuous covariate", {
  data(data_3class)
  out1 <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                   name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_vus(out1, newdata = data.frame(X1 = 0.5)),
    class = "clus_vus")
})

test_that("clus_vus works in case of no covariate", {
  data(data_3class)
  out1 <- clus_lme(fixed_formula = Y ~ 1, name_class = "D",
                   name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_vus(out1, newdata = data.frame(X1 = 0.5)),
    class = "clus_vus")
})

test_that("clus_vus works in case of no covariate", {
  data(data_3class)
  out1 <- clus_lme(fixed_formula = Y ~ 1, name_class = "D",
                   name_clust = "id_Clus", data = data_3class)
  expect_s3_class(
    object = clus_vus(out1),
    class = "clus_vus")
})

test_that("clus_vus does not work in case of one contiuous covariate without newdata", {
  data(data_3class)
  out1 <- clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                   name_clust = "id_Clus", data = data_3class)
  expect_error(
    object = clus_vus(out1))
})
