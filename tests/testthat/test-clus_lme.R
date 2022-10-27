test_that("clus_mle works with one contiuous covariate", {
  data(data_3class)
  expect_s3_class(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      name_clust = "id_Clus", data = data_3class),
    class = "clus_lme")
})

test_that("clus_mle works with two contiuous covariates", {
  data(data_3class)
  expect_s3_class(
    object = clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
                      name_clust = "id_Clus", data = data_3class),
    class = "clus_lme")
})

test_that("clus_mle works with one covariate and subset", {
  data(data_3class)
  expect_s3_class(
    object = clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
                      name_clust = "id_Clus", data = data_3class,
                      subset = id_Clus %in% sample(1:30, 20, replace = FALSE)),
    class = "clus_lme")
})

test_that("clus_mle works with one covariate and levl_class", {
  data(data_3class)
  expect_s3_class(
    object = clus_lme(fixed_formula = Y ~ X1 + X2, name_class = "D",
                      name_clust = "id_Clus", data = data_3class,
                      levl_class = c("1", "3", "2")),
    class = "clus_lme")
})


test_that("clus_mle works with one covariate and boxcox transformation", {
  data(data_3class_bcx)
  expect_s3_class(
    object = clus_lme(fixed_formula = Y ~ X, name_class = "D",
                      name_clust = "id_Clus", data = data_3class_bcx,
                      boxcox = TRUE),
    class = "clus_lme")
})


test_that("clus_mle does not work without input fixed_formula", {
  data(data_3class)
  expect_error(
    object = clus_lme(name_class = "D", name_clust = "id_Clus",
                      data = data_3class),
    regexp = "agrument \"fixed_formula\" must be a formula of the form \"resp ~ pred\"")
})

test_that("clus_mle does not work if input fixed_formula as dot in predictors", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ . , name_class = "D",
                      name_clust = "id_Clus", data = data_3class))
})

test_that("clus_mle does not work if input fixed_formula as dot in response", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = . ~ X1 , name_class = "D",
                      name_clust = "id_Clus", data = data_3class))
})

test_that("clus_mle does not work if input fixed_formula as dot", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = . ~ . , name_class = "D",
                      name_clust = "id_Clus", data = data_3class))
})

test_that("clus_mle does not work if input wrong fixed_formula", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y, name_class = "D",
                      name_clust = "id_Clus", data = data_3class))
})

test_that("clus_mle does not work if input wrong fixed_formula", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X, name_class = "D",
                      name_clust = "id_Clus", data = data_3class))
})

test_that("clus_mle does not work if input wrong fixed_formula", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Z ~ X1, name_class = "D",
                      name_clust = "id_Clus", data = data_3class))
})


test_that("clus_mle does not work if input missing name_class", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1,
                      name_clust = "id_Clus", data = data_3class),
    regexp = "agrument \"name_class\" was either missing or wrong name!")
})

test_that("clus_mle does not work if input missing name_class", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "",
                      name_clust = "id_Clus", data = data_3class),
    regexp = "agrument \"name_class\" was either missing or wrong name!")
})

test_that("clus_mle does not work if name_class is not character", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = D,
                      name_clust = "id_Clus", data = data_3class),
    regexp = "agrument \"name_class\" was either missing or wrong name!")
})

test_that("clus_mle does not work if input wrong name_class", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "X2",
                      name_clust = "id_Clus", data = data_3class),
    regexp = "agrument \"name_class\" must have 3 levels or classes!")
})

test_that("clus_mle does not work if input missing name_clust", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D", name_clust = "",
                      data = data_3class),
    regexp = "agrument \"name_clust\" was either missing or wrong name!")
})

test_that("clus_mle does not work if input missing name_clust", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      data = data_3class),
    regexp = "agrument \"name_clust\" was either missing or wrong name!")
})

test_that("clus_mle does not work if input missing name_clust", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      name_clust = "D", data = data_3class),
    regexp = "agrument \"name_clust\" cannot be the name of neither test, covariates nor classes!")
})

test_that("clus_mle does not work if input wrong levl_class", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      name_clust = "id_Clus", data = data_3class,
                      levl_class = 3))
})

test_that("clus_mle does not work if input wrong levl_class", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      name_clust = "id_Clus", data = data_3class,
                      levl_class = c("1", "2", NA)))
})

test_that("clus_mle does not work if input wrong levl_class", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      name_clust = "id_Clus", data = data_3class,
                      levl_class = c(1, 2, NA)))
})


test_that("clus_mle does not works for boxcox transformation in case of negative values", {
  data(data_3class)
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      name_clust = "id_Clus", data = data_3class,
                      boxcox = TRUE),
    regexp = "Cannot apply Box-Cox transform for negative values.")
})

test_that("clus_mle does not works without action for missing data", {
  data(data_3class)
  data_3class$X1[sample(1:30, 6, replace = FALSE)] <- NA
  expect_error(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      name_clust = "id_Clus", data = data_3class))
})

test_that("clus_mle does works with na.omit() for missing data", {
  data(data_3class)
  data_3class$X1[sample(1:30, 6, replace = FALSE)] <- NA
  expect_s3_class(
    object = clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                      name_clust = "id_Clus", data = data_3class,
                      na_action = na.omit),
    class = "clus_lme")
})

