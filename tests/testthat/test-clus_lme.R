test_that("clus_mle work with one contiuous covariate", {
  data(data_3class)
  expect_s3_class(clus_lme(fixed_formula = Y ~ X1, name_class = "D",
                           name_clust = "id_Clus", data = data_3class),
                  "clus_lme")
})


