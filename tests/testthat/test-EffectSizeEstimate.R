context("EffectSizeEstimate")

test_that("EffectSizeEstimate returns correct structure with MDS", {
  set.seed(42)
  Data = matrix(rnorm(500 * 20, mean = 0, sd = 1), nrow = 500, ncol = 20)
  sample_types = sample(c("A", "B"), 20, replace = TRUE)

  result = EffectSizeEstimate(Data, sample_types, dim.reduct.method = "MDS")

  expect_type(result, "list")
  expect_named(result, c("coef", "pval", "R2"))
  expect_true(is.numeric(result$coef))
  expect_true(is.numeric(result$pval))
  expect_true(is.numeric(result$R2))
  expect_gte(result$coef, 0)
  expect_gte(result$pval, 0)
  expect_lte(result$pval, 1)
  expect_gte(result$R2, 0)
  expect_lte(result$R2, 1)
})

test_that("EffectSizeEstimate returns correct structure with PCA", {
  set.seed(42)
  Data = matrix(rnorm(500 * 20, mean = 0, sd = 1), nrow = 500, ncol = 20)
  sample_types = sample(c("group1", "group2"), 20, replace = TRUE)

  result = EffectSizeEstimate(Data, sample_types, dim.reduct.method = "PCA")

  expect_type(result, "list")
  expect_named(result, c("coef", "pval", "R2"))
  expect_true(is.numeric(result$coef))
  expect_gte(result$coef, 0)
})
