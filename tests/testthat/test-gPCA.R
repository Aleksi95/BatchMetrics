context("gPCA")

test_that("gPCA.batchdetect returns correct structure", {
  set.seed(42)
  # gPCA.batchdetect expects samples x features
  Data = matrix(rnorm(20 * 100, mean = 0, sd = 1), nrow = 20, ncol = 100)
  batch = sample(1:2, 20, replace = TRUE)

  result = gPCA.batchdetect(Data, batch, nperm = 20, seed = 1)

  expect_type(result, "list")
  expect_true("delta" %in% names(result))
  expect_true("p.val" %in% names(result))
  expect_true("varPCu1" %in% names(result))
  expect_true("varPCg1" %in% names(result))
  expect_true(is.numeric(result$delta))
  expect_gte(result$delta, 0)
})

test_that("gPCA_percentage returns a numeric value", {
  set.seed(42)
  # gPCA_percentage expects features x samples
  Data = matrix(rnorm(100 * 20, mean = 0, sd = 1), nrow = 100, ncol = 20)
  sample_types = sample(c("A", "B"), 20, replace = TRUE)

  pct = gPCA_percentage(Data, sample_types, nperm = 20)

  expect_true(is.numeric(pct))
  expect_length(pct, 1)
})
