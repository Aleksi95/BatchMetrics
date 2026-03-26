context("Ftest2")

test_that("Ftest2 with AOV returns a single numeric F-statistic", {
  set.seed(42)
  # rows are groups, columns are observations
  data = matrix(rnorm(11 * 20, mean = 0, sd = 1), nrow = 11, ncol = 20)

  result = Ftest2(data, test = "AOV")

  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_gte(result, 0)
})

test_that("Ftest2 with JT returns a named vector with JT and pvalue", {
  set.seed(42)
  data = matrix(rnorm(11 * 20, mean = 0, sd = 1), nrow = 11, ncol = 20)

  result = Ftest2(data, test = "JT", alternative = "increasing")

  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true("JT" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_gte(result["pvalue"], 0)
  expect_lte(result["pvalue"], 1)
})
