context("evalBatchEffect")

test_that("evalBatchEffect returns a named numeric vector with bio, batch, ratio", {
  set.seed(42)
  Data = matrix(rnorm(500 * 30, mean = 0, sd = 1), nrow = 500, ncol = 30)
  sample_types = sample(c("A", "B"), 30, replace = TRUE)
  batch = sample(c("batch1", "batch2"), 30, replace = TRUE)

  result = evalBatchEffect(Data, sample_types, batch, metric = "F-score")

  expect_true(is.numeric(result))
  expect_length(result, 3)
  expect_named(result, c("bio", "batch", "ratio"))
  expect_gte(result["bio"], 0)
  expect_gte(result["batch"], 0)
})

test_that("runMetrics returns a list with results for each metric", {
  set.seed(42)
  Data = matrix(rnorm(500 * 30, mean = 0, sd = 1), nrow = 500, ncol = 30)
  sample_types = sample(c("A", "B"), 30, replace = TRUE)
  batch = sample(c("batch1", "batch2"), 30, replace = TRUE)

  results = runMetrics(Data, sample_types, batch,
                       metrics = c("F-score", "Davies-Bouldin"))

  expect_type(results, "list")
  expect_length(results, 2)
  expect_named(results, c("F-score", "Davies-Bouldin"))
  for (r in results) {
    expect_true(is.numeric(r))
    expect_length(r, 3)
    expect_named(r, c("bio", "batch", "ratio"))
  }
})
