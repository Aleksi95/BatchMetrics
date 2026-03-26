context("ADS_A and diluteSeries2")

test_that("ADS_A returns a matrix of the same dimensions", {
  set.seed(42)
  Data = matrix(rnorm(200 * 20, mean = 0, sd = 1), nrow = 200, ncol = 20)
  samples = sample(c("A", "B"), 20, replace = TRUE)

  result = ADS_A(Data, samples, p = 0.5)

  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(Data))
})

test_that("ADS_A with p=0 returns unchanged data", {
  set.seed(42)
  Data = matrix(rnorm(200 * 20, mean = 0, sd = 1), nrow = 200, ncol = 20)
  samples = sample(c("A", "B"), 20, replace = TRUE)

  result = ADS_A(Data, samples, p = 0)

  expect_equal(result, Data)
})

test_that("diluteSeries2 returns a list of the correct length", {
  set.seed(42)
  Data = matrix(rnorm(200 * 20, mean = 0, sd = 1), nrow = 200, ncol = 20)
  samples = sample(c("A", "B"), 20, replace = TRUE)
  perc = seq(0, 1, by = 0.25)

  series = diluteSeries2(Data, samples, perc = perc)

  expect_type(series, "list")
  expect_length(series, length(perc))
})

test_that("diluteSeries2 returns matrices of correct dimensions", {
  set.seed(42)
  Data = matrix(rnorm(200 * 20, mean = 0, sd = 1), nrow = 200, ncol = 20)
  samples = sample(c("A", "B"), 20, replace = TRUE)
  perc = c(0, 0.5, 1)

  series = diluteSeries2(Data, samples, perc = perc)

  for (d in series) {
    expect_equal(dim(d), dim(Data))
  }
})
