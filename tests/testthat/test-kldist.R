context("kldist")

test_that("kldist returns a single numeric value", {
  set.seed(42)
  # kldist expects samples x features
  Data = matrix(rnorm(20 * 50, mean = 0, sd = 1), nrow = 20, ncol = 50)
  batch = factor(sample(c("batch1", "batch2"), 20, replace = TRUE))

  kl = kldist(Data, batch)

  expect_true(is.numeric(kl))
  expect_length(kl, 1)
  expect_gte(kl, 0)
})

test_that("kldistTwo returns a single numeric value", {
  set.seed(42)
  xb1 = matrix(rnorm(10 * 20, mean = 0, sd = 1), nrow = 10, ncol = 20)
  xb2 = matrix(rnorm(10 * 20, mean = 1, sd = 1), nrow = 10, ncol = 20)

  kl_two = kldistTwo(xb1, xb2)

  expect_true(is.numeric(kl_two))
  expect_length(kl_two, 1)
  expect_gte(kl_two, 0)
})
