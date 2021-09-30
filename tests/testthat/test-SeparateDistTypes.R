library(testthat)

context("SeparateDistTypes")

test_that("diagonal is -1", {
  expect_true(all(diag(SeparateDistTypes(sample(1:2, 10, replace = TRUE))) == -1))
})
