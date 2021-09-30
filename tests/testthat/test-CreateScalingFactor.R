context("CreateScalingFactor")

test_that("inputs match", {
  expect_error(CreateScalingFactor(rep(1,10), rep(1, 20)), 'Error in CreateScalingFactor. Inputs do not match')
})
