context("ReAdjustData")

test_that("inputs match in size",{
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  expect_error(ReAdjustData(data1, rep(1, 100), rep(1, 200)), 'Inputs for RescaleData do not match in size')
})
