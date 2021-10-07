library(testthat)

context('plotMDS_andLegend')


test_that("data and selected scores match", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  batch1  = sample(1:2, size = 10, replace = TRUE)
  expect_error(plotMDS_andLegend(data1, batch1, SelectScores = rep(1, 200)), 'data and SelectScores do not match each other')
})
