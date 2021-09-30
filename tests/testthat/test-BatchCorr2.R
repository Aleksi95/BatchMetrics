context("BatchCorr2")

test_that("harman requires sample types", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  batch1  = sample(1:2, size = 10, replace = TRUE)
  batch2  = sample(1:2, size = 10, replace = TRUE)
  expect_error(BatchCorr2(data1, batch1, batch2, method = "harman"), "Error: harman requires sample types")
})

test_that("samplecenter requires sample types", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  batch1  = sample(1:2, size = 10, replace = TRUE)
  batch2  = sample(1:2, size = 10, replace = TRUE)
  expect_error(BatchCorr2(data1, batch1, batch2, method = "samplescenter"), "Error: samplescenter requires sample types")
})

test_that("subset correction only defined on combat and harman", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  batch1  = sample(1:2, size = 10, replace = TRUE)
  batch2  = sample(1:2, size = 10, replace = TRUE)
  subset1  = sample(1:10, size = 5, replace = FALSE)
  subset2 = sample(1:10, size = 5, replace = FALSE)
  expect_error(BatchCorr2(data1, batch1, batch2, subset1 = subset1, subset2 = subset2, method = "clustercenter"), "Error: subset correction defined only for combat and harman")
})
