context("QC_bootstrap")

test_that("error message correct nrows", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  data2 = matrix(rnorm(10*200,mean=0,sd=1), ncol=10, nrow=200)
  batch  = sample(1:2, size = 10, replace = TRUE)
  groups = sample(1:2, size = 10, replace = TRUE)
  expect_error(QC_bootstrap(list(data1, data2), biol.groups = groups, batches = batch), "Error: datasets should have equal number of rows")
})

test_that("error message NA values", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  data1[1,1] <- NA
  batch  = sample(1:2, size = 10, replace = TRUE)
  groups = sample(1:2, size = 10, replace = TRUE)
  expect_error(QC_bootstrap(list(data1), biol.groups = groups, batches = batch), "Error: datasets should not have NA values")
})

test_that("error message finite values", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  data1[1,1] <- Inf
  batch  = sample(1:2, size = 10, replace = TRUE)
  groups = sample(1:2, size = 10, replace = TRUE)
  expect_error(QC_bootstrap(list(data1), biol.groups = groups, batches = batch), "Error: datasets should contain only finite values")
})

test_that("datalist and method names length matches", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  batch  = sample(1:2, size = 10, replace = TRUE)
  groups = sample(1:2, size = 10, replace = TRUE)
  expect_error(QC_bootstrap(list(data1), biol.groups = groups, batches = batch, method_names = c("1", "2")), "The list of method names has to be of the same length as the data list")
})


