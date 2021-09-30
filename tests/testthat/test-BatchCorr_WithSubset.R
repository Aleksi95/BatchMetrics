context("BatchCorr_WithSubset")

test_that("harman requires sample types", {
  data1 = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
  batch1  = sample(1:2, size = 10, replace = TRUE)
  subset  = sample(1:10, size = 5, replace = FALSE)
  expect_error(BatchCorr_WithSubset(data1, batch1, SubsetIndex = subset, CorrMethod = "harman"), "Error: harman requires sample types")
})

