context("GetDistMatrix")


test_that("data contains nas", {
  expect_error(GetDistMatrix(matrix(c(1,2,3,NA), nrow=2)), "Data should not contain NAs!")
})
