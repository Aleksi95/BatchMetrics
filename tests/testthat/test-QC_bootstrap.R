context("QC_bootstrap")

test_that("DataNorm is first parameter (not data_list)", {
  expect_true("DataNorm" %in% names(formals(QC_bootstrap)))
  expect_false("data_list" %in% names(formals(QC_bootstrap)))
})

test_that("function signature contains expected parameters", {
  params = names(formals(QC_bootstrap))
  expect_true("biol.groups" %in% params)
  expect_true("batches" %in% params)
  expect_true("method_names" %in% params)
  expect_true("metrics" %in% params)
  expect_true("iters" %in% params)
  expect_true("corrMethod" %in% params)
  expect_true("dilute_samples" %in% params)
  expect_true("parallel" %in% params)
  expect_true("usePCA" %in% params)
  expect_true("nPCs" %in% params)
})

test_that("function signature no longer contains removed parameters", {
  params = names(formals(QC_bootstrap))
  expect_false("Fscore_method" %in% params)
  expect_false("savefile" %in% params)
  expect_false("filename" %in% params)
  expect_false("plot" %in% params)
})

test_that("default parameter values match documentation", {
  f = formals(QC_bootstrap)
  expect_null(f$method_names)
  expect_equal(as.character(f$dist_method), "pearson")
  expect_equal(as.logical(f$scaledF), FALSE)
  expect_equal(as.integer(f$iters), 50L)
  expect_equal(as.character(f$corrMethod), "ComBat")
  expect_equal(as.logical(f$dilute_samples), FALSE)
  expect_equal(as.logical(f$parallel), TRUE)
  expect_equal(as.integer(f$nCores), 16L)
  expect_equal(as.logical(f$zeroRows), FALSE)
  expect_equal(as.logical(f$usePCA), FALSE)
  expect_equal(as.integer(f$nPCs), 50L)
})


