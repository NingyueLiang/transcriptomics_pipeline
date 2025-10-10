library(testthat)
library(transcriptomicsPipeline)

test_that("setup_config is properly configured", {
  expect_true("cran_packages" %in% names(setup_config))
  expect_true("bioc_packages" %in% names(setup_config))
  expect_true(is.character(setup_config$cran_packages))
  expect_true(is.character(setup_config$bioc_packages))
})

test_that("install_single_package validates inputs", {
  expect_error(install_single_package(NULL), "package_name cannot be NULL")
  expect_error(install_single_package("", "invalid_source"), "package_name must be a non-empty character string")
})

test_that("setup_packages validates config", {
  expect_error(setup_packages(NULL))
  expect_error(setup_packages(list()))
  expect_error(setup_packages(list(cran_packages = "test")))
})