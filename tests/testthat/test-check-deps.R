# tests/testthat/test-check-deps.R
# -----------------------------------------------------------
# Tests for check_optional_dependencies(). No model needed.
# -----------------------------------------------------------

test_that("check_optional_dependencies returns a named logical vector", {
  result <- check_optional_dependencies()
  expect_type(result, "logical")
  expect_true(length(result) > 0)
  expect_true(all(c("hdf5r", "uwot", "digest") %in% names(result)))
})
