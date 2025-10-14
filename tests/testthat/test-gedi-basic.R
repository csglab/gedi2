# Basic tests for GEDI package
# Tests package loading, object creation, and basic functionality

test_that("GEDI package loads successfully", {
  expect_true("gedi" %in% loadedNamespaces())
})

test_that("CreateGEDIObject function exists", {
  expect_true(exists("CreateGEDIObject"))
  expect_true(is.function(CreateGEDIObject))
})

test_that("GEDI R6 class exists", {
  expect_true(exists("GEDI"))
  expect_true(R6::is.R6Class(GEDI))
})

test_that("CreateGEDIObject creates valid object with minimal input", {
  skip_if_not_installed("Matrix")

  # Create minimal test data
  set.seed(123)
  n_genes <- 50
  n_cells <- 30
  K <- 5

  # Sparse count matrix
  M <- Matrix::Matrix(
    rpois(n_genes * n_cells, lambda = 2),
    nrow = n_genes,
    ncol = n_cells,
    sparse = TRUE
  )

  # Sample labels
  Samples <- rep(c("Sample1", "Sample2"), each = n_cells / 2)

  # Create GEDI object
  expect_silent({
    model <- CreateGEDIObject(
      Samples = Samples,
      M = M,
      K = K,
      verbose = 0
    )
  })

  # Check object structure
  expect_true(R6::is.R6(model))
  expect_s3_class(model, "GEDI")

  # Check object has expected methods
  expect_true("setup" %in% names(model))
  expect_true("initialize_lvs" %in% names(model))
  expect_true("optimize" %in% names(model))
  expect_true("train" %in% names(model))
})

test_that("GEDI object setup method works", {
  skip_if_not_installed("Matrix")

  # Create minimal test data
  set.seed(456)
  n_genes <- 30
  n_cells <- 20
  K <- 3

  M <- Matrix::Matrix(
    rpois(n_genes * n_cells, lambda = 3),
    nrow = n_genes,
    ncol = n_cells,
    sparse = TRUE
  )

  Samples <- rep(c("A", "B"), each = n_cells / 2)

  model <- CreateGEDIObject(
    Samples = Samples,
    M = M,
    K = K,
    verbose = 0
  )

  # Setup should work without error
  expect_silent({
    model$setup()
  })

  # After setup, model should be in ready state
  expect_true(model$setup_complete)
})

test_that("Input validation works for CreateGEDIObject", {
  skip_if_not_installed("Matrix")

  M <- Matrix::Matrix(rpois(100, 2), 10, 10, sparse = TRUE)
  Samples <- rep("Sample1", 10)

  # Test K validation
  expect_error(
    CreateGEDIObject(Samples = Samples, M = M, K = 0),
    regexp = "K"
  )

  expect_error(
    CreateGEDIObject(Samples = Samples, M = M, K = -1),
    regexp = "K"
  )

  # Test Samples length mismatch
  expect_error(
    CreateGEDIObject(Samples = rep("Sample1", 5), M = M, K = 3),
    regexp = "length|dimension|ncol"
  )
})

test_that("Helper functions are available internally", {
  # These should be internal to the package
  # Just verify the exported functions work
  expect_true(exists("CreateGEDIObject"))
})
