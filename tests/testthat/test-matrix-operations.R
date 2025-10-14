# Tests for matrix operation helper functions
# These functions are internal but critical for GEDI functionality

test_that("compute_s_0_dense calculates cell scaling factors correctly", {
  skip_if_not_installed("Matrix")

  # Create test data
  J <- 10
  N <- 5
  J_vec <- rep(1, J)
  Y <- matrix(rnorm(J * N, mean = 5, sd = 2), nrow = J, ncol = N)

  # Compute scaling factors
  s_0 <- compute_s_0_dense(J_vec, Y, J)

  # Check dimensions
  expect_length(s_0, N)

  # Check values are finite and positive
  expect_true(all(is.finite(s_0)))
  expect_true(all(s_0 > 0))

  # Check against manual calculation
  expected_s_0 <- colMeans(Y)
  expect_equal(s_0, expected_s_0, tolerance = 1e-10)
})

test_that("compute_s_0 calculates cell scaling factors from sparse matrices", {
  skip_if_not_installed("Matrix")

  # Create sparse test data
  set.seed(123)
  J <- 20
  N <- 10
  M <- Matrix::Matrix(
    rpois(J * N, lambda = 2),
    nrow = J,
    ncol = N,
    sparse = TRUE
  )

  J_vec <- rep(1, J)

  # Compute scaling factors
  s_0 <- compute_s_0(J_vec, M, J)

  # Check dimensions
  expect_length(s_0, N)

  # Check values are finite and positive
  expect_true(all(is.finite(s_0)))
  expect_true(all(s_0 > 0))

  # Check against manual calculation (with constant)
  expected_s_0 <- (Matrix::colSums(M) + 0.01) / J
  expect_equal(s_0, as.vector(expected_s_0), tolerance = 1e-10)
})

test_that("compute_o_0_dense calculates gene offsets correctly", {
  skip_if_not_installed("Matrix")

  # Create test data
  J <- 15
  N <- 8
  Yp <- matrix(rnorm(J * N, mean = 0, sd = 1), nrow = J, ncol = N)
  N_vec <- rep(1, N)

  # Compute offsets
  o_0 <- compute_o_0_dense(Yp, N_vec, N)

  # Check dimensions
  expect_length(o_0, J)

  # Check values are finite
  expect_true(all(is.finite(o_0)))

  # Check against manual calculation
  expected_o_0 <- rowMeans(Yp)
  expect_equal(o_0, expected_o_0, tolerance = 1e-10)
})

test_that("compute_o_0 calculates gene offsets from sparse matrices", {
  skip_if_not_installed("Matrix")

  # Create sparse test data
  set.seed(456)
  J <- 25
  N <- 12
  Mp <- Matrix::Matrix(
    rnorm(J * N, mean = 1, sd = 0.5),
    nrow = J,
    ncol = N,
    sparse = TRUE
  )

  N_vec <- rep(1, N)

  # Compute offsets
  o_0 <- compute_o_0(Mp, N_vec, N)

  # Check dimensions
  expect_length(o_0, J)

  # Check values are finite
  expect_true(all(is.finite(o_0)))

  # Check against manual calculation (with constant)
  expected_o_0 <- (Matrix::rowSums(Mp) + 1e-5) / N
  expect_equal(o_0, as.vector(expected_o_0), tolerance = 1e-10)
})

test_that("compute_Yp_dense normalizes expression correctly", {
  skip_if_not_installed("Matrix")

  # Create test data
  J <- 10
  N <- 5
  Y <- matrix(rnorm(J * N, mean = 5, sd = 2), nrow = J, ncol = N)
  J_vec <- rep(1, J)
  s_0 <- colMeans(Y)

  # Normalize
  Yp <- compute_Yp_dense(Y, J_vec, s_0)

  # Check dimensions
  expect_equal(dim(Yp), c(J, N))

  # Check values are finite
  expect_true(all(is.finite(Yp)))

  # Check normalization: should have mean close to 0 for each cell
  cell_means <- colMeans(Yp)
  expect_true(all(abs(cell_means) < 1e-10))
})

test_that("compute_Mp normalizes sparse count matrix correctly", {
  skip_if_not_installed("Matrix")

  # Create sparse test data
  set.seed(789)
  J <- 20
  N <- 10
  M <- Matrix::Matrix(
    rpois(J * N, lambda = 5),
    nrow = J,
    ncol = N,
    sparse = TRUE
  )

  J_vec <- rep(1, J)
  s_0 <- (Matrix::colSums(M) + 0.01) / J

  # Normalize
  Mp <- compute_Mp(M, J_vec, s_0)

  # Check dimensions
  expect_equal(dim(Mp), c(J, N))

  # Check it's still sparse
  expect_s4_class(Mp, "sparseMatrix")

  # Check values are finite where non-zero
  expect_true(all(is.finite(Mp@x)))
})

test_that("VecVecProduct computes outer product correctly", {
  # Create test vectors
  a <- c(1, 2, 3)
  b <- c(4, 5)

  # Compute outer product
  result <- VecVecProduct(a, b)

  # Check dimensions
  expect_equal(dim(result), c(3, 2))

  # Check against R's outer product
  expected <- outer(a, b)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("MatVecProduct computes matrix-vector product correctly", {
  # Create test data
  A <- matrix(1:12, nrow = 3, ncol = 4)
  b <- c(1, 2, 3, 4)

  # Compute product
  result <- MatVecProduct(A, b)

  # Check dimensions
  expect_length(result, 3)

  # Check against R's matrix multiplication
  expected <- as.vector(A %*% b)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("eigenSparseMatVecProduct computes sparse matrix-vector product", {
  skip_if_not_installed("Matrix")

  # Create sparse matrix
  set.seed(111)
  A <- Matrix::Matrix(
    sample(0:5, 50, replace = TRUE),
    nrow = 10,
    ncol = 5,
    sparse = TRUE
  )
  b <- rnorm(5)

  # Compute product
  result <- eigenSparseMatVecProduct(A, b)

  # Check dimensions
  expect_length(result, 10)

  # Check values are finite
  expect_true(all(is.finite(result)))

  # Check against R's matrix multiplication
  expected <- as.vector(A %*% b)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("Matrix operations handle edge cases properly", {
  skip_if_not_installed("Matrix")

  # Test with very small matrix
  J <- 2
  N <- 2
  M <- Matrix::Matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, sparse = TRUE)
  J_vec <- rep(1, J)

  expect_silent({
    s_0 <- compute_s_0(J_vec, M, J)
    expect_length(s_0, N)
    expect_true(all(s_0 > 0))
  })

  # Test with zero matrix
  M_zero <- Matrix::Matrix(0, nrow = 5, ncol = 5, sparse = TRUE)
  J_vec_zero <- rep(1, 5)

  expect_silent({
    s_0_zero <- compute_s_0(J_vec_zero, M_zero, 5)
    # Should still be positive due to added constant
    expect_true(all(s_0_zero > 0))
  })
})
