# tests/testthat/test-diff-contrast-length.R
# -----------------------------------------------------------------------------
# Regression test for the contrast-length bug in the differential expression
# API (model$diffO / $diffQ / $diffExp).
#
# Before the fix:
#   - R guards (R/gedi_projections.R::compute_diffQ / compute_diffO) and the
#     C++ guards (src/gedi_differential.cpp::getDiff{O,Q,Exp}_cpp) all enforced
#     length(contrast) == L.
#   - The C++ formula computes Ro %*% H_rotation %*% contrast, where
#     H_rotation is L x num_covariates. So contrast must be num_covariates-long.
#   - Whenever L < num_covariates (rank-deficient H), the only call path the
#     guards permitted triggered an Eigen out-of-bounds read (silent UB).
#
# After the fix:
#   - Both guards check length(contrast) == ncol(H.rotation) (= num_covariates).
#   - A length-num_covariates contrast is accepted; a length-L contrast is
#     rejected with a clear error (when L < num_covariates).
#   - The returned values match a pure-R replica of the documented formula.
# -----------------------------------------------------------------------------

# Fixture: 4 samples, 3 covariates that are linearly dependent
#   cov3 = cov1 + cov2  =>  rank(H^T) = 2  =>  L = 2 < num_covariates = 3
make_diff_fixture <- function() {
  set.seed(1234)
  n_genes   <- 60
  n_cells   <- 80
  n_samples <- 4

  M <- Matrix::Matrix(
    matrix(stats::rpois(n_genes * n_cells, lambda = 5), n_genes, n_cells),
    sparse = TRUE
  )
  rownames(M) <- paste0("Gene_", seq_len(n_genes))
  colnames(M) <- paste0("Cell_", seq_len(n_cells))
  samples <- factor(rep(paste0("S", seq_len(n_samples)), each = n_cells / n_samples))

  H <- matrix(
    c(1, 1, 0, 0,
      0, 0, 1, 1,
      1, 1, 1, 1),
    nrow = 3, byrow = TRUE
  )
  rownames(H) <- c("batch", "condition", "combined")
  colnames(H) <- paste0("S", seq_len(n_samples))

  model <- CreateGEDIObject(
    Samples = samples, M = M, K = 3, H = H, verbose = 0
  )
  model$train(iterations = 10, track_interval = 5)

  list(model = model, H = H, num_cov = nrow(H))
}

# -----------------------------------------------------------------------------
test_that("diffO/diffQ accept length-num_covariates contrast (rank-deficient H)", {
  fx <- make_diff_fixture()
  m  <- fx$model

  # Sanity: fixture really has L < num_covariates
  expect_lt(m$aux$L, fx$num_cov)
  expect_equal(ncol(m$priors$H.rotation), fx$num_cov)

  # length-num_covariates contrast must be accepted by all three methods
  contrast_full <- c(batch = 1, condition = -1, combined = 0)
  expect_no_error(diffO_v   <- m$diffO(contrast_full))
  expect_no_error(diffQ_m   <- m$diffQ(contrast_full))
  expect_no_error(diffExp_m <- m$diffExp(contrast_full))

  expect_length(diffO_v, length(m$metadata$geneIDs))
  expect_equal(dim(diffQ_m),   c(length(m$metadata$geneIDs),
                                 length(m$metadata$cellIDs)))
  expect_equal(dim(diffExp_m), c(length(m$metadata$geneIDs),
                                 length(m$metadata$cellIDs)))
})

# -----------------------------------------------------------------------------
test_that("diffO/diffQ reject length-L contrast when L < num_covariates", {
  fx <- make_diff_fixture()
  m  <- fx$model

  # length-L contrast is the OLD (buggy) contract; must now be rejected with
  # a message that mentions the correct expected length.
  short <- rep(1, m$aux$L)
  expect_error(m$diffO(short),   regexp = as.character(fx$num_cov))
  expect_error(m$diffQ(short),   regexp = as.character(fx$num_cov))
  expect_error(m$diffExp(short), regexp = as.character(fx$num_cov))
})

# -----------------------------------------------------------------------------
test_that("diffO matches the documented formula Ro %*% H.rotation %*% contrast", {
  fx <- make_diff_fixture()
  m  <- fx$model

  contrast <- c(0.5, -0.3, 0.7)
  Ro       <- m$params$Ro            # J x L
  H_rot    <- m$priors$H.rotation    # L x num_cov

  expected <- as.vector(Ro %*% H_rot %*% contrast)
  observed <- as.vector(m$diffO(contrast))

  expect_equal(observed, expected, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
test_that("diffQ matches sum_k Rk %*% H.rotation %*% contrast for each gene-factor pair", {
  fx <- make_diff_fixture()
  m  <- fx$model

  contrast <- c(0.5, -0.3, 0.7)
  H_rot    <- m$priors$H.rotation    # L x num_cov
  H_c      <- as.vector(H_rot %*% contrast)   # length L
  Rk_list  <- m$params$Rk
  D        <- m$params$D
  Bi_list  <- m$params$Bi
  K        <- m$aux$K
  J        <- length(m$metadata$geneIDs)
  N        <- length(m$metadata$cellIDs)

  # Rebuild diffExp = sum_k (Rk %*% H_c) %o% (D[k] * B[k, :])
  B_full <- do.call(cbind, Bi_list)              # K x N
  DB     <- diag(D, nrow = K) %*% B_full         # K x N

  expected <- matrix(0, nrow = J, ncol = N)
  for (k in seq_len(K)) {
    eff_k    <- Rk_list[[k]] %*% H_c             # J x 1
    expected <- expected + eff_k %*% DB[k, , drop = FALSE]
  }

  observed <- m$diffExp(contrast)
  expect_equal(unname(observed), unname(expected), tolerance = 1e-10)
})
