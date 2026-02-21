# tests/testthat/test-CreateGEDIObject.R
# -----------------------------------------------------------
# Tests for CreateGEDIObject(): object creation, validation,
# training, and basic output structure.
# -----------------------------------------------------------

# ---- shared fixture ----
dat <- make_synthetic_data(n_genes = 200, n_cells = 400,
                           n_samples = 4, seed = 42)

# ===========================================================
# 1. Object creation — correct class and metadata
# ===========================================================
test_that("CreateGEDIObject returns an R6 object with correct class", {
  model <- CreateGEDIObject(
    Samples = dat$sample_labels,
    M       = dat$M,
    K       = 10,
    mode    = "Bsphere",
    verbose = 0
  )
  expect_s3_class(model, "R6")
  # print method should work without error
  expect_output(print(model))
})

test_that("CreateGEDIObject preserves gene and cell IDs before training", {
  model <- CreateGEDIObject(
    Samples = dat$sample_labels,
    M       = dat$M,
    K       = 10,
    mode    = "Bsphere",
    verbose = 0
  )
  meta <- model$metadata
  expect_equal(meta$geneIDs, rownames(dat$M))
  expect_equal(meta$cellIDs, colnames(dat$M))
  expect_equal(levels(meta$sampleIDs), levels(dat$sample_labels))
})

# ===========================================================
# 2. Input validation — bad inputs should error
# ===========================================================
test_that("CreateGEDIObject errors on mismatched sample labels length", {
  bad_labels <- factor(rep("S1", 10))
  expect_error(
    CreateGEDIObject(Samples = bad_labels, M = dat$M, K = 10, verbose = 0)
  )
})

test_that("CreateGEDIObject errors when K >= min(n_genes, n_cells)", {
  expect_error(
    CreateGEDIObject(Samples = dat$sample_labels, M = dat$M, K = 500,
                     verbose = 0)
  )
})

test_that("CreateGEDIObject errors on missing data matrix", {
  expect_error(
    CreateGEDIObject(Samples = dat$sample_labels, K = 10, verbose = 0)
  )
})

# ===========================================================
# 3. Mode selection
# ===========================================================
test_that("Both Bsphere and Bl2 modes create models without error", {
  expect_no_error(
    CreateGEDIObject(Samples = dat$sample_labels, M = dat$M,
                     K = 5, mode = "Bsphere", verbose = 0)
  )
  expect_no_error(
    CreateGEDIObject(Samples = dat$sample_labels, M = dat$M,
                     K = 5, mode = "Bl2", verbose = 0)
  )
})

# ===========================================================
# 4. Training runs and produces expected output
# ===========================================================
test_that("Model trains without error and populates parameters", {
  model <- CreateGEDIObject(
    Samples = dat$sample_labels,
    M       = dat$M,
    K       = 5,
    mode    = "Bsphere",
    verbose = 0
  )
  expect_no_error(model$train(iterations = 10))

  # Core parameter matrices should exist with correct shapes
  Z <- model$params$Z
  expect_true(is.matrix(Z))
  expect_equal(nrow(Z), 200)
  expect_equal(ncol(Z), 5)

  D <- model$params$D
  expect_true(is.numeric(D))
  expect_length(D, 5)
})

test_that("Metadata dimensions are correct after training", {
  model <- CreateGEDIObject(
    Samples = dat$sample_labels,
    M       = dat$M,
    K       = 10,
    mode    = "Bsphere",
    verbose = 0
  )
  model$train(iterations = 10)

  meta <- model$metadata
  expect_equal(meta$n_genes,   200)
  expect_equal(meta$n_cells,   400)
  expect_equal(meta$n_samples, 4)

  expect_equal(model$aux$K,    10)
  expect_equal(model$aux$mode, "Bsphere")
})

test_that("Training populates convergence tracking", {
  model <- CreateGEDIObject(
    Samples = dat$sample_labels,
    M       = dat$M,
    K       = 5,
    mode    = "Bsphere",
    verbose = 0
  )
  model$train(iterations = 20, track_interval = 5)

  tracking <- model$tracking
  expect_true(is.list(tracking))
  expect_true(length(tracking) > 0)
  # Each tracked vector should have entries
  expect_true(length(tracking[[1]]) > 0)
})

# ===========================================================
# 5. Projections (lazy accessors)
# ===========================================================
test_that("Projections ZDB and DB have correct dimensions", {
  model <- make_trained_model(K = 5, iterations = 10)

  zdb <- model$projections$ZDB
  expect_true(is.matrix(zdb) || inherits(zdb, "Matrix"))
  expect_equal(nrow(zdb), 200)  # genes
  expect_equal(ncol(zdb), 400)  # cells

  db <- model$projections$DB
  expect_equal(nrow(db), 5)     # K
  expect_equal(ncol(db), 400)   # cells
})

# ===========================================================
# 6. Imputation accessor
# ===========================================================
test_that("Imputed Y matrix has correct dimensions", {
  model <- make_trained_model(K = 5, iterations = 10)

  Y_imp <- model$imputed$Y()
  expect_true(is.matrix(Y_imp) || inherits(Y_imp, "Matrix"))
  expect_equal(nrow(Y_imp), 200)
  expect_equal(ncol(Y_imp), 400)
})

# ===========================================================
# 7. Reproducibility with same seed
# ===========================================================
test_that("Same seed produces identical synthetic data", {
  dat1 <- make_synthetic_data(seed = 123)
  dat2 <- make_synthetic_data(seed = 123)

  expect_identical(as.matrix(dat1$M), as.matrix(dat2$M))
  expect_identical(dat1$sample_labels, dat2$sample_labels)
})

# ===========================================================
# 8. Original user scenario: 2000 cells, 1000 genes, K=10
# ===========================================================
test_that("Full scenario: NB(1000x2000), 10 samples, K=10, Bsphere", {
  big <- make_synthetic_data(n_genes = 1000, n_cells = 2000,
                             n_samples = 10, seed = 99)

  model <- CreateGEDIObject(
    Samples = big$sample_labels,
    M       = big$M,
    K       = 10,
    mode    = "Bsphere",
    verbose = 0
  )

  # Train for a few iterations
  model$train(iterations = 15)

  # Verify dimensions after training
  expect_equal(model$metadata$n_genes,   1000)
  expect_equal(model$metadata$n_cells,   2000)
  expect_equal(model$metadata$n_samples, 10)
  expect_equal(model$aux$K,              10)
  expect_equal(model$aux$mode,           "Bsphere")

  Z <- model$params$Z
  expect_equal(dim(Z), c(1000, 10))

  # Projections should be computable
  zdb <- model$projections$ZDB
  expect_equal(dim(zdb), c(1000, 2000))
})