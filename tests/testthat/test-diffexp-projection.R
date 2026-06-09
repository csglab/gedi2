# tests/testthat/test-diffexp-projection.R
# -----------------------------------------------------------------------------
# Tests for the "diffexp" and "diffadb" projection types wired into
# model$diffADB() and plot_features().
#
# Fixture: 40 genes, 60 cells, 2 samples, K=3 latent factors.
#   - H design matrix: 1 covariate (condition), 2 samples =>
#       ncol(H.rotation) == 1 (full rank, so L == num_cov == 1)
#   - C pathway matrix: 5 pathways, 0/1 gene membership
# -----------------------------------------------------------------------------

make_diff_proj_fixture <- function() {
  set.seed(777)
  n_genes     <- 40
  n_cells     <- 60
  n_samples   <- 2
  n_pathways  <- 5
  K           <- 3

  M <- Matrix::Matrix(
    matrix(stats::rpois(n_genes * n_cells, lambda = 4), n_genes, n_cells),
    sparse = TRUE
  )
  rownames(M) <- paste0("Gene_", seq_len(n_genes))
  colnames(M) <- paste0("Cell_", seq_len(n_cells))

  samples <- factor(rep(paste0("S", seq_len(n_samples)),
                        each = n_cells / n_samples))

  # Single-covariate H: condition (0 = ctrl, 1 = treated) per sample
  H <- matrix(c(0, 1), nrow = 1,
              dimnames = list("condition",
                              paste0("S", seq_len(n_samples))))

  # C pathway matrix: 0/1 membership, each gene in exactly one pathway
  set.seed(42)
  C <- matrix(0L, nrow = n_genes, ncol = n_pathways,
              dimnames = list(rownames(M),
                              paste0("PW", seq_len(n_pathways))))
  membership <- sample(seq_len(n_pathways), n_genes, replace = TRUE)
  for (i in seq_len(n_genes)) C[i, membership[i]] <- 1L

  model <- CreateGEDIObject(
    Samples = samples,
    M       = M,
    K       = K,
    H       = H,
    C       = C,
    verbose = 0
  )
  model$train(iterations = 10, track_interval = 5)

  list(
    model      = model,
    H          = H,
    C          = C,
    num_cov    = nrow(H),       # == 1
    n_pathways = n_pathways,
    n_genes    = n_genes,
    n_cells    = n_cells
  )
}

# ---- model$diffADB() return value ------------------------------------------

test_that("diffADB returns a num_pathways x N matrix with correct dimnames", {
  fx      <- make_diff_proj_fixture()
  m       <- fx$model
  contrast <- c(condition = 1)  # length 1 == ncol(H.rotation)

  result <- m$diffADB(contrast)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(fx$n_pathways, fx$n_cells))
  expect_equal(rownames(result), paste0("PW", seq_len(fx$n_pathways)))
  expect_equal(colnames(result), m$metadata$cellIDs)
})

test_that("diffADB values are numeric and finite", {
  fx      <- make_diff_proj_fixture()
  m       <- fx$model
  contrast <- c(condition = 1)

  result <- m$diffADB(contrast)
  expect_true(all(is.finite(result)))
})

# ---- plot_features: projection = "diffexp" ----------------------------------

test_that("plot_features returns ggplot for projection='diffexp'", {
  fx       <- make_diff_proj_fixture()
  m        <- fx$model
  contrast <- c(condition = 1)
  genes    <- m$metadata$geneIDs[1:3]
  set.seed(1)
  emb <- matrix(rnorm(fx$n_cells * 2), ncol = 2)

  p <- plot_features(m,
                     features    = genes,
                     embedding   = emb,
                     projection  = "diffexp",
                     contrast    = contrast)
  expect_s3_class(p, "ggplot")
})

test_that("plot_features diffexp with include_O=TRUE returns ggplot", {
  fx       <- make_diff_proj_fixture()
  m        <- fx$model
  contrast <- c(condition = 1)
  genes    <- m$metadata$geneIDs[1:2]
  set.seed(1)
  emb <- matrix(rnorm(fx$n_cells * 2), ncol = 2)

  p <- plot_features(m,
                     features    = genes,
                     embedding   = emb,
                     projection  = "diffexp",
                     contrast    = contrast,
                     include_O   = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_features diffexp works for a single gene (F=1 edge case)", {
  fx       <- make_diff_proj_fixture()
  m        <- fx$model
  contrast <- c(condition = 1)
  set.seed(1)
  emb <- matrix(rnorm(fx$n_cells * 2), ncol = 2)

  p <- plot_features(m,
                     features   = m$metadata$geneIDs[1],
                     embedding  = emb,
                     projection = "diffexp",
                     contrast   = contrast)
  expect_s3_class(p, "ggplot")
})

# ---- plot_features: projection = "diffadb" ----------------------------------

test_that("plot_features returns ggplot for projection='diffadb'", {
  fx       <- make_diff_proj_fixture()
  m        <- fx$model
  contrast <- c(condition = 1)
  pathways <- paste0("PW", 1:3)
  set.seed(1)
  emb <- matrix(rnorm(fx$n_cells * 2), ncol = 2)

  p <- plot_features(m,
                     features   = pathways,
                     embedding  = emb,
                     projection = "diffadb",
                     contrast   = contrast)
  expect_s3_class(p, "ggplot")
})

test_that("plot_features diffadb works for a single pathway (F=1 edge case)", {
  fx       <- make_diff_proj_fixture()
  m        <- fx$model
  contrast <- c(condition = 1)
  set.seed(1)
  emb <- matrix(rnorm(fx$n_cells * 2), ncol = 2)

  p <- plot_features(m,
                     features   = "PW1",
                     embedding  = emb,
                     projection = "diffadb",
                     contrast   = contrast)
  expect_s3_class(p, "ggplot")
})

# ---- Error: contrast = NULL ------------------------------------------------

test_that("plot_features errors when contrast=NULL for projection='diffexp'", {
  fx  <- make_diff_proj_fixture()
  m   <- fx$model
  set.seed(1)
  emb <- matrix(rnorm(fx$n_cells * 2), ncol = 2)

  expect_error(
    plot_features(m, features = m$metadata$geneIDs[1:2],
                  embedding  = emb,
                  projection = "diffexp",
                  contrast   = NULL),
    regexp = "requires a `contrast`"
  )
})

test_that("plot_features errors when contrast=NULL for projection='diffadb'", {
  fx  <- make_diff_proj_fixture()
  m   <- fx$model
  set.seed(1)
  emb <- matrix(rnorm(fx$n_cells * 2), ncol = 2)

  expect_error(
    plot_features(m, features = "PW1",
                  embedding  = emb,
                  projection = "diffadb",
                  contrast   = NULL),
    regexp = "requires a `contrast`"
  )
})

# ---- Error: wrong-length contrast ------------------------------------------

test_that("plot_features errors on wrong-length contrast for diffexp", {
  fx  <- make_diff_proj_fixture()
  m   <- fx$model
  set.seed(1)
  emb <- matrix(rnorm(fx$n_cells * 2), ncol = 2)

  # correct length is 1 (ncol(H.rotation)); supply length 3
  expect_error(
    plot_features(m, features = m$metadata$geneIDs[1:2],
                  embedding  = emb,
                  projection = "diffexp",
                  contrast   = c(1, 0, 0)),
    regexp = "ncol\\(model\\$priors\\$H\\.rotation\\)"
  )
})

test_that("model$diffADB errors on wrong-length contrast", {
  fx  <- make_diff_proj_fixture()
  m   <- fx$model

  expect_error(m$diffADB(c(1, 0)),
               regexp = "number of original H covariates")
})

# ---- Error: diffadb on model with no C prior (P=0) -------------------------

test_that("plot_features errors for diffadb when model has P=0", {
  # Build a model WITHOUT a C matrix (P=0)
  set.seed(99)
  n_genes  <- 20
  n_cells  <- 30
  M <- Matrix::Matrix(
    matrix(stats::rpois(n_genes * n_cells, lambda = 3), n_genes, n_cells),
    sparse = TRUE
  )
  rownames(M) <- paste0("G_", seq_len(n_genes))
  colnames(M) <- paste0("C_", seq_len(n_cells))
  samples <- factor(rep(c("S1", "S2"), each = n_cells / 2))
  H <- matrix(c(0, 1), nrow = 1,
              dimnames = list("cond", c("S1", "S2")))

  m_noC <- CreateGEDIObject(Samples = samples, M = M, K = 2, H = H, verbose = 0)
  m_noC$train(iterations = 5)

  set.seed(1)
  emb <- matrix(rnorm(n_cells * 2), ncol = 2)

  expect_error(
    plot_features(m_noC, features = 1L,
                  embedding  = emb,
                  projection = "diffadb",
                  contrast   = c(1)),
    regexp = "P=0"
  )
})

test_that("model$diffADB errors when model has P=0", {
  set.seed(99)
  n_genes  <- 20
  n_cells  <- 30
  M <- Matrix::Matrix(
    matrix(stats::rpois(n_genes * n_cells, lambda = 3), n_genes, n_cells),
    sparse = TRUE
  )
  rownames(M) <- paste0("G_", seq_len(n_genes))
  colnames(M) <- paste0("C_", seq_len(n_cells))
  samples <- factor(rep(c("S1", "S2"), each = n_cells / 2))
  H <- matrix(c(0, 1), nrow = 1,
              dimnames = list("cond", c("S1", "S2")))

  m_noC <- CreateGEDIObject(Samples = samples, M = M, K = 2, H = H, verbose = 0)
  m_noC$train(iterations = 5)

  expect_error(m_noC$diffADB(c(1)),
               regexp = "gene-level prior")
})
