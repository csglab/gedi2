# tests/testthat/test-plot-features.R
# -----------------------------------------------------------
# Tests for plot_features() and plot_feature_ratio().
# Require a trained model + custom embedding matrix.
# -----------------------------------------------------------

test_that("plot_features returns a ggplot with integer indices", {
  model <- make_trained_model(K = 5, iterations = 10)
  # Use a random 2-column matrix as embedding to avoid uwot dependency
  set.seed(1)
  emb <- matrix(rnorm(400 * 2), ncol = 2)

  p <- plot_features(model, features = c(1, 2, 3), embedding = emb)
  expect_s3_class(p, "ggplot")
})

test_that("plot_features returns a ggplot with gene names", {
  model <- make_trained_model(K = 5, iterations = 10)
  set.seed(1)
  emb <- matrix(rnorm(400 * 2), ncol = 2)
  genes <- model$metadata$geneIDs[1:2]

  p <- plot_features(model, features = genes, embedding = emb)
  expect_s3_class(p, "ggplot")
})

test_that("plot_features errors on invalid gene name", {
  model <- make_trained_model(K = 5, iterations = 10)
  set.seed(1)
  emb <- matrix(rnorm(400 * 2), ncol = 2)
  expect_error(plot_features(model, features = "NONEXISTENT_GENE", embedding = emb))
})

test_that("plot_feature_ratio returns a ggplot with integer indices", {
  model <- make_trained_model(K = 5, iterations = 10)
  set.seed(1)
  emb <- matrix(rnorm(400 * 2), ncol = 2)

  p <- plot_feature_ratio(model, gene1 = 1, gene2 = 2, embedding = emb)
  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_ratio errors on untrained model", {
  dat <- make_synthetic_data()
  model <- CreateGEDIObject(
    Samples = dat$sample_labels, M = dat$M,
    K = 5, mode = "Bsphere", verbose = 0
  )
  set.seed(1)
  emb <- matrix(rnorm(400 * 2), ncol = 2)
  expect_error(plot_feature_ratio(model, gene1 = 1, gene2 = 2, embedding = emb))
})
