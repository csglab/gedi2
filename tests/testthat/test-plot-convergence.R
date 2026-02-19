# tests/testthat/test-plot-convergence.R
# -----------------------------------------------------------
# Tests for plot_convergence(). Requires a trained model.
# -----------------------------------------------------------

test_that("plot_convergence returns a ggplot (faceted layout)", {
  model <- make_trained_model(K = 5, iterations = 20)
  p <- plot_convergence(model, layout = "faceted")
  expect_s3_class(p, "ggplot")
})

test_that("plot_convergence returns a list of ggplots (separate layout)", {
  model <- make_trained_model(K = 5, iterations = 20)
  plots <- plot_convergence(model, layout = "separate")
  expect_type(plots, "list")
  expect_true(length(plots) > 0)
  # Each element should be a ggplot
  for (p in plots) {
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_convergence errors on untrained model", {
  dat <- make_synthetic_data()
  model <- CreateGEDIObject(
    Samples = dat$sample_labels, M = dat$M,
    K = 5, mode = "Bsphere", verbose = 0
  )
  expect_error(plot_convergence(model))
})
