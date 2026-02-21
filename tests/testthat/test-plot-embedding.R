# tests/testthat/test-plot-embedding.R
# -----------------------------------------------------------
# Tests for plot_embedding() using the backwards-compatible
# (matrix) API — no model or optional packages needed.
# -----------------------------------------------------------

test_that("plot_embedding works with matrix + continuous color (old API)", {
  set.seed(1)
  emb  <- matrix(rnorm(200), ncol = 2)
  vals <- rnorm(100)
  p <- plot_embedding(emb, color = vals)
  expect_s3_class(p, "ggplot")
})

test_that("plot_embedding works with matrix + discrete color (old API)", {
  set.seed(1)
  emb    <- matrix(rnorm(200), ncol = 2)
  groups <- rep(c("A", "B"), each = 50)
  p <- plot_embedding(emb, color = groups)
  expect_s3_class(p, "ggplot")
})

test_that("plot_embedding works with no color (old API)", {
  set.seed(1)
  emb <- matrix(rnorm(200), ncol = 2)
  p <- plot_embedding(emb)
  expect_s3_class(p, "ggplot")
})

test_that("plot_embedding errors on wrong embedding dimensions", {
  emb <- matrix(rnorm(300), ncol = 3)
  expect_error(plot_embedding(emb))
})

test_that("plot_embedding errors on mismatched color length", {
  emb <- matrix(rnorm(200), ncol = 2)
  expect_error(plot_embedding(emb, color = 1:5))
})
