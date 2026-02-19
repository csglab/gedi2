# tests/testthat/test-plot-utils.R
# -----------------------------------------------------------
# Tests for theme_gedi(), scale_color_gedi_*, scale_fill_gedi_*,
# and compute_color_limits(). No model needed.
# -----------------------------------------------------------

# ===========================================================
# 1. theme_gedi
# ===========================================================
test_that("theme_gedi returns a ggplot2 theme", {
  th <- theme_gedi()
  expect_s3_class(th, "theme")
  expect_s3_class(th, "gg")
})

test_that("theme_gedi respects base_size", {
  th8  <- theme_gedi(base_size = 8)
  th14 <- theme_gedi(base_size = 14)
  expect_s3_class(th8, "theme")
  expect_s3_class(th14, "theme")
})

# ===========================================================
# 2. scale_color_gedi_diverging
# ===========================================================
test_that("scale_color_gedi_diverging returns a ggplot2 scale", {
  sc <- scale_color_gedi_diverging()
  expect_s3_class(sc, "ScaleContinuous")
})

test_that("scale_color_gedi_diverging accepts custom limits", {
  sc <- scale_color_gedi_diverging(limits = c(-5, 5), name = "Test")
  expect_s3_class(sc, "ScaleContinuous")
})

# ===========================================================
# 3. scale_fill_gedi_diverging
# ===========================================================
test_that("scale_fill_gedi_diverging returns a ggplot2 scale", {
  sc <- scale_fill_gedi_diverging()
  expect_s3_class(sc, "ScaleContinuous")
})

# ===========================================================
# 4. scale_color_gedi_discrete
# ===========================================================
test_that("scale_color_gedi_discrete returns a ggplot2 scale", {
  sc <- scale_color_gedi_discrete()
  expect_s3_class(sc, "ScaleDiscrete")
})

# ===========================================================
# 5. compute_color_limits
# ===========================================================
test_that("compute_color_limits returns symmetric limits", {
  vals <- c(-3, -1, 0, 1, 5)
  lim  <- compute_color_limits(vals, symmetric = TRUE, quantile = 1.0)
  expect_length(lim, 2)
  expect_equal(lim[1], -lim[2])  # symmetric about 0
})
 
test_that("compute_color_limits returns asymmetric limits", {
  vals <- c(2, 5, 10, 20)
  lim  <- compute_color_limits(vals, symmetric = FALSE)
  expect_length(lim, 2)
  expect_equal(lim[[1]], 2)
  expect_equal(lim[[2]], 20)
})

test_that("compute_color_limits handles NA values", {
  vals <- c(-1, NA, 3, NA, 5)
  lim  <- compute_color_limits(vals, symmetric = FALSE)
  expect_length(lim, 2)
  expect_false(any(is.na(lim)))
})
