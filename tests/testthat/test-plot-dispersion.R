# tests/testthat/test-plot-dispersion.R
# -----------------------------------------------------------
# Tests for plot_dispersion(). Can use a hand-crafted data.frame.
# -----------------------------------------------------------

test_that("plot_dispersion returns a ggplot with valid data.frame", {
  disp_df <- data.frame(
    Expected_Var = runif(100, 0.1, 10),
    Observed_Var = runif(100, 0.1, 10),
    Sample       = rep(c("S1", "S2"), each = 50)
  )
  p <- plot_dispersion(disp_df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_dispersion errors on missing columns", {
  bad_df <- data.frame(x = 1:10, y = 1:10)
  expect_error(plot_dispersion(bad_df))
})

test_that("plot_dispersion show_identity parameter works", {
  disp_df <- data.frame(
    Expected_Var = runif(50, 0.1, 5),
    Observed_Var = runif(50, 0.1, 5),
    Sample       = rep("S1", 50)
  )
  p1 <- plot_dispersion(disp_df, show_identity = TRUE)
  p2 <- plot_dispersion(disp_df, show_identity = FALSE)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})
