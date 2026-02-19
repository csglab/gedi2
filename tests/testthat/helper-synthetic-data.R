# tests/testthat/helper-synthetic-data.R
# Shared synthetic data generator for all test files.
# Sourced automatically by testthat before any test-*.R file.

#' Create a small synthetic negative-binomial count matrix with sample labels
#'
#' @param n_genes Number of genes (rows)
#' @param n_cells Number of cells (columns)
#' @param n_samples Number of biological samples
#' @param seed Random seed for reproducibility
#' @return A list with components: M (sparse dgCMatrix), sample_labels (factor)
make_synthetic_data <- function(n_genes = 200,
                                n_cells = 400,
                                n_samples = 4,
                                seed = 42) {
  set.seed(seed)
  cells_per_sample <- n_cells / n_samples
  stopifnot(n_cells %% n_samples == 0)

  sample_labels <- factor(rep(paste0("S", seq_len(n_samples)),
                              each = cells_per_sample))

  gene_means  <- runif(n_genes, min = 0.5, max = 8)
  dispersion  <- 0.5

  counts <- matrix(
    stats::rnbinom(n_genes * n_cells,
                   size = dispersion,
                   mu   = rep(gene_means, n_cells)),
    nrow = n_genes,
    ncol = n_cells
  )
  rownames(counts) <- paste0("Gene_", seq_len(n_genes))
  colnames(counts) <- paste0("Cell_", seq_len(n_cells))

  M <- Matrix::Matrix(counts, sparse = TRUE)

  list(M = M, sample_labels = sample_labels)
}

#' Train a small GEDI model on synthetic data (cached per session)
#'
#' @param K Latent factors
#' @param iterations Training iterations
#' @return Trained GEDI model object
make_trained_model <- function(K = 5, iterations = 20) {
  dat   <- make_synthetic_data()
  model <- CreateGEDIObject(
    Samples = dat$sample_labels,
    M       = dat$M,
    K       = K,
    mode    = "Bsphere",
    verbose = 0
  )
  model$train(iterations = iterations)
  model
}
