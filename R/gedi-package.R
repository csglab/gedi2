#' @keywords internal
"_PACKAGE"

#' gedi: Gene Expression Data Integration
#'
#' @name gedi-package
#' @description
#' A memory-efficient implementation for integrating gene expression data from
#' single-cell RNA sequencing experiments. GEDI v2 uses a high-performance C++
#' backend with thin R wrappers to enable analysis of large-scale single-cell
#' datasets with minimal memory overhead.
#'
#' @details
#' ## Key Features:
#' * **Memory-efficient**: All data lives in C++ backend; R objects are ~1 KB
#' * **Multiple data modalities**: Count matrices (M), paired data (CITE-seq),
#'   binary indicators (X), or pre-processed expression (Y)
#' * **Latent variable model**: Dimensionality reduction with batch effect correction
#' * **High performance**: OpenMP parallelization with optimized C++ backend
#' * **Sparse matrix support**: Efficiently handles sparse single-cell data
#'
#' ## Main Function:
#' The primary interface is \code{\link{CreateGEDIObject}}, which creates a GEDI model
#' from expression data.
#'
#' ## Workflow:
#' \enumerate{
#'   \item Create model: \code{model <- CreateGEDIObject(Samples, M, K)}
#'   \item Train model: \code{model$train(iterations = 50)}
#'   \item Access results: \code{Z <- model$Z}, \code{params <- model$params}
#' }
#'
#' @section Architecture:
#' GEDI v2 implements a three-layer architecture:
#' \itemize{
#'   \item \strong{C++ Core}: Stateful GEDI class with full optimization
#'   \item \strong{R6 Wrapper}: Thin R6 class exposing methods and active bindings
#'   \item \strong{Factory Function}: \code{CreateGEDIObject()} for user-friendly creation
#' }
#'
#' @section Computational Requirements:
#' \itemize{
#'   \item \strong{R}: >= 4.0.0
#'   \item \strong{C++ Compiler}: C++14 or later (default in R >= 4.0)
#'   \item \strong{Eigen}: >= 3.3.0 (linear algebra library)
#'   \item \strong{OpenMP}: Optional, for parallelization
#' }
#'
#' @author
#' Computational and Statistical Genomics Laboratory, McGill University
#'
#' @references
#' Add your publication reference here when available.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{CreateGEDIObject}}: Create a GEDI model
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data("pbmc_small", package = "Seurat")
#'
#' # Create GEDI model
#' model <- CreateGEDIObject(
#'   Samples = pbmc_small@meta.data$orig.ident,
#'   M = pbmc_small@assays$RNA@counts,
#'   K = 10,
#'   verbose = 1
#' )
#'
#' # Train model
#' model$train(iterations = 50, track_interval = 5)
#'
#' # Access latent representation
#' Z <- model$Z
#'
#' # View model summary
#' print(model)
#' }
#'
#' @useDynLib gedi2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import R6
#' @import Matrix
#' @importFrom methods as is
#' @importFrom stats coef median rnorm runif var
#' @importFrom utils object.size setTxtProgressBar txtProgressBar
NULL
