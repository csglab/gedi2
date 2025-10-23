# ==============================================================================
# GEDI Dimensionality Reduction - R Wrappers
# Lazy-computed embedding methods with caching
# ==============================================================================

#' Compute Factorized SVD with Caching
#'
#' @description
#' Computes factorized SVD preserving GEDI structure: SVD(Z) × SVD(middle) × SVD(DB).
#' This approach maintains biological interpretability by respecting the decomposition
#' structure. Results are cached automatically for subsequent access.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param force_recompute Logical, if TRUE bypasses cache and recomputes
#'
#' @return List with three components:
#'   \itemize{
#'     \item d: Singular values (length K vector)
#'     \item u: Left singular vectors (J × K matrix) with gene names
#'     \item v: Right singular vectors (N × K matrix) with cell names
#'   }
#'
#' @keywords internal
#' @noRd
compute_svd_factorized <- function(self, private, force_recompute = FALSE) {
  
  # Check cache
  if (!force_recompute && !is.null(private$.cache$svd)) {
    if (private$.verbose > 0) {
      cat("Using cached SVD\n")
    }
    return(private$.cache$svd)
  }
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Call C++ function
  svd_result <- compute_svd_factorized_cpp(
    Z = self$params$Z,
    D = self$params$D,
    Bi_list = self$params$Bi,
    verbose = private$.verbose
  )
  
  # Add dimension names
  names(svd_result$d) <- paste0("LV", seq_len(self$aux$K))
  
  rownames(svd_result$u) <- private$.geneIDs
  colnames(svd_result$u) <- paste0("LV", seq_len(self$aux$K))
  
  rownames(svd_result$v) <- private$.cellIDs
  colnames(svd_result$v) <- paste0("LV", seq_len(self$aux$K))
  
  # Cache result
  if (is.null(private$.cache)) {
    private$.cache <- list()
  }
  private$.cache$svd <- svd_result
  
  return(svd_result)
}


#' Compute PCA Coordinates with Caching
#'
#' @description
#' Computes PCA coordinates as V × diag(d) from the factorized SVD.
#' Uses cached SVD result if available.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param force_recompute Logical, if TRUE bypasses cache and recomputes
#'
#' @return Matrix (N × K) with PCA coordinates, rows = cells, columns = PCs
#'
#' @keywords internal
#' @noRd
compute_pca <- function(self, private, force_recompute = FALSE) {
  
  # Check cache
  if (!force_recompute && !is.null(private$.cache$pca)) {
    if (private$.verbose > 0) {
      cat("Using cached PCA\n")
    }
    return(private$.cache$pca)
  }
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Get SVD (uses cache if available)
  svd_result <- compute_svd_factorized(self, private)
  
  if (private$.verbose >= 1) {
    cat("Computing PCA coordinates from SVD...\n")
  }
  
  # PCA = V × diag(d)
  pca_coords <- svd_result$v %*% diag(svd_result$d, nrow = length(svd_result$d))
  
  # Set dimension names
  rownames(pca_coords) <- private$.cellIDs
  colnames(pca_coords) <- paste0("PC", seq_len(ncol(pca_coords)))
  
  # Cache result
  if (is.null(private$.cache)) {
    private$.cache <- list()
  }
  private$.cache$pca <- pca_coords
  
  if (private$.verbose >= 1) {
    cat("✓ PCA computed: ", nrow(pca_coords), " cells × ", 
        ncol(pca_coords), " PCs\n", sep = "")
  }
  
  return(pca_coords)
}


#' Compute UMAP Embedding
#'
#' @description
#' Computes UMAP embedding using the uwot package. Input can be PCA coordinates,
#' DB projection, or ZDB projection. UMAP results are NOT cached as parameters vary.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param input Character, input data type: "pca" (default), "db", or "zdb"
#' @param n_neighbors Integer, size of local neighborhood (default: 15)
#' @param min_dist Numeric, minimum distance between points in embedding (default: 0.1)
#' @param n_components Integer, dimensionality of output (default: 2)
#' @param metric Character, distance metric (default: "euclidean")
#' @param n_threads Integer, number of threads (default: 0 = auto)
#' @param ... Additional arguments passed to uwot::umap()
#'
#' @return Matrix (N × n_components) with UMAP coordinates, rows = cells
#'
#' @keywords internal
#' @noRd
compute_umap <- function(self, private, 
                         input = "pca",
                         n_neighbors = 15,
                         min_dist = 0.1,
                         n_components = 2,
                         metric = "euclidean",
                         n_threads = 0,
                         ...) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check uwot package
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' required for UMAP.\n",
         "  Install with: install.packages('uwot')",
         call. = FALSE)
  }
  
  # Get input data
  if (input == "pca") {
    if (private$.verbose >= 1) {
      cat("Using PCA coordinates as UMAP input\n")
    }
    X <- compute_pca(self, private)
    
  } else if (input == "db") {
    if (private$.verbose >= 1) {
      cat("Using DB projection as UMAP input\n")
    }
    X <- t(compute_DB(self, private))
    
  } else if (input == "zdb") {
    if (private$.verbose >= 1) {
      cat("Using ZDB projection as UMAP input\n")
    }
    X <- t(compute_ZDB(self, private))
    
  } else {
    stop("input must be 'pca', 'db', or 'zdb'", call. = FALSE)
  }
  
  # Run UMAP
  if (private$.verbose >= 1) {
    cat("Computing UMAP embedding...\n")
    cat("  Parameters: n_neighbors=", n_neighbors, 
        ", min_dist=", min_dist,
        ", n_components=", n_components, "\n", sep = "")
  }
  
  umap_result <- uwot::umap(
    X,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    n_components = n_components,
    metric = metric,
    n_threads = if (n_threads > 0) n_threads else NULL,
    ret_model = FALSE,
    verbose = private$.verbose > 0,
    ...
  )
  
  # Set dimension names
  rownames(umap_result) <- private$.cellIDs
  colnames(umap_result) <- paste0("UMAP", seq_len(n_components))
  
  if (private$.verbose >= 1) {
    cat("✓ UMAP computed: ", nrow(umap_result), " cells × ", 
        ncol(umap_result), " dimensions\n", sep = "")
  }
  
  return(umap_result)
}


#' Clear Dimensionality Reduction Cache
#'
#' @description
#' Clears cached embedding results (SVD, PCA). UMAP is not cached.
#' Useful after model parameters have been updated.
#'
#' @param private Reference to private environment
#' @param what Character vector specifying which cache entries to clear.
#'   Options: "svd", "pca", or NULL to clear all (default: NULL)
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
clear_dimred_cache <- function(private, what = NULL) {
  
  if (is.null(private$.cache)) {
    return(invisible(NULL))
  }
  
  if (is.null(what)) {
    # Clear all dimred cache
    private$.cache$svd <- NULL
    private$.cache$pca <- NULL
    if (private$.verbose > 0) {
      message("Dimensionality reduction cache cleared")
    }
  } else {
    # Clear specific entries
    valid_entries <- c("svd", "pca")
    invalid <- setdiff(what, valid_entries)
    if (length(invalid) > 0) {
      warning("Invalid cache entries: ", paste(invalid, collapse = ", "))
    }
    
    for (entry in intersect(what, valid_entries)) {
      private$.cache[[entry]] <- NULL
      if (private$.verbose > 0) {
        message(sprintf("Cache cleared: %s", entry))
      }
    }
  }
  
  invisible(NULL)
}


#' Get Dimensionality Reduction Cache Status
#'
#' @description
#' Returns information about which embeddings are currently cached.
#'
#' @param private Reference to private environment
#'
#' @return Named logical vector indicating which embeddings are cached
#'
#' @keywords internal
#' @noRd
get_dimred_cache_status <- function(private) {
  
  if (is.null(private$.cache)) {
    return(c(svd = FALSE, pca = FALSE))
  }
  
  status <- c(
    svd = !is.null(private$.cache$svd),
    pca = !is.null(private$.cache$pca)
  )
  
  return(status)
}


#' Get Dimensionality Reduction Cache Memory Usage
#'
#' @description
#' Estimates memory usage of cached embeddings in MB.
#'
#' @param private Reference to private environment
#'
#' @return Named numeric vector with memory usage in MB for each cached embedding
#'
#' @keywords internal
#' @noRd
get_dimred_cache_memory <- function(private) {
  
  if (is.null(private$.cache)) {
    return(c(svd = 0, pca = 0, Total = 0))
  }
  
  memory_mb <- c(
    svd = if (!is.null(private$.cache$svd)) {
      object.size(private$.cache$svd) / 1024^2
    } else 0,
    pca = if (!is.null(private$.cache$pca)) {
      object.size(private$.cache$pca) / 1024^2
    } else 0
  )
  
  memory_mb["Total"] <- sum(memory_mb)
  
  return(memory_mb)
}
