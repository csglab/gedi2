# ==============================================================================
# GEDI Dimensionality Reduction - R Wrappers
# Lazy-computed embedding methods with caching
# ==============================================================================

# ==============================================================================
# Embeddings Accessor R6 Class (for lazy evaluation)
# ==============================================================================

EmbeddingsAccessor <- R6Class(
  "EmbeddingsAccessor",
  public = list(
    initialize = function(gedi_self, gedi_private) {
      private$.gedi_self <- gedi_self
      private$.gedi_private <- gedi_private
    },

    compute_umap = function(...) {
      # Custom UMAP with parameters (not cached)
      args <- list(...)
      if (length(args) > 0) {
        return(compute_umap(private$.gedi_self, private$.gedi_private, ...))
      }
      # If no args, return cached default UMAP (same as $umap active binding)
      return(compute_umap_cached(private$.gedi_self, private$.gedi_private))
    },

    print = function() {
      cat("<GEDI Embeddings Accessor>\n")
      cat("\nAvailable embeddings (lazy-computed):\n")
      cat("  $svd   - Factorized SVD\n")
      cat("  $pca   - PCA coordinates\n")
      cat("  $umap  - UMAP embedding (default parameters, cached)\n")
      cat("\nFor custom UMAP parameters, use:\n")
      cat("  $compute_umap(n_neighbors=30, min_dist=0.5, ...)\n")
      cat("\nAccess with: model$embeddings$pca or model$embeddings$umap\n")
      invisible(self)
    }
  ),
  private = list(
    .gedi_self = NULL,
    .gedi_private = NULL
  ),
  active = list(
    svd = function(value) {
      if (!missing(value)) stop("svd is read-only", call. = FALSE)
      compute_svd_factorized(private$.gedi_self, private$.gedi_private)
    },

    pca = function(value) {
      if (!missing(value)) stop("pca is read-only", call. = FALSE)
      compute_pca(private$.gedi_self, private$.gedi_private)
    },

    umap = function(value) {
      if (!missing(value)) stop("umap is read-only", call. = FALSE)
      compute_umap_cached(private$.gedi_self, private$.gedi_private)
    }
  )
)


#' Compute Factorized SVD with Caching
#'
#' @description
#' Computes factorized SVD preserving GEDI structure: SVD(Z) x SVD(middle) x SVD(DB).
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
#'     \item u: Left singular vectors (J x K matrix) with gene names
#'     \item v: Right singular vectors (N x K matrix) with cell names
#'   }
#'
#' @keywords internal
#' @noRd
compute_svd_factorized <- function(self, private, force_recompute = FALSE) {
  
  # Check cache
  if (!force_recompute && !is.null(private$.cache$svd)) {
    log_cached("SVD", private$.verbose)
    return(private$.cache$svd)
  }

  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }

  log_start("Factorized SVD", private$.verbose)

  # Call C++ function (silent)
  svd_result <- compute_svd_factorized_cpp(
    Z = self$params$Z,
    D = self$params$D,
    Bi_list = self$params$Bi,
    verbose = 0  # C++ silent
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

  details <- sprintf("%d genes x %d cells x %d factors",
                     nrow(svd_result$u), nrow(svd_result$v), length(svd_result$d))
  log_complete("Factorized SVD", details, private$.verbose)

  return(svd_result)
}


#' Compute PCA Coordinates with Caching
#'
#' @description
#' Computes PCA coordinates as V x diag(d) from the factorized SVD.
#' Uses cached SVD result if available.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param force_recompute Logical, if TRUE bypasses cache and recomputes
#'
#' @return Matrix (N x K) with PCA coordinates, rows = cells, columns = PCs
#'
#' @keywords internal
#' @noRd
compute_pca <- function(self, private, force_recompute = FALSE) {
  
  # Check cache
  if (!force_recompute && !is.null(private$.cache$pca)) {
    log_cached("PCA", private$.verbose)
    return(private$.cache$pca)
  }

  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }

  log_start("PCA coordinates", private$.verbose)

  # Get SVD (uses cache if available)
  svd_result <- compute_svd_factorized(self, private)

  # PCA = V x diag(d)
  pca_coords <- svd_result$v %*% diag(svd_result$d, nrow = length(svd_result$d))

  # Set dimension names
  rownames(pca_coords) <- private$.cellIDs
  colnames(pca_coords) <- paste0("PC", seq_len(ncol(pca_coords)))

  # Cache result
  if (is.null(private$.cache)) {
    private$.cache <- list()
  }
  private$.cache$pca <- pca_coords

  details <- format_dims(nrow(pca_coords), ncol(pca_coords), "cells", "PCs")
  log_complete("PCA coordinates", details, private$.verbose)

  return(pca_coords)
}


#' Compute UMAP Embedding with Caching (Default Parameters Only)
#'
#' @description
#' Computes UMAP embedding with default parameters and caches the result.
#' Only the default parameter UMAP is cached to save memory.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param force_recompute Logical, if TRUE bypasses cache and recomputes
#'
#' @return Matrix (N x 2) with UMAP coordinates, rows = cells
#'
#' @keywords internal
#' @noRd
compute_umap_cached <- function(self, private, force_recompute = FALSE) {

  # Check cache
  if (!force_recompute && !is.null(private$.cache$umap)) {
    log_cached("UMAP", private$.verbose)
    return(private$.cache$umap)
  }

  # Compute UMAP with default parameters
  umap_result <- compute_umap(self, private,
                               input = "pca",
                               n_neighbors = 15,
                               min_dist = 0.1,
                               n_components = 2,
                               metric = "euclidean",
                               n_threads = 0)

  # Cache result
  if (is.null(private$.cache)) {
    private$.cache <- list()
  }
  private$.cache$umap <- umap_result

  return(umap_result)
}


#' Compute UMAP Embedding
#'
#' @description
#' Computes UMAP embedding using the uwot package. Input can be PCA coordinates,
#' DB projection, or ZDB projection. This function does NOT cache results,
#' allowing for custom parameters. Use compute_umap_cached() for cached default UMAP.
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
#' @return Matrix (N x n_components) with UMAP coordinates, rows = cells
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
  
  log_start("UMAP embedding", private$.verbose)

  # Get input data
  if (input == "pca") {
    log_info(sprintf("Input: PCA coordinates"), private$.verbose)
    X <- compute_pca(self, private)

  } else if (input == "db") {
    log_info(sprintf("Input: DB projection"), private$.verbose)
    X <- t(compute_DB(self, private))

  } else if (input == "zdb") {
    log_info(sprintf("Input: ZDB projection"), private$.verbose)
    X <- t(compute_ZDB(self, private))

  } else {
    stop("input must be 'pca', 'db', or 'zdb'", call. = FALSE)
  }

  # Log parameters
  log_info(sprintf("Parameters: n_neighbors=%d, min_dist=%.2f, n_components=%d",
                   n_neighbors, min_dist, n_components), private$.verbose)

  # Run UMAP
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

  details <- format_dims(nrow(umap_result), ncol(umap_result), "cells", "dimensions")
  log_complete("UMAP embedding", details, private$.verbose)

  return(umap_result)
}


#' Clear Dimensionality Reduction Cache
#'
#' @description
#' Clears cached embedding results (SVD, PCA, UMAP).
#' Useful after model parameters have been updated.
#'
#' @param private Reference to private environment
#' @param what Character vector specifying which cache entries to clear.
#'   Options: "svd", "pca", "umap", or NULL to clear all (default: NULL)
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
    private$.cache$umap <- NULL
    if (private$.verbose > 0) {
      message("Dimensionality reduction cache cleared")
    }
  } else {
    # Clear specific entries
    valid_entries <- c("svd", "pca", "umap")
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
    return(c(svd = FALSE, pca = FALSE, umap = FALSE))
  }

  status <- c(
    svd = !is.null(private$.cache$svd),
    pca = !is.null(private$.cache$pca),
    umap = !is.null(private$.cache$umap)
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
    return(c(svd = 0, pca = 0, umap = 0, Total = 0))
  }

  memory_mb <- c(
    svd = if (!is.null(private$.cache$svd)) {
      object.size(private$.cache$svd) / 1024^2
    } else 0,
    pca = if (!is.null(private$.cache$pca)) {
      object.size(private$.cache$pca) / 1024^2
    } else 0,
    umap = if (!is.null(private$.cache$umap)) {
      object.size(private$.cache$umap) / 1024^2
    } else 0
  )

  memory_mb["Total"] <- sum(memory_mb)

  return(memory_mb)
}
