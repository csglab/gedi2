# ==============================================================================
# GEDI Projections - R Wrappers
# Lazy-computed projection methods with caching
# ==============================================================================

#' Compute ZDB Projection with Caching
#'
#' @description
#' Computes the shared manifold projection ZDB = Z * diag(D) * B, where B is
#' the concatenation of all sample-specific Bi matrices. Results are cached
#' automatically for subsequent access.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param force_recompute Logical, if TRUE bypasses cache and recomputes
#'
#' @return Dense matrix (J × N) with row/column names set to geneIDs/cellIDs
#'
#' @keywords internal
#' @noRd
compute_ZDB <- function(self, private, force_recompute = FALSE) {
  
  # Check if cached and not forcing recompute
  if (!force_recompute && !is.null(private$.cache$ZDB)) {
    if (private$.verbose > 0) {
      cat("✓ Using cached ZDB\n")
    }
    return(private$.cache$ZDB)
  }
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Call C++ function
  ZDB <- compute_ZDB_cpp(
    Z = self$params$Z,
    D = self$params$D,
    Bi_list = self$params$Bi,
    verbose = private$.verbose
  )
  
  # Add row/column names
  rownames(ZDB) <- private$.geneIDs
  colnames(ZDB) <- private$.cellIDs
  
  # Cache result
  if (is.null(private$.cache)) {
    private$.cache <- list()
  }
  private$.cache$ZDB <- ZDB
  
  return(ZDB)
}


#' Compute DB Projection with Caching
#'
#' @description
#' Computes the latent factor embedding DB = diag(D) * B. This is the
#' K-dimensional representation of cells in latent factor space.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param force_recompute Logical, if TRUE bypasses cache and recomputes
#'
#' @return Dense matrix (K × N) with row names "LV1", "LV2", ... and
#'   column names set to cellIDs
#'
#' @keywords internal
#' @noRd
compute_DB <- function(self, private, force_recompute = FALSE) {
  
  # Check cache
  if (!force_recompute && !is.null(private$.cache$DB)) {
    if (private$.verbose > 0) {
      cat("✓ Using cached DB\n")
    }
    return(private$.cache$DB)
  }
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Call C++ function
  DB <- compute_DB_cpp(
    D = self$params$D,
    Bi_list = self$params$Bi,
    verbose = private$.verbose
  )
  
  # Add row/column names
  rownames(DB) <- paste0("LV", seq_len(self$aux$K))
  colnames(DB) <- private$.cellIDs
  
  # Cache result
  if (is.null(private$.cache)) {
    private$.cache <- list()
  }
  private$.cache$DB <- DB
  
  return(DB)
}


#' Compute ADB Projection with Caching
#'
#' @description
#' Computes the pathway activity projection ADB = C.rotation * A * diag(D) * B.
#' Only available when gene-level prior matrix C was provided during model setup.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param force_recompute Logical, if TRUE bypasses cache and recomputes
#'
#' @return Dense matrix (P × N) with row names from original C matrix column
#'   names and column names set to cellIDs
#'
#' @keywords internal
#' @noRd
compute_ADB <- function(self, private, force_recompute = FALSE) {
  
  # Check cache
  if (!force_recompute && !is.null(private$.cache$ADB)) {
    if (private$.verbose > 0) {
      cat("✓ Using cached ADB\n")
    }
    return(private$.cache$ADB)
  }
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if C prior was provided
  if (self$aux$P == 0) {
    stop("Cannot compute ADB: no gene-level prior (C) was provided during setup.", 
         call. = FALSE)
  }
  
  # Check if aux_static exists
  if (is.null(private$.aux_static)) {
    stop("Missing auxiliary data. Model may not be properly initialized.", 
         call. = FALSE)
  }
  
  # Call C++ function
  ADB <- compute_ADB_cpp(
    C_rotation = private$.aux_static$C.rotation,
    A = self$params$A,
    D = self$params$D,
    Bi_list = self$params$Bi,
    verbose = private$.verbose
  )
  
  # Add row/column names
  rownames(ADB) <- colnames(private$.aux_static$inputC)
  colnames(ADB) <- private$.cellIDs
  
  # Cache result
  if (is.null(private$.cache)) {
    private$.cache <- list()
  }
  private$.cache$ADB <- ADB
  
  return(ADB)
}


#' Compute Differential Expression (Cell-Space JxN)
#'
#' @description
#' Computes the cell-specific differential expression effect (J x N).
#' This is the R wrapper for the C++ function getDiffExp_cpp().
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param contrast Numeric vector of length L specifying the contrast
#' @param include_O Logical, if TRUE adds the global offset effect (diffO)
#'
#' @return Dense matrix (J × N) of differential expression values
#'
#' @keywords internal
#' @noRd
compute_diffQ <- function(self, private, contrast, include_O = FALSE) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if H prior was provided
  if (self$aux$L == 0) {
    stop("Cannot compute diffQ: no sample-level prior (H) was provided during setup.", 
         call. = FALSE)
  }
  
  # Validate contrast
  if (!is.numeric(contrast) || length(contrast) != self$aux$L) {
    stop("contrast must be a numeric vector of length L (", self$aux$L, ")", 
         call. = FALSE)
  }
  
  # Check if aux_static exists
  if (is.null(private$.aux_static)) {
    stop("Missing auxiliary data. Model may not be properly initialized.", 
         call. = FALSE)
  }
  
  # Call C++ function (renamed to match old getDiffExp.gedi logic)
  diffExp <- getDiffExp_cpp(
    Rk_list = self$params$Rk,
    H_rotation = private$.aux_static$H.rotation,
    contrast = contrast,
    D = self$params$D,
    Bi_list = self$params$Bi,
    verbose = private$.verbose
  )
  
  # Add diffO if requested
  if (include_O) {
    diffO_vec <- compute_diffO(self, private, contrast) # Calls the wrapper below
    
    # Add diffO to each column (broadcast)
    diffExp <- sweep(diffExp, 1, diffO_vec, "+")
  }
  
  # Add row/column names
  rownames(diffExp) <- private$.geneIDs
  colnames(diffExp) <- private$.cellIDs
  
  return(diffExp)
}


#' Compute Differential O (Global Offset) Effect
#'
#' @description
#' Computes the global differential effect on gene-specific offsets (J x 1).
#' This is the R wrapper for the C++ function getDiffO_cpp().
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param contrast Numeric vector of length L specifying the contrast
#'
#' @return Numeric vector of length J
#'
#' @keywords internal
#' @noRd
compute_diffO <- function(self, private, contrast) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if H prior was provided
  if (self$aux$L == 0) {
    stop("Cannot compute diffO: no sample-level prior (H) was provided during setup.", 
         call. = FALSE)
  }
  
  # Validate contrast
  if (!is.numeric(contrast) || length(contrast) != self$aux$L) {
    stop("contrast must be a numeric vector of length L (", self$aux$L, ")", 
         call. = FALSE)
  }
  
  # Check if aux_static exists
  if (is.null(private$.aux_static)) {
    stop("Missing auxiliary data. Model may not be properly initialized.", 
         call. = FALSE)
  }
  
  # Call C++ function
  diffO <- getDiffO_cpp(
    Ro = self$params$Ro,
    H_rotation = private$.aux_static$H.rotation,
    contrast = contrast,
    verbose = private$.verbose
  )
  
  # Add names
  names(diffO) <- private$.geneIDs
  
  return(diffO)
}


#' Compute Differential Expression (Alias)
#'
#' @description
#' Convenience alias for compute_diffQ.
#'
#' @keywords internal
#' @noRd
compute_diffExp <- function(self, private, contrast, include_O = FALSE) {
  # This is just an alias for compute_diffQ
  compute_diffQ(self, private, contrast, include_O)
}


#' Clear Projection Cache
#'
#' @description
#' Clears cached projection results.
#'
#' @keywords internal
#' @noRd
clear_projection_cache <- function(private, what = NULL) {
  
  if (is.null(private$.cache)) {
    return(invisible(NULL))
  }
  
  if (is.null(what)) {
    # Clear all cache
    private$.cache <- list()
    if (private$.verbose > 0) {
      message("Projection cache cleared")
    }
  } else {
    # Clear specific entries
    valid_entries <- c("ZDB", "DB", "ADB")
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


#' Get Projection Cache Status
#'
#' @description
#' Returns information about which projections are currently cached.
#'
#' @keywords internal
#' @noRd
get_cache_status <- function(private) {
  
  if (is.null(private$.cache)) {
    return(c(ZDB = FALSE, DB = FALSE, ADB = FALSE))
  }
  
  status <- c(
    ZDB = !is.null(private$.cache$ZDB),
    DB = !is.null(private$.cache$DB),
    ADB = !is.null(private$.cache$ADB)
  )
  
  return(status)
}


#' Get Projection Cache Memory Usage
#'
#' @description
#' Estimates memory usage of cached projections in MB.
#'
#' @keywords internal
#' @noRd
get_cache_memory <- function(private) {
  
  if (is.null(private$.cache)) {
    return(c(ZDB = 0, DB = 0, ADB = 0, Total = 0))
  }
  
  memory_mb <- c(
    ZDB = if (!is.null(private$.cache$ZDB)) {
      object.size(private$.cache$ZDB) / 1024^2
    } else 0,
    DB = if (!is.null(private$.cache$DB)) {
      object.size(private$.cache$DB) / 1024^2
    } else 0,
    ADB = if (!is.null(private$.cache$ADB)) {
      object.size(private$.cache$ADB) / 1024^2
    } else 0
  )
  
  memory_mb["Total"] <- sum(memory_mb)
  
  return(memory_mb)
}