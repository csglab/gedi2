# ==============================================================================
# GEDI Imputation Accessors - R6 Methods
# User-facing interface for imputation, variance, and dispersion analysis
# ==============================================================================

#' Create Imputation Accessor Object
#'
#' @description
#' Creates accessor object for imputation analysis with Y, variance, and
#' dispersion methods. This provides a clean user interface: model$imputed$Y(M)
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#'
#' @return S3 object of class "gedi_imputation"
#'
#' @keywords internal
#' @noRd
create_imputation_accessor <- function(self, private) {
  
  structure(
    list(
      Y = function(M = NULL, logScale = TRUE, rowCentre = TRUE) {
        get_imputed_Y(self, private, M, logScale, rowCentre)
      },
      
      variance = function(M) {
        get_Y_variance(self, private, M)
      },
      
      dispersion = function(M, subsample = 1e6) {
        get_dispersion(self, private, M, subsample)
      }
    ),
    class = "gedi_imputation"
  )
}


#' Get Imputed Y Expression (with caching)
#'
#' @description
#' Returns imputed gene expression with sample-specific effects removed.
#' Results are cached for faster repeated access.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param M Count matrix (optional - only needed for non-log scale output)
#' @param logScale Logical, if TRUE returns log-scale values (default)
#' @param rowCentre Logical, if TRUE removes global gene offset o (default)
#'
#' @return Dense matrix (J × N) with imputed expression values
#'
#' @keywords internal
#' @noRd
get_imputed_Y <- function(self, private, M = NULL, logScale = TRUE, rowCentre = TRUE) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  obs_type <- self$aux$obs.type
  
  # Create cache key
  cache_key <- paste0("Y_imputed_", logScale, "_", rowCentre)
  
  # Check if cached
  if (!is.null(private$.imputation_cache) &&
      !is.null(private$.imputation_cache[[cache_key]])) {
    log_cached("Imputed Y", private$.verbose)
    return(private$.imputation_cache[[cache_key]])
  }
  
  # For non-log scale output, need M for proper transformation
  if (!logScale && is.null(M)) {
    if (obs_type %in% c("M", "M_paired")) {
      stop("M must be provided when logScale = FALSE for observation type '", 
           obs_type, "'", call. = FALSE)
    }
  }
  
  # Validate M if provided
  if (!is.null(M)) {
    validate_M_identity(M, private$.M_fingerprint)
  }
  
  # Split by sample
  numSamples <- self$aux$numSamples
  Samples <- private$.samples
  sampleNames <- private$.sampleNames
  
  cells_by_sample <- split(seq_along(Samples), Samples)[sampleNames]
  
  # Compute imputed Y for each sample
  Y_imputed_list <- vector("list", numSamples)

  if (private$.verbose == 1) {
    cat(sprintf("Computing imputed Y (%d samples)\n", numSamples))
    cat("|", rep(" ", 50), "| 0%\r", sep = "")
    flush.console()

    for (i in 1:numSamples) {
      # Reconstruct Yi_fitted from parameters (NO M needed!)
      Yi_fitted <- compute_Y_fitted(self, private, i)

      # Remove sample effects
      Y_imputed_list[[i]] <- compute_Y_imputed(
        Yi_fitted = Yi_fitted,
        params = self$params,
        sample_idx = i,
        rowCentre = rowCentre
      )

      # Update progress bar
      pct <- round(i / numSamples * 100)
      n_filled <- round(50 * i / numSamples)
      cat("|", rep("=", n_filled), rep(" ", 50 - n_filled), "| ", pct, "%\r", sep = "")
      flush.console()
    }
    cat("\n")  # Final newline

  } else {
    # Silent or debug mode - no progress bar
    for (i in 1:numSamples) {
      # Reconstruct Yi_fitted from parameters (NO M needed!)
      Yi_fitted <- compute_Y_fitted(self, private, i)

      # Remove sample effects
      Y_imputed_list[[i]] <- compute_Y_imputed(
        Yi_fitted = Yi_fitted,
        params = self$params,
        sample_idx = i,
        rowCentre = rowCentre
      )
    }
  }

  # Concatenate all samples
  Y_imputed <- do.call(cbind, Y_imputed_list)

  # Add dimension names
  rownames(Y_imputed) <- private$.geneIDs
  colnames(Y_imputed) <- private$.cellIDs
  
  # Transform to non-log scale if requested
  if (!logScale) {
    if (obs_type == "M") {
      Y_imputed <- exp(Y_imputed)
      
    } else if (obs_type == "M_paired") {
      Y_imputed <- 1 / (1 + exp(-Y_imputed))
      
    } else if (obs_type == "Y") {
      Y_imputed <- exp(Y_imputed)
      
    } else {
      stop("Unrecognized observation type: ", obs_type, call. = FALSE)
    }
  }
  
  # Cache result
  if (is.null(private$.imputation_cache)) {
    private$.imputation_cache <- list()
  }
  private$.imputation_cache[[cache_key]] <- Y_imputed

  details <- format_dims(nrow(Y_imputed), ncol(Y_imputed), "genes", "cells")
  log_complete("Imputed Y", details, private$.verbose)
  log_info(sprintf("logScale = %s, rowCentre = %s", logScale, rowCentre), private$.verbose)
  
  return(Y_imputed)
}


#' Get Y Variance (with caching)
#'
#' @description
#' Returns posterior variance of imputed Y. Results are cached.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param M Count matrix (required - must match original)
#'
#' @return Dense matrix (J × N) with variance values, or NULL if not applicable
#'
#' @keywords internal
#' @noRd
get_Y_variance <- function(self, private, M) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if M provided
  if (missing(M) || is.null(M)) {
    stop("M must be provided for variance calculation", call. = FALSE)
  }
  
  # Check cache
  if (!is.null(private$.imputation_cache) &&
      !is.null(private$.imputation_cache$Y_var)) {
    log_cached("Y variance", private$.verbose)
    return(private$.imputation_cache$Y_var)
  }
  
  # Call core computation function
  Y_var <- compute_Y_variance(self, private, M)
  
  # Cache result
  if (!is.null(Y_var)) {
    if (is.null(private$.imputation_cache)) {
      private$.imputation_cache <- list()
    }
    private$.imputation_cache$Y_var <- Y_var

    details <- format_dims(nrow(Y_var), ncol(Y_var), "genes", "cells")
    log_complete("Y variance", details, private$.verbose)
  }
  
  return(Y_var)
}


#' Get Dispersion Analysis
#'
#' @description
#' Analyzes dispersion (variance vs mean relationship) for count data.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param M Count matrix (required - must match original)
#' @param subsample Integer, maximum number of nonzero positions to sample
#'
#' @return Data frame with dispersion statistics
#'
#' @keywords internal
#' @noRd
get_dispersion <- function(self, private, M, subsample = 1e6) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if M provided
  if (missing(M) || is.null(M)) {
    stop("M must be provided for dispersion analysis", call. = FALSE)
  }
  
  # Call core computation function
  dispersion_df <- compute_dispersion(self, private, M, subsample)

  details <- sprintf("%d bins across %d samples", nrow(dispersion_df), self$aux$numSamples)
  log_complete("Dispersion analysis", details, private$.verbose)

  return(dispersion_df)
}


#' Print Method for Imputation Accessor
#'
#' @param x Object of class gedi_imputation
#' @param ... Additional arguments (ignored)
#'
#' @keywords internal
#' @export
print.gedi_imputation <- function(x, ...) {
  cat("<GEDI Imputation Interface>\n")
  cat("\nAvailable methods:\n\n")
  
  cat("Imputed Expression:\n")
  cat("  $Y(M = NULL, logScale = TRUE, rowCentre = TRUE)\n")
  cat("    Get imputed gene expression with sample effects removed\n")
  cat("    - M: Optional count matrix (needed for non-log scale)\n")
  cat("    - logScale: Return log-transformed values (default: TRUE)\n")
  cat("    - rowCentre: Remove global gene offset (default: TRUE)\n\n")
  
  cat("Variance Analysis:\n")
  cat("  $variance(M)\n")
  cat("    Get posterior variance of imputed Y\n")
  cat("    - M: Count matrix (required, must match original)\n\n")
  
  cat("Dispersion Analysis:\n")
  cat("  $dispersion(M, subsample = 1e6)\n")
  cat("    Analyze variance vs mean relationship\n")
  cat("    - M: Count matrix (required, must match original)\n")
  cat("    - subsample: Max nonzero positions to sample per sample\n\n")
  
  cat("Example usage:\n")
  cat("  # Get imputed expression (no M needed!)\n")
  cat("  Y_imputed <- model$imputed$Y()\n")
  cat("  Y_raw_scale <- model$imputed$Y(M, logScale = FALSE)\n\n")
  
  cat("  # Variance and dispersion (M required)\n")
  cat("  Y_var <- model$imputed$variance(M)\n")
  cat("  disp <- model$imputed$dispersion(M)\n")
  
  invisible(x)
}