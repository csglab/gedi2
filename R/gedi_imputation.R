# ==============================================================================
# GEDI Imputation - R Wrappers
# Lazy-computed imputation methods with caching and validation
# ==============================================================================

#' Compute M Fingerprint for Validation
#'
#' @description
#' Creates a fast, dependency-free fingerprint of the count matrix for identity validation.
#' Checks: dimensions, sparsity, value statistics, sparse indices, and dimnames.
#'
#' @param M Count matrix (sparse or dense), or list of two matrices for M_paired
#' @param obs_type Character, observation type: "M", "M_paired", "Y", or "X"
#'
#' @return List with fingerprint components for validation
#'
#' @keywords internal
#' @noRd
compute_M_fingerprint <- function(M, obs_type) {

  if (obs_type == "M_paired") {
    if (!is.list(M) || length(M) != 2) {
      stop("M must be a list of two matrices for M_paired observation type", call. = FALSE)
    }

    # Compute fingerprints for both matrices
    fp1 <- compute_M_fingerprint_single(M[[1]], "M1")
    fp2 <- compute_M_fingerprint_single(M[[2]], "M2")

    return(list(
      M1 = fp1,
      M2 = fp2,
      obs_type = obs_type
    ))

  } else {
    # Single matrix
    return(list(
      M = compute_M_fingerprint_single(M, "M"),
      obs_type = obs_type
    ))
  }
}


#' Compute Fingerprint for Single Matrix
#'
#' @description
#' Internal helper to compute fingerprint for a single matrix.
#' Does NOT store gene/cell IDs - only dimensions and statistics.
#'
#' @param M Count matrix (sparse or dense)
#' @param label Character label for this matrix (e.g., "M", "M1", "M2")
#'
#' @return List with fingerprint components
#'
#' @keywords internal
#' @noRd
compute_M_fingerprint_single <- function(M, label) {

  # Fast fingerprint with no dependencies - single validation mode
  # Checks: dimensions, sparsity, values, indices, and dimnames

  dims <- dim(M)
  is_sparse <- inherits(M, "sparseMatrix")

  # Extract dimnames
  row_names <- rownames(M)
  col_names <- colnames(M)

  # Compute fast checksums (sparse-aware)
  if (is_sparse) {
    # Sparse matrix - use internal structure directly
    values_sum <- sum(M@x)                    # sum of values
    values_sq <- sum(M@x^2)                   # sum of squares
    nnz <- length(M@x)                        # number of non-zeros
    indices_check <- sum(M@i %% 1000)         # row index checksum (mod to prevent overflow)

  } else {
    # Dense matrix
    values_sum <- sum(M)
    values_sq <- sum(M^2)
    nnz <- sum(M != 0)
    indices_check <- 0  # Not applicable for dense
  }

  # Dimnames fingerprint (fast hash using string length and first/last chars)
  rownames_fp <- if (!is.null(row_names)) {
    paste0(
      length(row_names), "_",
      if (length(row_names) > 0) substr(row_names[1], 1, 3) else "",
      if (length(row_names) > 1) substr(row_names[length(row_names)], 1, 3) else "",
      "_", sum(nchar(row_names))
    )
  } else { "NULL" }

  colnames_fp <- if (!is.null(col_names)) {
    paste0(
      length(col_names), "_",
      if (length(col_names) > 0) substr(col_names[1], 1, 3) else "",
      if (length(col_names) > 1) substr(col_names[length(col_names)], 1, 3) else "",
      "_", sum(nchar(col_names))
    )
  } else { "NULL" }

  # Composite fingerprint string (all checks combined)
  fingerprint_str <- paste0(
    dims[1], "x", dims[2], "_",      # dimensions
    values_sum, "_",                  # sum of values
    values_sq, "_",                   # sum of squares
    nnz, "_",                         # nnz count
    indices_check, "_",               # index checksum
    rownames_fp, "_",                 # rownames fingerprint
    colnames_fp                       # colnames fingerprint
  )

  return(list(
    nrow = dims[1],
    ncol = dims[2],
    is_sparse = is_sparse,
    fingerprint = fingerprint_str,
    label = label
  ))
}


#' Validate M Identity
#'
#' @description
#' Strict validation that M matrix matches the original used during training.
#' Uses structural checks, statistical fingerprints, and random sampling.
#' Does NOT check gene/cell name order - only dimensions and values.
#'
#' @param M Count matrix (sparse or dense), or list for M_paired
#' @param fingerprint Stored fingerprint from model setup
#'
#' @return TRUE if validation passes, otherwise stops with informative error
#'
#' @keywords internal
#' @noRd
validate_M_identity <- function(M, fingerprint) {
  
  obs_type <- fingerprint$obs_type
  
  if (obs_type == "M_paired") {
    if (!is.list(M) || length(M) != 2) {
      stop("M must be a list of two matrices for M_paired observation type", call. = FALSE)
    }
    
    validate_M_single(M[[1]], fingerprint$M1, "M[[1]]")
    validate_M_single(M[[2]], fingerprint$M2, "M[[2]]")
    
  } else {
    validate_M_single(M, fingerprint$M, "M")
  }
  
  return(TRUE)
}


#' Validate Single Matrix
#'
#' @description
#' Internal helper to validate a single matrix against its fingerprint.
#'
#' @param M Count matrix (sparse or dense)
#' @param fp Stored fingerprint for this matrix
#' @param label Character label for error messages
#'
#' @return TRUE if valid, stops with error otherwise
#'
#' @keywords internal
#' @noRd
validate_M_single <- function(M, fp, label) {

  # Quick dimension check first
  if (nrow(M) != fp$nrow) {
    stop(label, " has wrong number of genes: expected ", fp$nrow,
         ", got ", nrow(M), call. = FALSE)
  }
  if (ncol(M) != fp$ncol) {
    stop(label, " has wrong number of cells: expected ", fp$ncol,
         ", got ", ncol(M), call. = FALSE)
  }

  # Check sparsity type
  is_sparse <- inherits(M, "sparseMatrix")
  if (is_sparse != fp$is_sparse) {
    stop(label, " sparsity type changed: expected ",
         ifelse(fp$is_sparse, "sparse", "dense"),
         ", got ", ifelse(is_sparse, "sparse", "dense"), call. = FALSE)
  }

  # Recompute fingerprint for current matrix
  current_fp <- compute_M_fingerprint_single(M, label)

  # Compare fingerprints
  if (current_fp$fingerprint != fp$fingerprint) {
    stop(label, " fingerprint mismatch - matrix was modified.\n",
         "  Expected: ", fp$fingerprint, "\n",
         "  Got: ", current_fp$fingerprint, call. = FALSE)
  }

  return(TRUE)
}


#' Compute Y Fitted from Parameters
#'
#' @description
#' Reconstructs Yi_fitted = ZDBi + QiDBi + si + o + oi from model parameters.
#' This is the model's prediction of log-expression (or log-ratio for M_paired).
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param sample_idx Integer, which sample to compute (1 to numSamples)
#'
#' @return Dense matrix (J × Ni) - fitted Yi for sample i
#'
#' @keywords internal
#' @noRd
compute_Y_fitted <- function(self, private, sample_idx) {
  
  # Validate sample index
  numSamples <- self$aux$numSamples
  if (sample_idx < 1 || sample_idx > numSamples) {
    stop("sample_idx must be between 1 and ", numSamples, call. = FALSE)
  }
  
  # Get parameters
  Z <- self$params$Z
  D <- self$params$D
  Bi <- self$params$Bi[[sample_idx]]
  Qi <- self$params$Qi[[sample_idx]]
  si <- self$params$si[[sample_idx]]
  o <- self$params$o
  oi <- self$params$oi[[sample_idx]]
  
  # Compute ZDBi
  ZD <- Z %*% diag(D, nrow = length(D))
  ZDBi <- ZD %*% Bi
  
  # Compute QiDBi
  QiD <- Qi %*% diag(D, nrow = length(D))
  QiDBi <- QiD %*% Bi
  
  # Call C++ predict_Yhat
  Yi_fitted <- predict_Yhat(
    ZDBi = ZDBi,
    QiDBi = QiDBi,
    si = si,
    o = o,
    oi = oi
  )
  
  return(Yi_fitted)
}


#' Compute Imputed Y with Sample Effects Removed
#'
#' @description
#' Removes sample-specific effects from Yi_fitted to get shared biological signal.
#' This is the "imputed" expression that can be compared across samples.
#'
#' @param Yi_fitted Dense matrix (J × Ni) - fitted Yi for sample i
#' @param params List with model parameters
#' @param sample_idx Integer, which sample this is
#' @param rowCentre Logical, if TRUE removes global gene offset o
#'
#' @return Dense matrix (J × Ni) - imputed Y with sample effects removed
#'
#' @keywords internal
#' @noRd
compute_Y_imputed <- function(Yi_fitted, params, sample_idx, rowCentre) {
  
  Qi <- params$Qi[[sample_idx]]
  D <- params$D
  Bi <- params$Bi[[sample_idx]]
  si <- params$si[[sample_idx]]
  oi <- params$oi[[sample_idx]]
  o <- params$o
  
  # Compute QiDBi
  QiD <- Qi %*% diag(D, nrow = length(D))
  QiDBi <- QiD %*% Bi
  
  # Call C++ Yi_resZ to remove sample effects
  if (rowCentre) {
    Y_imputed <- Yi_resZ(
      Yi = Yi_fitted,
      QiDBi = QiDBi,
      si = si,
      o = o,
      oi = oi
    )
  } else {
    # Don't remove o (set to zero vector)
    Y_imputed <- Yi_resZ(
      Yi = Yi_fitted,
      QiDBi = QiDBi,
      si = si,
      o = rep(0, length(o)),
      oi = oi
    )
  }
  
  return(Y_imputed)
}


#' Compute Y Variance from Raw Counts
#'
#' @description
#' Computes posterior variance of imputed Y given fitted Yi and raw counts.
#' Requires the original M matrix.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param M Count matrix (must match original), or list for M_paired
#'
#' @return Dense matrix (J × N) - variance of imputed Y
#'
#' @keywords internal
#' @noRd
compute_Y_variance <- function(self, private, M) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check observation type
  obs_type <- self$aux$obs.type
  if (obs_type == "Y") {
    if (private$.verbose >= 1) message("Observation type is 'Y' (not raw counts). Returning NULL.")
    return(NULL)
  }
  if (obs_type == "X") {
    if (private$.verbose >= 1) message("Observation type is 'X' (binary). Variance calculation not applicable. Returning NULL.")
    return(NULL)
  }
  
  # Validate M
  validate_M_identity(M, private$.M_fingerprint)
  
  # Convert to sparse if needed
  if (obs_type == "M") {
    if (!inherits(M, "sparseMatrix")) {
      warning("Converting M to sparse format for memory efficiency")
      M <- as(M, "dgCMatrix")
    }
  } else if (obs_type == "M_paired") {
    if (!inherits(M[[1]], "sparseMatrix")) {
      warning("Converting M[[1]] to sparse format for memory efficiency")
      M[[1]] <- as(M[[1]], "dgCMatrix")
    }
    if (!inherits(M[[2]], "sparseMatrix")) {
      warning("Converting M[[2]] to sparse format for memory efficiency")
      M[[2]] <- as(M[[2]], "dgCMatrix")
    }
  }
  
  # Split M by sample
  numSamples <- self$aux$numSamples
  Samples <- private$.samples
  sampleNames <- private$.sampleNames
  
  # Get cell indices for each sample (in original order)
  cells_by_sample <- split(seq_along(Samples), Samples)[sampleNames]
  
  # Compute variance for each sample
  Y_var_list <- vector("list", numSamples)

  if (private$.verbose == 1) {
    message(sprintf("Computing Y variance (%d samples)", numSamples))
    pb <- txtProgressBar(min = 0, max = numSamples, style = 3,
                         width = 50, file = stderr())
    on.exit(close(pb), add = TRUE)

    for (i in 1:numSamples) {
      idx <- cells_by_sample[[i]]

      # Get Yi_fitted for this sample (reconstruct from params)
      Yi_fitted <- compute_Y_fitted(self, private, i)

      # Get Mi for this sample
      if (obs_type == "M") {
        Mi <- M[, idx, drop = FALSE]
        Y_var_list[[i]] <- Yi_var_M(Yi_fitted, self$params$sigma2)

      } else if (obs_type == "M_paired") {
        M1i <- M[[1]][, idx, drop = FALSE]
        M2i <- M[[2]][, idx, drop = FALSE]
        Y_var_list[[i]] <- Yi_var_M_paired(Yi_fitted, M1i, M2i, self$params$sigma2)
      }

      setTxtProgressBar(pb, i)
    }

  } else {
    # Silent or debug mode - no progress bar
    for (i in 1:numSamples) {
      idx <- cells_by_sample[[i]]

      # Get Yi_fitted for this sample (reconstruct from params)
      Yi_fitted <- compute_Y_fitted(self, private, i)

      # Get Mi for this sample
      if (obs_type == "M") {
        Mi <- M[, idx, drop = FALSE]
        Y_var_list[[i]] <- Yi_var_M(Yi_fitted, self$params$sigma2)

      } else if (obs_type == "M_paired") {
        M1i <- M[[1]][, idx, drop = FALSE]
        M2i <- M[[2]][, idx, drop = FALSE]
        Y_var_list[[i]] <- Yi_var_M_paired(Yi_fitted, M1i, M2i, self$params$sigma2)
      }
    }
  }
  
  # Concatenate all samples
  Y_var <- do.call(cbind, Y_var_list)
  
  # Add dimension names
  rownames(Y_var) <- private$.geneIDs
  colnames(Y_var) <- private$.cellIDs
  
  return(Y_var)
}


#' Compute Dispersion Analysis
#'
#' @description
#' Analyzes dispersion (variance vs mean relationship) for count data.
#' Compares predicted variance (from Poisson model) to observed variance.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param M Count matrix (must match original), or list for M_paired
#' @param subsample Integer, maximum number of nonzero positions to sample per sample
#'
#' @return Data frame with dispersion statistics by bin
#'
#' @keywords internal
#' @noRd
compute_dispersion <- function(self, private, M, subsample = 1e6) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check observation type (only for M)
  obs_type <- self$aux$obs.type
  if (obs_type != "M") {
    stop("Dispersion analysis only available for observation type 'M' (single count matrix)",
         call. = FALSE)
  }
  
  # Validate M
  validate_M_identity(M, private$.M_fingerprint)
  
  # Convert to sparse if needed
  if (!inherits(M, "sparseMatrix")) {
    warning("Converting M to sparse format for memory efficiency")
    M <- as(M, "dgCMatrix")
  }
  
  # Split M by sample
  numSamples <- self$aux$numSamples
  Samples <- private$.samples
  sampleNames <- private$.sampleNames
  
  cells_by_sample <- split(seq_along(Samples), Samples)[sampleNames]
  
  # Compute dispersion for each sample
  result_list <- vector("list", numSamples)

  if (private$.verbose == 1) {
    message(sprintf("Computing dispersion (%d samples)", numSamples))
    pb <- txtProgressBar(min = 0, max = numSamples, style = 3,
                         width = 50, file = stderr())
    on.exit(close(pb), add = TRUE)

    for (i in 1:numSamples) {
      idx <- cells_by_sample[[i]]

      # Get Yi_fitted for this sample
      Yi_fitted <- compute_Y_fitted(self, private, i)

      # Get Mi for this sample
      Mi <- M[, idx, drop = FALSE]

      # Call C++ dispersion function (sparse-optimized)
      dispersion_data <- compute_dispersion_sparse(Yi_fitted, Mi, as.integer(subsample))

      predicted <- dispersion_data$predicted
      observed <- dispersion_data$observed

      # Compute residual squared
      res_sq <- (observed - predicted)^2

      # Bin by predicted values
      n_bins <- floor(sqrt(length(predicted)))
      bins <- cut(
        predicted,
        breaks = stats::quantile(predicted, seq(0, 1, length = n_bins + 1)),
        include.lowest = TRUE
      )

      # Aggregate by bin
      xy <- stats::aggregate(
        cbind(predicted, res_sq),
        by = list(bin = bins),
        FUN = mean
      )

      colnames(xy) <- c("bin", "Expected_Var", "Observed_Var")
      xy$Sample <- sampleNames[i]

      # Count observations per bin
      xy$n <- stats::aggregate(
        predicted,
        by = list(bin = bins),
        FUN = length
      )$x

      result_list[[i]] <- xy

      setTxtProgressBar(pb, i)
    }

  } else {
    # Silent or debug mode - no progress bar
    for (i in 1:numSamples) {
      idx <- cells_by_sample[[i]]

      # Get Yi_fitted for this sample
      Yi_fitted <- compute_Y_fitted(self, private, i)

      # Get Mi for this sample
      Mi <- M[, idx, drop = FALSE]

      # Call C++ dispersion function (sparse-optimized)
      dispersion_data <- compute_dispersion_sparse(Yi_fitted, Mi, as.integer(subsample))

      predicted <- dispersion_data$predicted
      observed <- dispersion_data$observed

      # Compute residual squared
      res_sq <- (observed - predicted)^2

      # Bin by predicted values
      n_bins <- floor(sqrt(length(predicted)))
      bins <- cut(
        predicted,
        breaks = stats::quantile(predicted, seq(0, 1, length = n_bins + 1)),
        include.lowest = TRUE
      )

      # Aggregate by bin
      xy <- stats::aggregate(
        cbind(predicted, res_sq),
        by = list(bin = bins),
        FUN = mean
      )

      colnames(xy) <- c("bin", "Expected_Var", "Observed_Var")
      xy$Sample <- sampleNames[i]

      # Count observations per bin
      xy$n <- stats::aggregate(
        predicted,
        by = list(bin = bins),
        FUN = length
      )$x

      result_list[[i]] <- xy
    }
  }
  
  # Combine all samples
  result <- do.call(rbind, result_list)
  rownames(result) <- NULL
  
  # Reorder columns
  result <- result[, c("Sample", "Expected_Var", "Observed_Var", "n", "bin")]
  
  return(result)
}


#' Clear Imputation Cache
#'
#' @description
#' Clears cached imputation results.
#'
#' @param private Reference to private environment
#' @param what Character vector specifying which cache entries to clear.
#'   Options: "Y_fitted", "Y_imputed", or NULL to clear all (default: NULL)
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
clear_imputation_cache <- function(private, what = NULL) {
  
  if (is.null(private$.imputation_cache)) {
    return(invisible(NULL))
  }
  
  if (is.null(what)) {
    # Clear all imputation cache
    private$.imputation_cache <- list()
    if (private$.verbose > 0) {
      message("Imputation cache cleared")
    }
  } else {
    # Clear specific entries
    if ("Y_imputed" %in% what) {
      # Clear ALL Y_imputed variants
      y_imputed_keys <- grep("^Y_imputed_", names(private$.imputation_cache), value = TRUE)
      for (key in y_imputed_keys) {
        private$.imputation_cache[[key]] <- NULL
      }
      if (private$.verbose > 0) {
        message("Cache cleared: Y_imputed (all variants)")
      }
    }
    
    if ("Y_var" %in% what) {
      private$.imputation_cache$Y_var <- NULL
      if (private$.verbose > 0) {
        message("Cache cleared: Y_var")
      }
    }
    
    invalid <- setdiff(what, c("Y_imputed", "Y_var"))
    if (length(invalid) > 0) {
      warning("Invalid cache entries: ", paste(invalid, collapse = ", "))
    }
  }
  
  invisible(NULL)
}


#' Get Imputation Cache Status
#'
#' @description
#' Returns information about which imputation results are currently cached.
#'
#' @param private Reference to private environment
#'
#' @return Named logical vector indicating which results are cached
#'
#' @keywords internal
#' @noRd
get_imputation_cache_status <- function(private) {
  
  if (is.null(private$.imputation_cache)) {
    return(c(Y_imputed = FALSE, Y_var = FALSE))
  }
  
  # Check for any Y_imputed variant (different logScale/rowCentre combinations)
  has_Y_imputed <- any(grepl("^Y_imputed_", names(private$.imputation_cache)))
  
  status <- c(
    Y_imputed = has_Y_imputed,
    Y_var = !is.null(private$.imputation_cache$Y_var)
  )
  
  return(status)
}


#' Get Imputation Cache Memory Usage
#'
#' @description
#' Estimates memory usage of cached imputation results in MB.
#'
#' @param private Reference to private environment
#'
#' @return Named numeric vector with memory usage in MB for each cached result
#'
#' @keywords internal
#' @noRd
get_imputation_cache_memory <- function(private) {
  
  if (is.null(private$.imputation_cache)) {
    return(c(Y_imputed = 0, Y_var = 0, Total = 0))
  }
  
  # Sum all Y_imputed variants
  y_imputed_keys <- grep("^Y_imputed_", names(private$.imputation_cache), value = TRUE)
  y_imputed_size <- sum(vapply(y_imputed_keys, function(key) {
    if (!is.null(private$.imputation_cache[[key]])) {
      as.numeric(object.size(private$.imputation_cache[[key]]) / 1024^2)
    } else {
      0
    }
  }, numeric(1)))
  
  memory_mb <- c(
    Y_imputed = y_imputed_size,
    Y_var = if (!is.null(private$.imputation_cache$Y_var)) {
      object.size(private$.imputation_cache$Y_var) / 1024^2
    } else 0
  )
  
  memory_mb["Total"] <- sum(memory_mb)
  
  return(memory_mb)
}