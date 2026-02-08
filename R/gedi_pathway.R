# ==============================================================================
# GEDI Pathway Analysis Functions
# Compute pathway-factor associations in dense and sparse forms
# ==============================================================================

#' Compute Dense Pathway-Factor Associations
#'
#' @description
#' Computes pathway-factor association matrix A in the original pathway space
#' by rotating the learned A matrix from the reduced SVD space back to the
#' original pathway annotations.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#'
#' @return Dense matrix (num_pathways × K) where rows correspond to original
#'   pathways and columns to latent factors. Positive values indicate pathways
#'   that are enriched in a factor, negative values indicate depletion.
#'
#' @details
#' The model internally learns A in a reduced P-dimensional space (where P is
#' determined by SVD of the C matrix). This function maps A back to the original
#' pathway space using the C.rotation matrix:
#' 
#' \deqn{A_{full} = C.rotation \times A_{reduced}}
#' 
#' This provides biological interpretability by showing how each original
#' pathway/gene set relates to each latent factor.
#'
#' @keywords internal
#' @noRd
compute_dense_A <- function(self, private) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if C prior was provided
  if (self$aux$P == 0) {
    stop("Cannot compute pathway associations: no gene-level prior (C) was provided during setup.",
         call. = FALSE)
  }
  
  # Check if aux_static exists
  if (is.null(private$.aux_static)) {
    stop("Missing auxiliary data. Model may not be properly initialized.",
         call. = FALSE)
  }
  
  # Check if inputC exists (need original pathway names)
  if (is.null(private$.aux_static$inputC)) {
    stop("Original C matrix not found. Cannot determine pathway names.",
         call. = FALSE)
  }
  
  if (private$.verbose >= 1) {
    message("Computing dense pathway-factor associations...")
  }
  
  # Compute A in original pathway space: C.rotation %*% A
  # C.rotation is num_pathways × P
  # A is P × K
  # Result is num_pathways × K
  A_full <- private$.aux_static$C.rotation %*% self$params$A
  
  # Add row and column names
  rownames(A_full) <- colnames(private$.aux_static$inputC)
  colnames(A_full) <- paste0("LV", seq_len(self$aux$K))
  
  if (private$.verbose >= 1) {
    message(sprintf("  Computed %d pathways x %d factors", 
                nrow(A_full), ncol(A_full)))
  }
  
  return(A_full)
}


#' Compute Sparse Pathway-Factor Associations
#'
#' @description
#' Re-fits pathway-factor associations using LASSO regression (L1 regularization)
#' to obtain sparse associations. This is useful for interpretation as it
#' identifies the most important pathways for each factor while setting many
#' coefficients to exactly zero.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param C Optional gene × pathway matrix. If NULL (default), uses the original
#'   C matrix provided during model setup. Can provide a different pathway
#'   database for post-hoc interpretation.
#'
#' @return Sparse matrix (num_pathways × K) with many zero entries. Non-zero
#'   values indicate pathways that are strongly associated with each factor.
#'
#' @details
#' For each latent factor k, this function:
#' \enumerate{
#'   \item Extracts the gene expression pattern Z[, k]
#'   \item Performs L1-regularized (LASSO) regression: Z[, k] ~ C
#'   \item Uses cross-validation to select optimal regularization strength
#'   \item Extracts non-zero pathway coefficients
#' }
#' 
#' The result is typically much sparser than the dense A from model training,
#' making it easier to identify key biological processes.
#' 
#' Requires the \code{glmnet} package to be installed.
#'
#' @examples
#' \dontrun{
#' # Use original pathways (sparse version)
#' A_sparse <- model$pathway_associations$sparse()
#' 
#' # Test with KEGG pathways post-hoc
#' kegg_C <- create_kegg_matrix(rownames(counts))
#' A_kegg <- model$pathway_associations$sparse(C = kegg_C)
#' }
#'
#' @keywords internal
#' @noRd
compute_sparse_A <- function(self, private, C = NULL) {
  
  # Check if glmnet is available
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for sparse pathway associations.\n",
         "  Install it with: install.packages('glmnet')",
         call. = FALSE)
  }
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Handle C matrix
  if (is.null(C)) {
    # Use original C if available
    if (self$aux$P == 0) {
      stop("No gene-level prior (C) was provided during setup.\n",
           "  You must provide a C matrix to compute sparse associations.",
           call. = FALSE)
    }
    
    if (is.null(private$.aux_static) || is.null(private$.aux_static$inputC)) {
      stop("Original C matrix not found. Cannot compute sparse associations.",
           call. = FALSE)
    }
    
    C <- private$.aux_static$inputC
    
    if (private$.verbose >= 1) {
      message("Computing sparse pathway-factor associations (using original C)...")
    }
    
  } else {
    # User provided custom C
    C <- as.matrix(C)
    
    # Validate dimensions
    if (nrow(C) != self$aux$J) {
      stop(sprintf("Custom C matrix must have J rows (got %d, expected %d)",
                   nrow(C), self$aux$J),
           call. = FALSE)
    }
    
    if (ncol(C) == 0) {
      stop("Custom C matrix must have at least one column (pathway)",
           call. = FALSE)
    }
    
    if (private$.verbose >= 1) {
      message(sprintf("Computing sparse pathway-factor associations (custom C: %d pathways)...",
                  ncol(C)))
    }
  }
  
  # Get Z matrix
  Z <- self$params$Z
  K <- self$aux$K
  num_pathways <- ncol(C)
  
  # Initialize result matrix
  A_sparse <- matrix(NA, nrow = num_pathways, ncol = K)
  
  # For each latent factor, run LASSO regression
  if (private$.verbose >= 1) {
    message(sprintf("  Running LASSO for %d factors...", K))
  }
  
  for (k in 1:K) {
    if (private$.verbose >= 2) {
      message(sprintf("    Factor %d/%d", k, K))
    }
    
    # Run cross-validated LASSO
    # Note: intercept = FALSE because Z columns should be centered
    tryCatch({
      fit <- glmnet::cv.glmnet(
        x = C,
        y = Z[, k],
        intercept = FALSE,
        standardize = TRUE,
        nfolds = 10
      )
      
      # Extract coefficients at lambda.1se (more regularized, sparser)
      # This gives more interpretable results than lambda.min
      coefs <- as.vector(coef(fit, s = "lambda.1se"))
      
      # Remove intercept (first element) even though we set intercept=FALSE
      # glmnet still returns it
      A_sparse[, k] <- coefs[-1]
      
    }, error = function(e) {
      warning(sprintf("LASSO failed for factor %d: %s", k, e$message),
              call. = FALSE)
      A_sparse[, k] <- rep(0, num_pathways)
    })
  }
  
  # Add row and column names
  if (!is.null(colnames(C))) {
    rownames(A_sparse) <- colnames(C)
  } else {
    rownames(A_sparse) <- paste0("Pathway", seq_len(num_pathways))
  }
  colnames(A_sparse) <- paste0("LV", seq_len(K))
  
  # Report sparsity
  if (private$.verbose >= 1) {
    nonzero_per_factor <- colSums(A_sparse != 0)
    total_nonzero <- sum(A_sparse != 0)
    total_possible <- num_pathways * K
    sparsity_pct <- 100 * (1 - total_nonzero / total_possible)
    
    message(sprintf("  Sparsity: %.1f%% (%d / %d nonzero)",
                sparsity_pct, total_nonzero, total_possible))
    message(sprintf("  Nonzero per factor: min=%d, median=%d, max=%d",
                min(nonzero_per_factor),
                median(nonzero_per_factor),
                max(nonzero_per_factor)))
  }
  
  return(A_sparse)
}


#' Create Pathway Associations Accessor Object
#'
#' @description
#' Creates a list-like object with methods to access pathway-factor associations
#' in both dense and sparse forms. This is used internally by the
#' pathway_associations active binding in the GEDI class.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#'
#' @return S3 object of class "gedi_pathway_associations" with two methods:
#'   \itemize{
#'     \item \code{$dense()} - Returns dense pathway associations
#'     \item \code{$sparse(C = NULL)} - Returns sparse pathway associations
#'   }
#'
#' @keywords internal
#' @noRd
create_pathway_associations_accessor <- function(self, private) {
  
  structure(
    list(
      dense = function() {
        compute_dense_A(self, private)
      },
      
      sparse = function(C = NULL) {
        compute_sparse_A(self, private, C)
      }
    ),
    class = "gedi_pathway_associations"
  )
}


#' Print Method for Pathway Associations Accessor
#'
#' @param x Object of class gedi_pathway_associations
#' @param ... Additional arguments (ignored)
#'
#' @keywords internal
#' @export
print.gedi_pathway_associations <- function(x, ...) {
  cat("<GEDI Pathway Associations>\n")
  cat("\nAvailable methods:\n")
  cat("  $dense()       - Compute dense pathway-factor associations\n")
  cat("                   (from model training, in original pathway space)\n")
  cat("\n")
  cat("  $sparse(C)     - Compute sparse pathway-factor associations\n")
  cat("                   (LASSO regression for interpretability)\n")
  cat("                   C: optional custom pathway matrix (default: uses original C)\n")
  cat("\nExample usage:\n")
  cat("  A_dense <- model$pathway_associations$dense()\n")
  cat("  A_sparse <- model$pathway_associations$sparse()\n")
  cat("  A_kegg <- model$pathway_associations$sparse(C = kegg_pathways)\n")
  
  invisible(x)
}
