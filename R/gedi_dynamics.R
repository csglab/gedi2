# ==============================================================================
# GEDI Dynamics Analysis
# Vector fields, activity gradients, and trajectory analysis
# ==============================================================================

#' Compute Differential Q in Z-space (JxK)
#'
#' @description
#' Computes sample-variable effects on Qi, returning a J × K matrix.
#' This is the R wrapper for the C++ function getDiffQ_cpp().
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param contrast Contrast vector (length L)
#'
#' @return Matrix J × K showing differential effect in Z-space
#'
#' @keywords internal
#' @noRd
get_diffQ_Z_space <- function(self, private, contrast) {
  
  # Validate contrast
  if (!is.numeric(contrast) || length(contrast) != self$aux$L) {
    stop("contrast must be a numeric vector of length L (", self$aux$L, ")", call. = FALSE)
  }
  
  # Call C++ function (implements the old getDiffQ.gedi logic)
  R <- getDiffQ_cpp(
    Rk_list = self$params$Rk,
    H_rotation = private$.aux_static$H.rotation,
    contrast = contrast,
    verbose = private$.verbose
  )
  
  # Set dimension names
  rownames(R) <- private$.geneIDs
  colnames(R) <- paste0("LV", 1:self$aux$K)
  
  return(R)
}


#' Compute Vector Field SVD (Internal Helper)
#'
#' @description
#' Internal helper to calculate SVD for cells and their vector fields.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param start.cond Contrast vector for start condition (length L)
#' @param end.cond Contrast vector for end condition (length L)
#' @param scale_cond_vector Scaling factor for vector field magnitude (default: 1)
#'
#' @return List with SVD results and indices
#'
#' @keywords internal
#' @noRd
run_vector_field_svd <- function(self, private, start.cond, end.cond, scale_cond_vector = 1) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if H prior was provided
  if (self$aux$L == 0) {
    stop("Cannot compute vector field: no sample-level prior (H) was provided during setup.",
         call. = FALSE)
  }
  
  # Validate contrasts
  if (!is.numeric(start.cond) || length(start.cond) != self$aux$L) {
    stop("start.cond must be a numeric vector of length L (", self$aux$L, ")", call. = FALSE)
  }
  if (!is.numeric(end.cond) || length(end.cond) != self$aux$L) {
    stop("end.cond must be a numeric vector of length L (", self$aux$L, ")", call. = FALSE)
  }
  
  # Get Z matrix
  Z <- self$params$Z
  K <- self$aux$K
  
  # Calculate Z-space differential expression (now calls C++)
  startQ_Z <- get_diffQ_Z_space(self, private, start.cond)
  endQ_Z <- get_diffQ_Z_space(self, private, 
                                     start.cond + (end.cond - start.cond) * scale_cond_vector)
  
  # Calculate coordinate vectors for start and end conditions
  startQ <- startQ_Z + Z
  endQ <- endQ_Z + Z
  
  # Map vectors onto Z coordinate system
  ZtZ <- crossprod(Z)
  rotStartQ <- solve(ZtZ, crossprod(Z, startQ))
  rotEndQ <- solve(ZtZ, crossprod(Z, endQ))
  
  # Get DB projection (uses cache)
  DB <- self$projections$DB
  
  # Concatenate start and end positions: [rotStartQ %*% DB, rotEndQ %*% DB]
  projDB <- cbind(rotStartQ %*% DB, rotEndQ %*% DB)
  
  # Perform factorized SVD (now calls C++)
  result <- run_factorized_svd_cpp(Z, projDB, private$.verbose)
  
  # Set dimension names
  names(result$d) <- paste0("LV", 1:K)
  colnames(result$u) <- paste0("LV", 1:K)
  colnames(result$v) <- paste0("LV", 1:K)
  rownames(result$u) <- private$.geneIDs
  rownames(result$v) <- c(
    paste0(private$.cellIDs, ".start"),
    paste0(private$.cellIDs, ".end")
  )
  
  # Return structured result
  structure(
    list(
      d = result$d,
      u = result$u,
      v = result$v,
      indices = list(
        embedding = 1:self$aux$N,
        vector_field = 1:(2 * self$aux$N)
      ),
      metadata = list(
        method = "vector_field",
        start_condition = start.cond,
        end_condition = end.cond,
        scale_factor = scale_cond_vector
      )
    ),
    class = "gedi_dynamics_svd"
  )
}


#' Compute Activity Gradient SVD (Internal Helper)
#'
#' @description
#' Internal helper to calculate SVD showing pathway activity gradient direction.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param C_index Index of pathway in C matrix
#' @param scale_gradient Scaling factor for gradient magnitude (NA = auto-compute)
#'
#' @return List with SVD results and indices
#'
#' @keywords internal
#' @noRd
run_gradient_svd <- function(self, private, C_index, scale_gradient = NA) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if C prior was provided
  if (self$aux$P == 0) {
    stop("Cannot compute activity gradient: no gene-level prior (C) was provided during setup.",
         call. = FALSE)
  }
  
  # Validate C_index
  if (!is.numeric(C_index) || length(C_index) != 1) {
    stop("C_index must be a single integer", call. = FALSE)
  }
  
  num_pathways <- ncol(private$.aux_static$inputC)
  if (C_index < 1 || C_index > num_pathways) {
    stop(sprintf("C_index must be between 1 and %d", num_pathways), call. = FALSE)
  }
  
  # Get matrices
  Z <- self$params$Z
  K <- self$aux$K
  
  # Compute A in original pathway space and extract gradient
  A_full <- private$.aux_static$C.rotation %*% self$params$A
  gradient <- A_full[C_index, ]
  
  # Get DB projection (uses cache)
  DB <- self$projections$DB
  
  # Auto-compute scale if needed
  if (is.na(scale_gradient)) {
    scale_gradient <- sqrt(var(as.vector(DB)) / mean(gradient^2)) * 0.2
  }
  
  # Concatenate: [DB, DB + gradient]
  projDB <- cbind(DB, DB + gradient * scale_gradient)
  
  # Perform factorized SVD (now calls C++)
  result <- run_factorized_svd_cpp(Z, projDB, private$.verbose)
  
  # Set dimension names
  pathway_name <- colnames(private$.aux_static$inputC)[C_index]
  names(result$d) <- paste0("LV", 1:K)
  colnames(result$u) <- paste0("LV", 1:K)
  colnames(result$v) <- paste0("LV", 1:K)
  rownames(result$u) <- private$.geneIDs
  rownames(result$v) <- c(
    paste0(private$.cellIDs, ".start"),
    paste0(private$.cellIDs, ".gradEnd")
  )
  
  # Return structured result
  structure(
    list(
      d = result$d,
      u = result$u,
      v = result$v,
      indices = list(
        embedding = 1:self$aux$N,
        gradient = 1:(2 * self$aux$N)
      ),
      metadata = list(
        method = "activity_gradient",
        pathway_index = C_index,
        pathway_name = pathway_name,
        scale_factor = scale_gradient
      )
    ),
    class = "gedi_dynamics_svd"
  )
}


#' Compute Joint SVD (Internal Helper)
#'
#' @description
#' Internal helper to create embedding for cells with both condition-based
#' vector field and pathway activity gradient.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#' @param start.cond Contrast vector for start condition
#' @param end.cond Contrast vector for end condition
#' @param C_index Index of pathway in C matrix
#' @param scale_cond_vector Scaling for condition vector field (default: 1)
#' @param scale_gradient Scaling for gradient (NA = auto-compute)
#'
#' @return List with SVD results and indices
#'
#' @keywords internal
#' @noRd
run_joint_svd <- function(self, private, start.cond, end.cond, C_index,
                               scale_cond_vector = 1, scale_gradient = NA) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if both H and C priors were provided
  if (self$aux$L == 0) {
    stop("Cannot compute joint analysis: no sample-level prior (H) was provided during setup.",
         call. = FALSE)
  }
  if (self$aux$P == 0) {
    stop("Cannot compute joint analysis: no gene-level prior (C) was provided during setup.",
         call. = FALSE)
  }
  
  # Validate contrasts
  if (!is.numeric(start.cond) || length(start.cond) != self$aux$L) {
    stop("start.cond must be a numeric vector of length L (", self$aux$L, ")", call. = FALSE)
  }
  if (!is.numeric(end.cond) || length(end.cond) != self$aux$L) {
    stop("end.cond must be a numeric vector of length L (", self$aux$L, ")", call. = FALSE)
  }
  
  # Validate C_index
  num_pathways <- ncol(private$.aux_static$inputC)
  if (C_index < 1 || C_index > num_pathways) {
    stop(sprintf("C_index must be between 1 and %d", num_pathways), call. = FALSE)
  }
  
  # Get matrices
  Z <- self$params$Z
  K <- self$aux$K
  
  # Compute A and extract gradient
  A_full <- private$.aux_static$C.rotation %*% self$params$A
  gradient <- A_full[C_index, ]
  
  # Calculate Z-space differential expression (now calls C++)
  startQ_Z <- get_diffQ_Z_space(self, private, start.cond)
  endQ_Z <- get_diffQ_Z_space(self, private, 
                                     start.cond + (end.cond - start.cond) * scale_cond_vector)
  
  # Calculate coordinate vectors for start and end conditions
  startQ <- startQ_Z + Z
  endQ <- endQ_Z + Z
  
  # Map vectors onto Z coordinate system
  ZtZ <- crossprod(Z)
  rotStartQ <- solve(ZtZ, crossprod(Z, startQ))
  rotEndQ <- solve(ZtZ, crossprod(Z, endQ))
  
  # Get DB projection (uses cache)
  DB <- self$projections$DB
  
  # Auto-compute gradient scale if needed
  if (is.na(scale_gradient)) {
    scale_gradient <- sqrt(var(as.vector(DB)) / mean(gradient^2)) * 0.2
  }
  
  # Concatenate: [rotStartQ*DB, rotEndQ*DB, rotStartQ*(DB+gradient)]
  projDB <- cbind(
    rotStartQ %*% DB,
    rotEndQ %*% DB,
    rotStartQ %*% (DB + gradient * scale_gradient)
  )
  
  # Perform factorized SVD (now calls C++)
  result <- run_factorized_svd_cpp(Z, projDB, private$.verbose)
  
  # Set dimension names
  pathway_name <- colnames(private$.aux_static$inputC)[C_index]
  names(result$d) <- paste0("LV", 1:K)
  colnames(result$u) <- paste0("LV", 1:K)
  colnames(result$v) <- paste0("LV", 1:K)
  rownames(result$u) <- private$.geneIDs
  rownames(result$v) <- c(
    paste0(private$.cellIDs, ".start"),
    paste0(private$.cellIDs, ".condEnd"),
    paste0(private$.cellIDs, ".gradEnd")
  )
  
  # Return structured result
  structure(
    list(
      d = result$d,
      u = result$u,
      v = result$v,
      indices = list(
        embedding = 1:self$aux$N,
        vector_field = 1:(2 * self$aux$N),
        gradient = c(1:self$aux$N, (2 * self$aux$N + 1):(3 * self$aux$N))
      ),
      metadata = list(
        method = "joint_vector_field_gradient",
        start_condition = start.cond,
        end_condition = end.cond,
        pathway_index = C_index,
        pathway_name = pathway_name,
        scale_cond_vector = scale_cond_vector,
        scale_gradient = scale_gradient
      )
    ),
    class = "gedi_dynamics_svd"
  )
}


#' Compute Pathway Gradients in Gene Space (Internal Helper)
#'
#' @description
#' Returns the gradient of all pathways projected into gene expression space.
#' This logic remains in R as it is just two fast matrix multiplications.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#'
#' @return Matrix (J × num_pathways) of gene expression gradients
#'
#' @keywords internal
#' @noRd
get_pathway_gradients <- function(self, private) {
  
  # Validate model state
  if (is.null(private$.lastResult)) {
    stop("No results yet. Run $train() first.", call. = FALSE)
  }
  
  # Check if C prior was provided
  if (self$aux$P == 0) {
    stop("Cannot compute pathway gradients: no gene-level prior (C) was provided during setup.",
         call. = FALSE)
  }
  
  # Compute A in original pathway space
  A_full <- private$.aux_static$C.rotation %*% self$params$A
  
  # Gradient is transpose of A
  gradient <- t(A_full)  # K × num_pathways
  
  # Project to gene space
  Z <- self$params$Z
  gradient <- Z %*% gradient  # J × num_pathways
  
  # Set dimension names
  rownames(gradient) <- private$.geneIDs
  colnames(gradient) <- colnames(private$.aux_static$inputC)
  
  return(gradient)
}


#' Create Dynamics Accessor Object
#'
#' @description
#' Creates accessor object for dynamics analysis with vector field and gradient methods.
#'
#' @param self Reference to GEDI R6 object
#' @param private Reference to private environment
#'
#' @return S3 object of class "gedi_dynamics"
#'
#' @keywords internal
#'@noRd
create_dynamics_accessor <- function(self, private) {
  
  structure(
    list(
      vector_field = function(start.cond, end.cond, scale = 1) {
        run_vector_field_svd(self, private, start.cond, end.cond, scale)
      },
      
      activity_gradient = function(C_index, scale = NA) {
        run_gradient_svd(self, private, C_index, scale)
      },
      
      joint = function(start.cond, end.cond, C_index, scale_cond = 1, scale_grad = NA) {
        run_joint_svd(self, private, start.cond, end.cond, C_index, scale_cond, scale_grad)
      },
      
      pathway_gradients = function() {
        get_pathway_gradients(self, private)
      }
    ),
    class = "gedi_dynamics"
  )
}


#' Print Method for Dynamics Accessor
#'
#' @param x Object of class gedi_dynamics
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns \code{x}.
#' @keywords internal
#' @export
print.gedi_dynamics <- function(x, ...) {
  cat("<GEDI Dynamics Analysis>\n")
  cat("\nAvailable methods:\n\n")
  
  cat("Vector Field Analysis:\n")
  cat("  $vector_field(start.cond, end.cond, scale = 1)\n")
  cat("    Compute SVD showing cell state transitions between conditions\n\n")
  
  cat("Pathway Activity:\n")
  cat("  $activity_gradient(C_index, scale = NA)\n")
  cat("    Compute SVD showing pathway activity gradient direction\n\n")
  
  cat("  $pathway_gradients()\n")
  cat("    Get all pathway gradients in gene expression space\n\n")
  
  cat("Joint Analysis:\n")
  cat("  $joint(start.cond, end.cond, C_index, scale_cond = 1, scale_grad = NA)\n")
  cat("    Combine vector field and gradient in single embedding\n\n")
  
  cat("Example usage:\n")
  cat("  vf <- model$dynamics$vector_field(c(0,0), c(1,0))\n")
  cat("  ag <- model$dynamics$activity_gradient(C_index = 5)\n")
  cat("  grads <- model$dynamics$pathway_gradients()\n")
  
  invisible(x)
}


#' Print Method for Dynamics SVD Results
#'
#' @param x Object of class gedi_dynamics_svd
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns \code{x}.
#' @keywords internal
#' @export
print.gedi_dynamics_svd <- function(x, ...) {
  cat("<GEDI Dynamics SVD Result>\n\n")
  cat(sprintf("Method: %s\n", x$metadata$method))
  cat(sprintf("Dimensions: %d genes × %d components\n", nrow(x$u), ncol(x$u)))
  cat(sprintf("Cells: %d\n", length(x$indices$embedding)))
  
  if (x$metadata$method == "vector_field") {
    cat(sprintf("\nVector field indices: %d:%d\n", 
                min(x$indices$vector_field), max(x$indices$vector_field)))
  } else if (x$metadata$method == "activity_gradient") {
    cat(sprintf("\nPathway: %s\n", x$metadata$pathway_name))
    cat(sprintf("Gradient scale: %.3f\n", x$metadata$scale_factor))
  } else if (x$metadata$method == "joint_vector_field_gradient") {
    cat(sprintf("\nPathway: %s\n", x$metadata$pathway_name))
    cat(sprintf("Vector field indices: %d:%d\n",
                min(x$indices$vector_field), max(x$indices$vector_field)))
    cat(sprintf("Gradient indices: %s\n", 
                paste(range(x$indices$gradient), collapse = ":")))
  }
  
  cat("\nSingular values:\n")
  print(head(x$d, 5))
  
  cat("\nComponents: $d, $u, $v, $indices, $metadata\n")
  
  invisible(x)
}