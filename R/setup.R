# ==============================================================================
# GEDI v2: Gene Expression Data Integration  
# Memory-efficient R6 wrapper around C++ core
# OPTION 2: C++ computes Yi from M (zero Yi copies in R)
# ==============================================================================

library(R6)
library(Matrix)

newLogger <- function(name, level = 1) {
  structure(
    list(
      name = name,
      level = level,
      info = function(msg) {
        if (level >= 1) message("[", name, "] ", msg)
      },
      debug = function(msg) {
        if (level >= 2) message("[", name, "] DEBUG: ", msg)
      },
      warn = function(msg) {
        warning("[", name, "] WARNING: ", msg, call. = FALSE)
      },
      error = function(msg) {
        stop("[", name, "] ERROR: ", msg, call. = FALSE)
      }
    ),
    class = "gedi_logger"
  )
}

generateRsvdO <- function(m, n, k, p = 10, sdist = "normal") {
  if (m < n) {
    temp <- m
    m <- n
    n <- temp
  }
  
  l <- min(k + p, n)
  
  O <- switch(sdist,
              normal = matrix(stats::rnorm(l * n), n, l),
              unif = matrix(stats::runif(l * n), n, l),
              rademacher = matrix(sample(c(-1, 1), l * n, replace = TRUE), n, l),
              stop("sdist must be: normal, unif, or rademacher", call. = FALSE))
  
  return(O)
}

GEDI <- R6Class(
  classname = "GEDI",
  
  private = list(
    .cppPtr = NULL,
    .lastResult = NULL,
    .verbose = 1,
    .logger = NULL,
    .isInitialized = FALSE,
    
    .geneIDs = NULL,
    .cellIDs = NULL,
    .sampleNames = NULL,
    
    prepareData = function(
      Samples, Y, X, M, colData, C, H, K, mode, adjustD, orthoZ,
      Z_shrinkage, A_shrinkage, Qi_shrinkage, Rk_shrinkage,
      oi_shrinkage, o_shrinkage, si_shrinkage, fixed_si,
      rsvd_p, rsvd_sdist
    ) {
      
      if (is.null(M) && is.null(Y) && is.null(X)) {
        private$.logger$error("One of Y, M, or X must be provided")
      }
      
      # For computing hyperparameters, we need Y_log_dense temporarily
      # But we won't pass it to C++ for M/M_paired/X cases
      Y_log_dense <- NULL
      obs_type <- NULL
      geneIDs <- NULL
      cellIDs <- NULL
      
      # ===================================================================
      # STEP 1: Determine observation type
      # ===================================================================
      
      if (!is.null(M)) {
        if (is.list(M)) {
          if (length(M) != 2) {
            private$.logger$error("M list must contain exactly 2 matrices")
          }
          if (!identical(dim(M[[1]]), dim(M[[2]]))) {
            private$.logger$error("M[[1]] and M[[2]] must have identical dimensions")
          }
          
          M1_dense <- as.matrix(M[[1]])
          M2_dense <- as.matrix(M[[2]])
          Y_log_dense <- log((M1_dense + 1) / (M2_dense + 1))
          
          obs_type <- "M_paired"
          geneIDs <- rownames(M[[1]])
          cellIDs <- colnames(M[[1]])
        } else {
          Y_log_dense <- as.matrix(log1p(M))
          
          obs_type <- "M"
          geneIDs <- rownames(M)
          cellIDs <- colnames(M)
        }
      } else if (!is.null(X)) {
        Y_log_dense <- as.matrix(X)
        Y_log_dense[which(X == 1)] <- 1
        Y_log_dense[which(X == 0)] <- -1
        Y_log_dense[which(is.na(X))] <- 0
        
        obs_type <- "X"
        geneIDs <- rownames(X)
        cellIDs <- colnames(X)
      } else {
        Y_log_dense <- as.matrix(Y)
        
        obs_type <- "Y"
        geneIDs <- rownames(Y)
        cellIDs <- colnames(Y)
      }
      
      # ===================================================================
      # STEP 2: Store dimensions
      # ===================================================================
      
      J <- nrow(Y_log_dense)
      J_vec <- rep(1, J)
      N <- ncol(Y_log_dense)
      N_vec <- rep(1, N)
      
      if (N != length(Samples)) {
        private$.logger$error("Number of cells doesn't match Samples length")
      }
      
      # ===================================================================
      # STEP 3: Process prior matrix C (gene-level)
      # ===================================================================
      
      inputC <- NULL
      C_matrix <- matrix(nrow = 0, ncol = 0)
      C_rotation <- NULL
      P <- 0
      
      if (!is.null(C)) {
        inputC <- as.matrix(C)
        if (J != nrow(inputC)) {
          private$.logger$error("C matrix must have same number of rows as data")
        }
        
        minCL2 <- min(colSums(inputC^2))
        svdC <- svd(inputC)
        P <- sum(cumsum(rev(svdC$d^2)) > minCL2 * 0.1)
        
        C_matrix <- svdC$u[, 1:P, drop = FALSE]
        C_rotation <- svdC$v[, 1:P, drop = FALSE] %*% diag(1 / svdC$d[1:P], nrow = P)
        
        if (private$.verbose >= 1) {
          private$.logger$info(sprintf("C matrix: %d pathways -> %d components", 
                                       ncol(inputC), P))
        }
      }
      
      # ===================================================================
      # STEP 4: Store K and related parameters
      # ===================================================================
      
      diag_K <- diag(K)
      BcolL2 <- 1
      BrowL2 <- 1
      
      # ===================================================================
      # STEP 5: Organize cells by sample
      # ===================================================================
      
      unique_samples <- unique(Samples)
      num_samples <- length(unique_samples)
      
      cells_by_sample <- split(seq_along(Samples), Samples)[unique_samples]
      Ni <- vapply(cells_by_sample, length, integer(1))
      Ni_vec <- lapply(Ni, function(n) rep(1, n))
      
      cells_all <- unlist(cells_by_sample, use.names = FALSE)
      cellIDs <- cellIDs[cells_all]
      
      if (!is.null(colData)) {
        colData <- colData[cells_all, , drop = FALSE]
      }
      
      # ===================================================================
      # STEP 6: Process prior matrix H (sample-level)
      # ===================================================================
      
      inputH <- NULL
      H_matrix <- matrix(nrow = 0, ncol = 0)
      H_rotation <- NULL
      L <- 0
      
      if (!is.null(H)) {
        inputH <- as.matrix(H)
        inputH <- inputH[, match(unique_samples, colnames(inputH)), drop = FALSE]
        
        if (ncol(inputH) != num_samples) {
          private$.logger$error("H matrix must have columns for all samples")
        }
        
        minHL2 <- min(rowSums(inputH^2))
        svdH <- svd(t(inputH))
        L <- sum(cumsum(rev(svdH$d^2)) > minHL2 * 0.1)
        
        H_matrix <- t(svdH$u[, 1:L, drop = FALSE])
        H_rotation <- t(svdH$v[, 1:L, drop = FALSE] %*% diag(1 / svdH$d[1:L], nrow = L))
        
        if (private$.verbose >= 1) {
          private$.logger$info(sprintf("H matrix: %d covariates -> %d components", 
                                       nrow(inputH), L))
        }
      }
      
      # ===================================================================
      # STEP 7: Initialize hyperparameters (compute s_0 and o_0)
      # ===================================================================
      
      s_0 <- NULL
      o_0 <- NULL
      
      if (is.null(M)) {
        if (is.null(X)) {
          # Y case - use dense functions
          s_0 <- compute_s_0_dense(J_vec, Y_log_dense, J)
          Yp <- compute_Yp_dense(Y_log_dense, J_vec, s_0)
          o_0 <- compute_o_0_dense(Yp, N_vec, N)
        } else {
          # X case
          s_0 <- rep(0, N)
          o_0 <- rep(0, J)
        }
      } else if (is.list(M)) {
        # Paired M case - use sparse functions
        s1_0 <- compute_s_0(J_vec, M[[1]], J)
        Mp1 <- compute_Mp(M[[1]], J_vec, s1_0)
        o1_0 <- compute_o_0(Mp1, N_vec, N)
        
        s2_0 <- compute_s_0(J_vec, M[[2]], J)
        Mp2 <- compute_Mp(M[[2]], J_vec, s2_0)
        o2_0 <- compute_o_0(Mp2, N_vec, N)
        
        s_0 <- log(s1_0 / s2_0)
        o_0 <- log(o1_0 / o2_0)
      } else {
        # Single M case - use sparse functions
        s_0_counts <- compute_s_0(J_vec, M, J)
        Mp <- compute_Mp(M, J_vec, s_0_counts)
        o_0_counts <- compute_o_0(Mp, N_vec, N)
        
        s_0 <- log(s_0_counts)
        o_0 <- log(o_0_counts)
      }
      
      # Split s_0 by sample
      si_0 <- split(s_0, Samples)[unique_samples]
      if (!is.na(fixed_si)) {
        si_0 <- lapply(Ni, function(n) rep(fixed_si, n))
      }
      is_si_fixed <- !is.na(fixed_si)
      
      # Initialize sigma2_0
      if (obs_type != "X") {
        Yp_full <- Y_log_dense - VecVecProduct(J_vec, s_0) - o_0
        sigma2_0 <- sum(Yp_full^2) / (N * J + 1)
      } else {
        sigma2_0 <- 1
      }
      
      # Other hyperparameters
      S_o <- 1 / N / o_shrinkage
      S_si <- 1 / J / si_shrinkage
      S_Z <- 1 / Z_shrinkage
      S_A <- 1 / A_shrinkage
      S_R <- 1 / Rk_shrinkage
      S_oi <- 1 / Ni / oi_shrinkage
      S_oi_mean <- num_samples / N / oi_shrinkage
      S_Qi <- N / Ni / Qi_shrinkage
      S_Qi_mean <- num_samples / Qi_shrinkage
      
      # Generate O matrix for randomized SVD initialization
      O_matrix <- generateRsvdO(m = J, n = N, k = K, p = rsvd_p, sdist = rsvd_sdist)
      
      if (private$.verbose >= 1) {
        private$.logger$info(sprintf("Initial sigma2 = %.6f", sigma2_0))
      }
      
      # ===================================================================
      # STEP 8: Initialize D based on mode
      # ===================================================================
      
      if (mode == "Bl2") {
        D_init <- rep(BrowL2, K)
      } else if (mode == "Bsphere") {
        D_init <- sqrt(BrowL2) * rep(sqrt(K / N), K)
      } else {
        private$.logger$error("mode must be 'Bl2' or 'Bsphere'")
      }
      
      # ===================================================================
      # STEP 9: Initialize model parameters
      # ===================================================================
      
      params <- list(
        Bi = lapply(Ni, function(n) matrix(0, nrow = K, ncol = n)),
        Qi = lapply(1:num_samples, function(i) matrix(0, nrow = J, ncol = K)),
        si = si_0,
        oi = lapply(1:num_samples, function(i) rep(0, J)),
        o = o_0,
        Z = matrix(0, nrow = J, ncol = K),
        U = matrix(0, nrow = J, ncol = K),
        S = rep(1, K),
        D = D_init,
        sigma2 = sigma2_0,
        A = if (P > 0) matrix(0, nrow = P, ncol = K) else matrix(nrow = 0, ncol = 0),
        Rk = if (L > 0) lapply(1:K, function(k) matrix(0, nrow = J, ncol = L)) else list(),
        Ro = if (L > 0) matrix(0, nrow = J, ncol = L) else matrix(nrow = 0, ncol = 0)
      )
      
      # ===================================================================
      # STEP 10: Initialize auxiliary variables
      # ===================================================================
      
      aux <- list(
        ZDBi = lapply(Ni, function(n) matrix(0, nrow = J, ncol = n)),
        QiDBi = lapply(Ni, function(n) matrix(0, nrow = J, ncol = n)),
        Qi_hat = if (L > 0) lapply(1:num_samples, function(i) matrix(0, nrow = J, ncol = K)) else list(),
        oi_hat = if (L > 0) lapply(1:num_samples, function(i) rep(0, J)) else list(),
        J = J,
        N = N,
        K = K,
        P = P,
        L = L,
        numSamples = num_samples,
        diag_K = diag_K,
        J_vec = J_vec,
        Ni = Ni,
        Ni_vec = Ni_vec,
        obs.type = obs_type,
        mode = mode,
        orthoZ = orthoZ,
        adjustD = adjustD,
        is_si_fixed = is_si_fixed,
        C = C_matrix,
        H = H_matrix,
        inputC = inputC,
        C.rotation = C_rotation,
        inputH = inputH,
        H.rotation = H_rotation,
        BcolL2 = BcolL2,
        BrowL2 = BrowL2,
        colData = colData,
        ite = 0
      )
      
      # ===================================================================
      # STEP 11: Prepare target matrices
      # OPTION 2: Don't compute Yi in R - C++ will compute it from M!
      # ===================================================================
      
      if (private$.verbose >= 1) {
        private$.logger$info("Preparing target matrices (C++ will compute Yi)...")
      }
      
      target <- list()
      
      if (obs_type == "M") {
        # Pass ONLY raw sparse M - C++ will compute Yi = log1p(M)
        target$Mi <- lapply(cells_by_sample, function(idx) M[, idx, drop = FALSE])
        target$Yi <- vector("list", num_samples)  # Empty - C++ fills
        target$M1i <- vector("list", num_samples)
        target$M2i <- vector("list", num_samples)
        target$Xi <- vector("list", num_samples)
        
      } else if (obs_type == "M_paired") {
        # Pass ONLY raw sparse M1 and M2 - C++ computes Yi = log((M1+1)/(M2+1))
        target$M1i <- lapply(cells_by_sample, function(idx) M[[1]][, idx, drop = FALSE])
        target$M2i <- lapply(cells_by_sample, function(idx) M[[2]][, idx, drop = FALSE])
        target$Yi <- vector("list", num_samples)  # Empty - C++ fills
        target$Mi <- vector("list", num_samples)
        target$Xi <- vector("list", num_samples)
        
      } else if (obs_type == "X") {
        # Pass X - C++ converts to Yi format
        target$Xi <- lapply(cells_by_sample, function(idx) X[, idx, drop = FALSE])
        target$Yi <- vector("list", num_samples)  # Empty - C++ fills
        target$Mi <- vector("list", num_samples)
        target$M1i <- vector("list", num_samples)
        target$M2i <- vector("list", num_samples)
        
      } else {
        # Y type - must compute Yi in R (no M available)
        target$Yi <- lapply(cells_by_sample, function(idx) as.matrix(Y[, idx, drop = FALSE]))
        target$Mi <- vector("list", num_samples)
        target$M1i <- vector("list", num_samples)
        target$M2i <- vector("list", num_samples)
        target$Xi <- vector("list", num_samples)
      }
      
      # Force garbage collection - Y_log_dense no longer needed
      rm(Y_log_dense)
      invisible(gc())
      
      # ===================================================================
      # STEP 12: Create hyperparams list (CORRECTED - no sigma2_0!)
      # ===================================================================
      
      hyperparams <- list(
        S_Qi = S_Qi,
        S_oi = S_oi,
        S_Z = S_Z,
        S_A = S_A,
        S_R = S_R,
        S_Qi_mean = S_Qi_mean,
        S_oi_mean = S_oi_mean,
        S_si = S_si,
        S_o = S_o,
        o_0 = o_0,
        si_0 = si_0,
        O = O_matrix
        # NO sigma2_0! It goes in params$sigma2
      )
      
      return(list(
        params = params,
        aux = aux,
        target = target,
        hyperparams = hyperparams,
        geneIDs = geneIDs,
        cellIDs = cellIDs,
        sampleNames = unique_samples
      ))
    }
  ),
  
  public = list(
    
    setup = function(
      Samples,
      M = NULL,
      Y = NULL,
      X = NULL,
      colData = NULL,
      C = NULL,
      H = NULL,
      K = 10,
      mode = "Bl2",
      adjustD = TRUE,
      orthoZ = TRUE,
      Z_shrinkage = 1,
      A_shrinkage = 1,
      Qi_shrinkage = 1,
      Rk_shrinkage = 1,
      oi_shrinkage = 1,
      o_shrinkage = 1,
      si_shrinkage = 1,
      fixed_si = NA,
      rsvd_p = 10,
      rsvd_sdist = "normal",
      verbose = 1,
      num_threads = 0
    ) {
      
      private$.verbose <- verbose
      private$.logger <- newLogger("GEDI", level = verbose)
      
      private$.logger$info("Setting up GEDI model...")
      
      data <- private$prepareData(
        Samples = Samples, Y = Y, X = X, M = M,
        colData = colData, C = C, H = H, K = K,
        mode = mode, adjustD = adjustD, orthoZ = orthoZ,
        Z_shrinkage = Z_shrinkage, A_shrinkage = A_shrinkage,
        Qi_shrinkage = Qi_shrinkage, Rk_shrinkage = Rk_shrinkage,
        oi_shrinkage = oi_shrinkage, o_shrinkage = o_shrinkage,
        si_shrinkage = si_shrinkage, fixed_si = fixed_si,
        rsvd_p = rsvd_p, rsvd_sdist = rsvd_sdist
      )
      
      private$.geneIDs <- data$geneIDs
      private$.cellIDs <- data$cellIDs
      private$.sampleNames <- data$sampleNames
      
      # Create C++ object (C++ will compute Yi from M)
      private$.cppPtr <- GEDI_new(
        params = data$params,
        aux = data$aux,
        target = data$target,
        hyperparams = data$hyperparams,
        verbose = private$.verbose,
        num_threads = num_threads
      )
      
      # CRITICAL: Force cleanup of R-side data to save memory
      data$target <- NULL
      data$params$Bi <- NULL
      data$params$Qi <- NULL
      data$aux$ZDBi <- NULL
      data$aux$QiDBi <- NULL
      rm(data)
      gc(full = TRUE)
      
      private$.logger$info("Setup complete.")
      
      invisible(self)
    },
    
    initialize_lvs = function(multimodal = FALSE) {
      if (is.null(private$.cppPtr)) {
        stop("Model not setup. Call CreateGEDIObject() first.", call. = FALSE)
      }
      
      private$.logger$info("Initializing latent variables...")
      
      private$.lastResult <- GEDI_initialize(private$.cppPtr, multimodal = multimodal)
      private$.isInitialized <- TRUE
      
      private$.logger$info("Initialization complete.")
      
      invisible(self)
    },
    
    optimize = function(iterations = 50, track_interval = 5) {
      if (!private$.isInitialized) {
        stop("Model not initialized. Call $initialize_lvs() first.", call. = FALSE)
      }
      
      private$.logger$info(sprintf("Running optimization (%d iterations)...", iterations))
      
      private$.lastResult <- GEDI_optimize(
        model_ptr = private$.cppPtr,
        iterations = iterations,
        track_interval = track_interval
      )
      
      private$.logger$info("Optimization complete.")
      
      return(private$.lastResult)
    },
    
    train = function(iterations = 50, track_interval = 5, multimodal = FALSE) {
      if (is.null(private$.cppPtr)) {
        stop("Model not setup. Call CreateGEDIObject() first.", call. = FALSE)
      }
      
      private$.logger$info("Training GEDI model...")
      
      private$.lastResult <- GEDI_train(
        model_ptr = private$.cppPtr,
        iterations = iterations,
        track_interval = track_interval,
        multimodal = multimodal
      )
      
      private$.isInitialized <- TRUE
      
      private$.logger$info("Training complete.")
      
      return(private$.lastResult)
    },
    
    print = function() {
      cat("<GEDI Model>\n")
      
      if (is.null(private$.cppPtr)) {
        cat("Status: Not setup\n")
        return(invisible(self))
      }
      
      if (!is.null(private$.lastResult)) {
        aux <- private$.lastResult$aux
        cat(sprintf("Dimensions: %d genes × %d cells\n", aux$J, aux$N))
        cat(sprintf("Samples: %d (%s)\n", 
                    aux$numSamples, 
                    paste(private$.sampleNames, collapse = ", ")))
        cat(sprintf("Latent factors: K = %d\n", aux$K))
        cat(sprintf("Mode: %s (orthoZ = %s)\n", aux$mode, aux$orthoZ))
        
        if (aux$P > 0) cat(sprintf("Prior C: %d pathways\n", aux$P))
        if (aux$L > 0) cat(sprintf("Prior H: %d covariates\n", aux$L))
      }
      
      if (private$.isInitialized) {
        sigma2 <- private$.lastResult$params$sigma2
        cat(sprintf("\nStatus: Fitted (sigma2 = %.6f)\n", sigma2))
      } else {
        cat("\nStatus: Setup complete (not yet initialized)\n")
      }
      
      invisible(self)
    }
  ),
  
  active = list(
    
    params = function(value) {
      if (!missing(value)) stop("params is read-only", call. = FALSE)
      if (is.null(private$.lastResult)) {
        stop("No results yet. Run $initialize_lvs() or $train() first.", call. = FALSE)
      }
      private$.lastResult$params
    },
    
    aux = function(value) {
      if (!missing(value)) stop("aux is read-only", call. = FALSE)
      if (is.null(private$.lastResult)) {
        stop("No results yet. Run $initialize_lvs() or $train() first.", call. = FALSE)
      }
      private$.lastResult$aux
    },
    
    hyperparams = function(value) {
      if (!missing(value)) stop("hyperparams is read-only", call. = FALSE)
      if (is.null(private$.lastResult)) {
        stop("No results yet. Run $initialize_lvs() or $train() first.", call. = FALSE)
      }
      private$.lastResult$hyperparams
    },
    
    tracking = function(value) {
      if (!missing(value)) stop("tracking is read-only", call. = FALSE)
      if (is.null(private$.lastResult)) {
        stop("No results yet. Run $optimize() or $train() first.", call. = FALSE)
      }
      private$.lastResult$tracking
    },
    
    Z = function(value) {
      if (!missing(value)) stop("Z is read-only", call. = FALSE)
      if (is.null(private$.lastResult)) {
        stop("No results yet. Run $initialize_lvs() or $train() first.", call. = FALSE)
      }
      private$.lastResult$params$Z
    },
    
    geneIDs = function(value) {
      if (!missing(value)) stop("geneIDs is read-only", call. = FALSE)
      private$.geneIDs
    },
    
    cellIDs = function(value) {
      if (!missing(value)) stop("cellIDs is read-only", call. = FALSE)
      private$.cellIDs
    },
    
    sampleNames = function(value) {
      if (!missing(value)) stop("sampleNames is read-only", call. = FALSE)
      private$.sampleNames
    }
  )
)

#' Create GEDI Object
#'
#' @description
#' Creates and configures a GEDI (Gene Expression Data Integration) model object.
#' This implementation uses Option 2 memory optimization: C++ computes Yi from M,
#' eliminating duplicate Yi storage in R.
#' 
#' @param Samples Factor or character vector indicating sample of origin for each cell
#' @param M Raw count matrix (sparse or dense), or list of two matrices for paired data.
#'   C++ will compute Yi = log(M+1) internally - no Yi copy stored in R!
#' @param Y Log-transformed expression matrix (optional if M provided)
#' @param X Binary indicator matrix (optional if M or Y provided)
#' @param colData Optional data.frame with cell metadata
#' @param C Gene-level prior matrix (genes × pathways)
#' @param H Sample-level covariate matrix (covariates × samples)
#' @param K Number of latent factors (default: 10)
#' @param mode Normalization mode: "Bl2" or "Bsphere" (default: "Bl2")
#' @param adjustD Whether to adjust D based on B row norms (default: TRUE)
#' @param orthoZ Whether Z columns should be orthogonal (default: TRUE)
#' @param Z_shrinkage,A_shrinkage,Qi_shrinkage,Rk_shrinkage,oi_shrinkage,o_shrinkage,si_shrinkage 
#'   Regularization strengths (default: 1)
#' @param fixed_si Fix cell library sizes at this value, or NA to optimize (default: NA)
#' @param rsvd_p Oversampling parameter for randomized SVD (default: 10)
#' @param rsvd_sdist Random distribution for rSVD: "normal", "unif", or "rademacher" 
#'   (default: "normal")
#' @param verbose Verbosity level: 0 (silent), 1 (info), 2 (debug) (default: 1)
#' @param num_threads Number of OpenMP threads (default: 0 = auto)
#'
#' @return GEDI R6 object with memory-efficient architecture
#' 
#' @examples
#' \dontrun{
#' # Basic usage - memory efficient!
#' model <- CreateGEDIObject(
#'   Samples = sample_labels,
#'   M = count_matrix,  # Only M stored in R; Yi computed in C++
#'   K = 15,
#'   num_threads = 32
#' )
#' 
#' # Train the model
#' result <- model$train(iterations = 50)
#' 
#' # Access results
#' Z <- model$Z
#' params <- model$params
#' }
#' 
#' @export
CreateGEDIObject <- function(
  Samples,
  M = NULL,
  Y = NULL,
  X = NULL,
  colData = NULL,
  C = NULL,
  H = NULL,
  K = 10,
  mode = "Bl2",
  adjustD = TRUE,
  orthoZ = TRUE,
  Z_shrinkage = 1,
  A_shrinkage = 1,
  Qi_shrinkage = 1,
  Rk_shrinkage = 1,
  oi_shrinkage = 1,
  o_shrinkage = 1,
  si_shrinkage = 1,
  fixed_si = NA,
  rsvd_p = 10,
  rsvd_sdist = "normal",
  verbose = 1,
  num_threads = 0
) {
  
  model <- GEDI$new()
  
  model$setup(
    Samples = Samples,
    M = M, Y = Y, X = X,
    colData = colData,
    C = C, H = H,
    K = K,
    mode = mode,
    adjustD = adjustD,
    orthoZ = orthoZ,
    Z_shrinkage = Z_shrinkage,
    A_shrinkage = A_shrinkage,
    Qi_shrinkage = Qi_shrinkage,
    Rk_shrinkage = Rk_shrinkage,
    oi_shrinkage = oi_shrinkage,
    o_shrinkage = o_shrinkage,
    si_shrinkage = si_shrinkage,
    fixed_si = fixed_si,
    rsvd_p = rsvd_p,
    rsvd_sdist = rsvd_sdist,
    verbose = verbose,
    num_threads = num_threads
  )
  
  return(model)
}
