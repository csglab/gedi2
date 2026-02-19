#' Write GEDI model to H5AD file
#'
#' Exports a trained GEDI model to H5AD (AnnData) format for interoperability
#' with Python tools like scanpy. The file contains expression data, embeddings,
#' metadata, and GEDI-specific parameters.
#'
#' @param model GEDI R6 object (must be trained)
#' @param file_path Character. Path where the H5AD file should be written.
#' @param X_slot Character. Which expression data to save in the main X slot:
#'   \itemize{
#'     \item "imputed": Imputed expression (default, requires M matrix)
#'     \item "projection": ZDB projection (always available)
#'     \item "original": Original M matrix (requires M parameter)
#'   }
#' @param M Optional. Original count matrix to save when X_slot="original" or
#'   include_raw=TRUE. Must match the dimensions used during model training.
#' @param include_embeddings Logical. Include PCA/UMAP in obsm if cached. Default TRUE.
#' @param include_raw Logical. If TRUE, saves original M in raw.X (requires M parameter).
#'   Default FALSE.
#' @param compression Integer. Gzip compression level (0-9). Default 6.
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return Invisibly returns the file path
#'
#' @details
#' The H5AD file structure contains:
#' \itemize{
#'   \item X: Main expression matrix (based on X_slot parameter)
#'   \item obs: Cell metadata (sample IDs, colData)
#'   \item var: Gene metadata (gene IDs)
#'   \item obsm: Cell embeddings (X_gedi, X_pca, X_umap)
#'   \item varm: Gene loadings (gedi_Z, gedi_Q_mean)
#'   \item uns: Model parameters and metadata
#'   \item raw.X: Original counts (if include_raw=TRUE)
#' }
#'
#' The function handles the technical details of HDF5/AnnData compatibility:
#' \itemize{
#'   \item Writes sparse matrices in CSR format
#'   \item Transposes dense matrices for Python's row-major layout
#'   \item Creates scalar string attributes (not arrays) for AnnData compatibility
#'   \item Validates M matrix identity using fingerprints
#' }
#'
#' @examples
#' \donttest{
#' # Train a GEDI model
#' model <- CreateGEDIObject(Samples = samples, M = counts, K = 15)
#' model$train(iterations = 50)
#'
#' # Save to H5AD with imputed expression
#' write_h5ad(model, "results.h5ad")
#'
#' # Save with original counts in raw slot
#' write_h5ad(model, "results.h5ad", M = counts, include_raw = TRUE)
#'
#' # Save original M matrix as main X slot
#' write_h5ad(model, "results.h5ad", X_slot = "original", M = counts)
#'
#' # Save projection instead of imputed values
#' write_h5ad(model, "results.h5ad", X_slot = "projection")
#'
#' # Read in Python:
#' # import scanpy as sc
#' # adata = sc.read_h5ad("results.h5ad")
#' # adata.obsm['X_gedi']  # GEDI embeddings
#' }
#'
#' @importFrom Matrix t
#' @export
write_h5ad <- function(model,
                       file_path,
                       X_slot = c("imputed", "projection", "original"),
                       M = NULL,
                       include_embeddings = TRUE,
                       include_raw = FALSE,
                       compression = 6,
                       verbose = TRUE) {

  if (verbose) message("[write_h5ad] Starting H5AD file writing process...")

  # Validate dependencies
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("[write_h5ad] ERROR: Package 'hdf5r' is required but not installed.\n",
         "  Install with: install.packages('hdf5r')")
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("[write_h5ad] ERROR: Package 'Matrix' is required but not installed.\n",
         "  Install with: install.packages('Matrix')")
  }

  # Validate model
  if (!inherits(model, "GEDI")) {
    stop("[write_h5ad] ERROR: model must be a GEDI R6 object")
  }

  if (is.null(model$params)) {
    stop("[write_h5ad] ERROR: Model not trained. Run model$train() first.")
  }

  X_slot <- match.arg(X_slot)

  # Get model components early to check obs type
  aux <- model$aux

  # Validate M if provided
  if (!is.null(M)) {
    # Check if model has M fingerprint (only for M/M_paired obs types)
    if (aux$obs.type %in% c("M", "M_paired")) {
      # Access private environment to get M fingerprint
      model_env <- environment(model$train)
      private_env <- model_env$private

      if (!is.null(private_env$.M_fingerprint)) {
        if (verbose) message("[write_h5ad] Validating M matrix identity...")
        validate_M_identity(M, private_env$.M_fingerprint)
        if (verbose) message("[write_h5ad]   M matrix validated successfully")
      } else {
        warning("[write_h5ad] No M fingerprint stored in model. Skipping validation.")
      }
    }
  }

  # Check if file exists
  if (file.exists(file_path)) {
    warning("[write_h5ad] File already exists and will be overwritten: ", file_path)
    file.remove(file_path)
  }

  # Create H5 file
  if (verbose) message("[write_h5ad] Creating H5AD file: ", file_path)
  h5file <- tryCatch({
    hdf5r::H5File$new(file_path, mode = "w")
  }, error = function(e) {
    stop("[write_h5ad] ERROR: Failed to create H5AD file.\n",
         "  File: ", file_path, "\n",
         "  Error: ", conditionMessage(e))
  })

  on.exit({
    if (!is.null(h5file)) {
      h5file$close_all()
      if (verbose) message("[write_h5ad] H5 file closed successfully")
    }
  }, add = TRUE)

  # Get model components
  metadata <- model$metadata
  params <- model$params
  aux <- model$aux

  J <- aux$J
  N <- aux$N
  K <- aux$K

  if (verbose) {
    message("[write_h5ad] Model dimensions: ", J, " genes x ", N, " cells")
    message("[write_h5ad] Latent factors: K = ", K)
  }

  # =========================================================================
  # STEP 1: Prepare and write main expression matrix (X)
  # =========================================================================

  if (verbose) message("[write_h5ad] Preparing main expression matrix (X)...")

  X_matrix <- NULL
  X_description <- ""

  if (X_slot == "imputed") {
    # Check if imputation is available
    if (aux$obs.type %in% c("M", "M_paired")) {
      if (verbose) message("[write_h5ad]   Computing imputed expression...")
      X_matrix <- model$imputed$Y()
      X_description <- "GEDI imputed expression"
    } else {
      warning("[write_h5ad] Imputed expression not available (no M matrix). Using projection instead.")
      X_matrix <- model$projections$ZDB
      X_description <- "GEDI projection (ZDB)"
    }
  } else if (X_slot == "projection") {
    if (verbose) message("[write_h5ad]   Using ZDB projection...")
    X_matrix <- model$projections$ZDB
    X_description <- "GEDI projection (ZDB)"
  } else if (X_slot == "original") {
    if (is.null(M)) {
      stop("[write_h5ad] ERROR: X_slot='original' requires M parameter.\n",
           "  Usage: write_h5ad(model, file, X_slot='original', M=count_matrix)")
    }

    if (verbose) message("[write_h5ad]   Using original M matrix...")

    # Handle paired M case
    if (is.list(M)) {
      if (verbose) message("[write_h5ad]   Detected paired M matrices - using M[[1]]...")
      X_matrix <- M[[1]]
      X_description <- "Original counts (M[[1]] from paired data)"
    } else {
      X_matrix <- M
      X_description <- "Original counts"
    }

    # Validate dimensions
    if (nrow(X_matrix) != J || ncol(X_matrix) != N) {
      stop("[write_h5ad] ERROR: M dimensions (", nrow(X_matrix), " x ", ncol(X_matrix),
           ") don't match model dimensions (", J, " x ", N, ")")
    }
  }

  # Transpose to cells x genes (AnnData format)
  X_matrix <- Matrix::t(X_matrix)

  if (verbose) message("[write_h5ad]   Writing X matrix (", nrow(X_matrix), " x ", ncol(X_matrix), ")...")
  .write_h5ad_sparse_matrix(h5file, "X", X_matrix, compression, verbose)

  # =========================================================================
  # STEP 2: Write cell metadata (obs)
  # =========================================================================

  if (verbose) message("[write_h5ad] Writing cell metadata (obs)...")

  obs_df <- data.frame(
    row.names = metadata$cellIDs,
    sample = metadata$Samples,
    stringsAsFactors = FALSE
  )

  # Add colData if available
  if (!is.null(metadata$colData)) {
    obs_df <- cbind(obs_df, metadata$colData)
  }

  .write_h5ad_dataframe(h5file, "obs", obs_df, verbose)

  # =========================================================================
  # STEP 3: Write gene metadata (var)
  # =========================================================================

  if (verbose) message("[write_h5ad] Writing gene metadata (var)...")

  var_df <- data.frame(
    row.names = metadata$geneIDs,
    gene_ids = metadata$geneIDs,
    stringsAsFactors = FALSE
  )

  .write_h5ad_dataframe(h5file, "var", var_df, verbose)

  # =========================================================================
  # STEP 4: Write multi-dimensional cell annotations (obsm)
  # =========================================================================

  if (verbose) message("[write_h5ad] Writing cell embeddings (obsm)...")

  obsm_group <- h5file$create_group("obsm")

  # Always include GEDI embeddings (DB)
  # DB is K x N, need cells x K for AnnData
  # HDF5 stores in column-major (R), Python reads row-major (C)
  # Write transposed so Python reads correct dimensions
  DB <- Matrix::t(model$projections$DB)  # K x N -> N x K
  .write_h5ad_dense_matrix(obsm_group, "X_gedi", t(as.matrix(DB)), verbose)

  # Include PCA and UMAP if cached and requested
  if (include_embeddings) {
    cache_status <- model$cache_status()

    if (cache_status["pca"]) {
      if (verbose) message("[write_h5ad]   Including PCA coordinates...")
      pca <- model$get_pca()
      # PCA needs transpose for Python compatibility
      .write_h5ad_dense_matrix(obsm_group, "X_pca", t(pca), verbose)
    }

    # Check if UMAP is available (would be in cache if computed)
    # Note: We don't have a direct cache flag for UMAP, so we skip for now
    # Users can compute UMAP separately and add it
  }

  # =========================================================================
  # STEP 5: Write multi-dimensional gene annotations (varm)
  # =========================================================================

  if (verbose) message("[write_h5ad] Writing gene loadings (varm)...")

  varm_group <- h5file$create_group("varm")

  # Write Z matrix (gene factor loadings)
  # Z is J x K, need to transpose for Python compatibility
  .write_h5ad_dense_matrix(varm_group, "gedi_Z", t(params$Z), verbose)

  # Write mean Q across samples
  Q_mean <- Reduce("+", params$Qi) / length(params$Qi)
  # Q is J x K, need to transpose for Python compatibility
  .write_h5ad_dense_matrix(varm_group, "gedi_Q_mean", t(Q_mean), verbose)

  # =========================================================================
  # STEP 6: Write unstructured metadata (uns)
  # =========================================================================

  if (verbose) message("[write_h5ad] Writing model parameters (uns)...")

  uns_group <- h5file$create_group("uns")
  gedi_group <- uns_group$create_group("gedi")

  # Model parameters
  gedi_group[["K"]] <- K
  gedi_group[["sigma2"]] <- params$sigma2
  gedi_group[["mode"]] <- aux$mode
  gedi_group[["orthoZ"]] <- aux$orthoZ
  gedi_group[["adjustD"]] <- aux$adjustD
  gedi_group[["obs_type"]] <- aux$obs.type

  # Write D and S vectors
  .write_h5ad_vector(gedi_group, "D", params$D)
  .write_h5ad_vector(gedi_group, "S", params$S)

  # Write o vector
  .write_h5ad_vector(gedi_group, "o", params$o)

  # Metadata
  gedi_group[["n_genes"]] <- J
  gedi_group[["n_cells"]] <- N
  gedi_group[["n_samples"]] <- aux$numSamples

  # Sample names
  if (length(metadata$sampleIDs) > 0) {
    .write_h5ad_string_array(gedi_group, "sample_names", metadata$sampleIDs)
  }

  # Store description of X slot
  uns_group[["X_description"]] <- X_description

  # =========================================================================
  # STEP 7: Write priors if available
  # =========================================================================

  priors <- model$priors
  if (!is.null(priors$inputC) || !is.null(priors$inputH)) {
    if (verbose) message("[write_h5ad] Writing prior matrices...")

    priors_group <- gedi_group$create_group("priors")

    if (!is.null(priors$inputC)) {
      # C is J x P, transpose for Python compatibility
      .write_h5ad_dense_matrix(priors_group, "C", t(priors$inputC), verbose)
      priors_group[["P"]] <- aux$P
    }

    if (!is.null(priors$inputH)) {
      # H is J x L, transpose for Python compatibility
      .write_h5ad_dense_matrix(priors_group, "H", t(priors$inputH), verbose)
      priors_group[["L"]] <- aux$L
    }
  }

  # =========================================================================
  # STEP 8: Write raw counts if requested
  # =========================================================================

  if (include_raw) {
    if (is.null(M)) {
      warning("[write_h5ad] include_raw=TRUE requires M parameter. Skipping raw.X.\n",
              "  Usage: write_h5ad(model, file, include_raw=TRUE, M=count_matrix)")
    } else {
      if (verbose) message("[write_h5ad] Writing raw counts (raw.X)...")

      # Create raw group
      raw_group <- h5file$create_group("raw")

      # Prepare raw matrix
      raw_matrix <- if (is.list(M)) M[[1]] else M

      # Validate dimensions
      if (nrow(raw_matrix) != J || ncol(raw_matrix) != N) {
        warning("[write_h5ad] M dimensions don't match model dimensions. Skipping raw.X.")
      } else {
        # Transpose to cells x genes (AnnData format)
        raw_matrix <- Matrix::t(raw_matrix)

        # Write raw.X
        .write_h5ad_sparse_matrix(raw_group, "X", raw_matrix, compression, verbose)

        # Write raw.var (same as main var)
        .write_h5ad_dataframe(raw_group, "var", var_df, verbose = FALSE)
      }
    }
  }

  # =========================================================================
  # STEP 9: Set top-level attributes
  # =========================================================================

  # Note: encoding-type and encoding-version attributes removed
  # They cause issues with AnnData due to being stored as arrays
  # AnnData can still read the file based on structure alone

  if (verbose) {
    file_size_mb <- file.info(file_path)$size / 1024^2
    message("[write_h5ad] SUCCESS: H5AD file written")
    message("[write_h5ad] File size: ", round(file_size_mb, 2), " MB")
    message("[write_h5ad] Path: ", file_path)
  }

  invisible(file_path)
}


#' Write sparse matrix to H5AD file (CSR format)
#' @keywords internal
#' @noRd
.write_h5ad_sparse_matrix <- function(h5file, name, matrix, compression = 6, verbose = TRUE) {

  # Convert to dgCMatrix if not already
  if (!inherits(matrix, "dgCMatrix")) {
    matrix <- Matrix::Matrix(matrix, sparse = TRUE)
  }

  # Convert to CSR format (AnnData standard)
  # dgCMatrix is CSC (column-compressed), we need CSR (row-compressed)
  # Transpose, convert, then adjust
  matrix_csr <- as(Matrix::t(as(Matrix::t(matrix), "CsparseMatrix")), "RsparseMatrix")

  # Create group for sparse matrix
  mat_group <- h5file$create_group(name)

  # Get CSR components
  # In R's dgRMatrix: x = data, j = column indices, p = row pointers
  data <- matrix_csr@x
  indices <- matrix_csr@j  # 0-based in R sparse matrices
  indptr <- matrix_csr@p   # 0-based

  # Write components
  mat_group$create_dataset(
    "data",
    robj = data,
    dtype = hdf5r::h5types$H5T_NATIVE_FLOAT,
    space = hdf5r::H5S$new("simple", dims = length(data), maxdims = length(data)),
    gzip_level = compression
  )

  mat_group$create_dataset(
    "indices",
    robj = as.integer(indices),
    dtype = hdf5r::h5types$H5T_NATIVE_INT32,
    space = hdf5r::H5S$new("simple", dims = length(indices), maxdims = length(indices)),
    gzip_level = compression
  )

  mat_group$create_dataset(
    "indptr",
    robj = as.integer(indptr),
    dtype = hdf5r::h5types$H5T_NATIVE_INT32,
    space = hdf5r::H5S$new("simple", dims = length(indptr), maxdims = length(indptr)),
    gzip_level = compression
  )

  # Write attributes for AnnData compatibility
  hdf5r::h5attr(mat_group, "shape") <- as.integer(c(nrow(matrix), ncol(matrix)))

  # Write string attributes as scalars for Python compatibility
  .write_h5ad_scalar_string_attr(mat_group, "encoding-type", "csr_matrix")
  .write_h5ad_scalar_string_attr(mat_group, "encoding-version", "0.1.0")

  if (verbose) {
    sparsity <- 100 * (1 - length(data) / (nrow(matrix) * ncol(matrix)))
    message("[write_h5ad]     Sparse matrix '", name, "': ",
            nrow(matrix), " x ", ncol(matrix),
            " (", round(sparsity, 1), "% sparse)")
  }
}


#' Write dense matrix to H5AD file
#' @keywords internal
#' @noRd
.write_h5ad_dense_matrix <- function(h5group, name, matrix, verbose = TRUE) {

  # Convert to regular matrix if sparse
  if (inherits(matrix, "sparseMatrix")) {
    matrix <- as.matrix(matrix)
  }

  # Ensure it's a matrix
  if (!is.matrix(matrix)) {
    matrix <- as.matrix(matrix)
  }

  # Create dataset
  h5group$create_dataset(
    name,
    robj = matrix,
    dtype = hdf5r::h5types$H5T_NATIVE_FLOAT
  )

  # Note: encoding attributes removed for AnnData compatibility

  if (verbose) {
    message("[write_h5ad]     Dense matrix '", name, "': ",
            nrow(matrix), " x ", ncol(matrix))
  }
}


#' Write data frame to H5AD file
#' @keywords internal
#' @noRd
.write_h5ad_dataframe <- function(h5file, name, df, verbose = TRUE) {

  df_group <- h5file$create_group(name)

  # Write _index (row names)
  if (!is.null(rownames(df))) {
    .write_h5ad_string_array(df_group, "_index", rownames(df))
  }

  # Write columns
  for (col_name in colnames(df)) {
    col_data <- df[[col_name]]

    if (is.factor(col_data)) {
      # Write as categorical
      col_group <- df_group$create_group(col_name)

      categories <- levels(col_data)
      codes <- as.integer(col_data) - 1L  # 0-based indexing

      .write_h5ad_string_array(col_group, "categories", categories)

      col_group$create_dataset(
        "codes",
        robj = codes,
        dtype = hdf5r::h5types$H5T_NATIVE_INT32
      )

      # Note: encoding attributes removed for AnnData compatibility
      hdf5r::h5attr(col_group, "ordered") <- FALSE

    } else if (is.character(col_data)) {
      .write_h5ad_string_array(df_group, col_name, col_data)

    } else if (is.numeric(col_data)) {
      df_group$create_dataset(
        col_name,
        robj = col_data,
        dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE
      )
      # Note: encoding attributes removed for AnnData compatibility

    } else {
      warning("[write_h5ad] Skipping column '", col_name,
              "' with unsupported type: ", class(col_data)[1])
    }
  }

  # Write dataframe encoding attributes for AnnData compatibility
  .write_h5ad_scalar_string_attr(df_group, "_index", "_index")
  .write_h5ad_scalar_string_attr(df_group, "encoding-type", "dataframe")
  .write_h5ad_scalar_string_attr(df_group, "encoding-version", "0.2.0")

  # Write column-order attribute
  if (length(colnames(df)) > 0) {
    # Column-order should be an attribute, not a dataset
    col_order_dtype <- hdf5r::H5T_STRING$new(size = Inf)
    col_order_dtype$set_cset(hdf5r::h5const$H5T_CSET_UTF8)
    space <- hdf5r::H5S$new("simple", dims = length(colnames(df)), maxdims = length(colnames(df)))
    attr <- df_group$create_attr("column-order", dtype = col_order_dtype, space = space)
    attr$write(colnames(df))
    attr$close()
  } else {
    # Empty column-order for dataframes with no columns (just index)
    hdf5r::h5attr(df_group, "column-order") <- numeric(0)
  }

  if (verbose) {
    message("[write_h5ad]   DataFrame '", name, "': ",
            nrow(df), " rows, ", ncol(df), " columns")
  }
}


#' Write string array to H5AD file
#' @keywords internal
#' @noRd
.write_h5ad_string_array <- function(h5group, name, strings) {

  # Convert to character
  strings <- as.character(strings)

  # Create variable-length string datatype
  str_dtype <- hdf5r::H5T_STRING$new(size = Inf)
  str_dtype$set_cset(hdf5r::h5const$H5T_CSET_UTF8)

  # Create dataset
  h5group$create_dataset(
    name,
    robj = strings,
    dtype = str_dtype
  )

  # Note: encoding attributes removed for AnnData compatibility
}


#' Write numeric vector to H5AD file
#' @keywords internal
#' @noRd
.write_h5ad_vector <- function(h5group, name, vec) {

  h5group$create_dataset(
    name,
    robj = as.numeric(vec),
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE
  )

  # Note: encoding attributes removed for AnnData compatibility
}


#' Write scalar string attribute to H5 object
#'
#' Writes a string attribute as a true scalar (not 1-element array) for Python compatibility.
#' By default, hdf5r writes string attributes as 1-element arrays which causes issues with AnnData.
#'
#' @keywords internal
#' @noRd
.write_h5ad_scalar_string_attr <- function(h5obj, attr_name, attr_value) {

  # Create variable-length string datatype
  str_type <- hdf5r::H5T_STRING$new(size = Inf)
  str_type$set_cset(hdf5r::h5const$H5T_CSET_UTF8)

  # Create scalar dataspace (this is key - not "simple" with dims=1)
  space <- hdf5r::H5S$new(type = "scalar")

  # Create and write attribute
  attr <- h5obj$create_attr(attr_name, dtype = str_type, space = space)
  attr$write(attr_value)
  attr$close()
}
