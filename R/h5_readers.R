#' Read H5AD file and convert to sparse matrix
#'
#' Reads an H5AD file (AnnData format) and extracts the expression matrix
#' as a sparse dgCMatrix suitable for use with gedi R6 object.
#'
#' @param file_path Character. Path to the H5AD file.
#' @param layer Character. The layer to extract from the H5AD file.
#'   Default is NULL, which reads from X (the main expression matrix).
#'   Common alternatives include "counts", "data", "scaled", etc.
#' @param use_raw Logical. If TRUE, reads from the raw.X slot instead of X.
#'   Default is FALSE.
#' @param transpose Logical. If TRUE, transposes the matrix so genes are rows
#'   and cells are columns (gedi format). Default is TRUE.
#' @param return_metadata Logical. If TRUE, returns a list with the expression
#'   matrix, cell metadata (obs), and gene metadata (var). If FALSE, returns
#'   only the expression matrix. Default is FALSE.
#' @param verbose Logical. If TRUE, prints detailed progress messages. Default is TRUE.
#'
#' @return If return_metadata = FALSE, returns a sparse matrix (dgCMatrix) with
#'   genes as rows and cells as columns. If return_metadata = TRUE, returns a
#'   list with components:
#'   \itemize{
#'     \item X: sparse expression matrix (genes x cells)
#'     \item obs: data.frame of cell metadata
#'     \item var: data.frame of gene metadata
#'   }
#'
#' @details
#' This function reads H5AD files, which are the standard format for AnnData
#' objects in Python. Compatible with gedi R6 class for seamless integration.
#'
#' @examples
#' \dontrun{
#' # Read H5AD file with default settings
#' expr_matrix <- read_h5ad("data.h5ad")
#'
#' # Use with gedi R6 object
#' library(gedi)
#' expr_matrix <- read_h5ad("data.h5ad")
#' samples <- factor(rep(c("sample1", "sample2"), each = ncol(expr_matrix)/2))
#' gedi_obj <- gedi$new(n_sample = 2)
#' gedi_obj$setup(Y = expr_matrix, sample_id = samples, K = 10)
#' }
#'
#' @importFrom hdf5r H5File
#' @importFrom Matrix sparseMatrix t Matrix
#' @export
read_h5ad <- function(file_path,
                      layer = NULL,
                      use_raw = FALSE,
                      transpose = TRUE,
                      return_metadata = FALSE,
                      verbose = TRUE) {

  if (verbose) message("[read_h5ad] Starting H5AD file reading process...")

  # Check dependencies
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("[read_h5ad] ERROR: Package 'hdf5r' is required but not installed.\n",
         "  Install with: install.packages('hdf5r')\n",
         "  Or use: gedi::install_h5_dependencies()")
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("[read_h5ad] ERROR: Package 'Matrix' is required but not installed.\n",
         "  Install with: install.packages('Matrix')")
  }

  # Validate file path
  if (!file.exists(file_path)) {
    stop("[read_h5ad] ERROR: File not found: ", file_path)
  }

  file_size_mb <- file.info(file_path)$size / 1024^2
  if (verbose) message("[read_h5ad] File size: ", round(file_size_mb, 2), " MB")

  # Open H5AD file with error handling
  if (verbose) message("[read_h5ad] Opening H5AD file...")
  h5file <- tryCatch({
    hdf5r::H5File$new(file_path, mode = "r")
  }, error = function(e) {
    stop("[read_h5ad] ERROR: Failed to open H5AD file.\n",
         "  File: ", file_path, "\n",
         "  Error: ", conditionMessage(e))
  })

  # Ensure file is closed on exit
  on.exit({
    if (!is.null(h5file)) {
      h5file$close()
      if (verbose) message("[read_h5ad] H5 file closed successfully")
    }
  }, add = TRUE)

  # Determine matrix path
  if (use_raw) {
    matrix_path <- "raw/X"
    var_path <- "raw/var"
    if (verbose) message("[read_h5ad] Reading from: raw.X (raw counts)")
  } else if (!is.null(layer)) {
    matrix_path <- paste0("layers/", layer)
    var_path <- "var"
    if (verbose) message("[read_h5ad] Reading from layer: ", layer)
  } else {
    matrix_path <- "X"
    var_path <- "var"
    if (verbose) message("[read_h5ad] Reading from: X (main expression matrix)")
  }

  # Check if matrix path exists
  if (!h5file$exists(matrix_path)) {
    available_paths <- tryCatch({
      h5file$ls(recursive = TRUE)$name
    }, error = function(e) {
      character(0)
    })
    stop("[read_h5ad] ERROR: Matrix path '", matrix_path, "' not found in H5AD file.\n",
         "  Available paths (first 20): ", paste(head(available_paths, 20), collapse = ", "), "\n",
         "  Use list_h5_structure('", file_path, "') to explore the file structure.")
  }

  # Read expression matrix
  if (verbose) message("[read_h5ad] Reading expression matrix from: ", matrix_path)
  expr_matrix <- tryCatch({
    .read_h5ad_matrix(h5file, matrix_path, verbose)
  }, error = function(e) {
    stop("[read_h5ad] ERROR: Failed to read expression matrix.\n",
         "  Matrix path: ", matrix_path, "\n",
         "  Error: ", conditionMessage(e))
  })

  if (verbose) {
    message("[read_h5ad] Matrix dimensions (before transpose): ",
            nrow(expr_matrix), " x ", ncol(expr_matrix))
    message("[read_h5ad] Matrix class: ", class(expr_matrix)[1])
    message("[read_h5ad] Sparsity: ",
            round(100 * (1 - Matrix::nnzero(expr_matrix) / (nrow(expr_matrix) * ncol(expr_matrix))), 2), "%")
  }

  # Transpose if needed (H5AD stores as cells x genes, gedi needs genes x cells)
  if (transpose) {
    if (verbose) message("[read_h5ad] Transposing matrix to genes x cells format...")
    expr_matrix <- Matrix::t(expr_matrix)
    if (verbose) message("[read_h5ad] New dimensions: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")
  }

  # Read metadata if requested
  if (return_metadata) {
    if (verbose) message("[read_h5ad] Reading metadata...")

    # Read cell metadata (obs)
    obs <- NULL
    if (h5file$exists("obs")) {
      if (verbose) message("[read_h5ad]   Reading cell metadata (obs)...")
      obs <- tryCatch({
        .read_h5ad_dataframe(h5file, "obs", verbose)
      }, error = function(e) {
        warning("[read_h5ad] WARNING: Failed to read obs metadata: ", conditionMessage(e))
        NULL
      })
      if (!is.null(obs) && verbose) {
        message("[read_h5ad]   Cell metadata: ", nrow(obs), " cells, ", ncol(obs), " columns")
      }
    }

    # Read gene metadata (var)
    var <- NULL
    if (h5file$exists(var_path)) {
      if (verbose) message("[read_h5ad]   Reading gene metadata (var)...")
      var <- tryCatch({
        .read_h5ad_dataframe(h5file, var_path, verbose)
      }, error = function(e) {
        warning("[read_h5ad] WARNING: Failed to read var metadata: ", conditionMessage(e))
        NULL
      })
      if (!is.null(var) && verbose) {
        message("[read_h5ad]   Gene metadata: ", nrow(var), " genes, ", ncol(var), " columns")
      }
    }

    # Set dimension names
    .set_h5ad_dimnames(expr_matrix, obs, var, verbose)

    if (verbose) {
      message("[read_h5ad] SUCCESS: H5AD file read complete")
      message("[read_h5ad] Final matrix: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")
    }

    return(list(
      X = expr_matrix,
      obs = obs,
      var = var
    ))
  } else {
    # Just set dimension names from _index datasets
    tryCatch({
      if (h5file$exists("var/_index")) {
        gene_names <- h5file[["var/_index"]]$read()
        if (length(gene_names) == nrow(expr_matrix)) {
          rownames(expr_matrix) <- as.character(gene_names)
          if (verbose) message("[read_h5ad] Set ", length(gene_names), " gene names")
        }
      }
      if (h5file$exists("obs/_index")) {
        cell_names <- h5file[["obs/_index"]]$read()
        if (length(cell_names) == ncol(expr_matrix)) {
          colnames(expr_matrix) <- as.character(cell_names)
          if (verbose) message("[read_h5ad] Set ", length(cell_names), " cell names")
        }
      }
    }, error = function(e) {
      if (verbose) warning("[read_h5ad] WARNING: Could not read dimension names: ", conditionMessage(e))
    })

    if (verbose) {
      message("[read_h5ad] SUCCESS: H5AD file read complete")
      message("[read_h5ad] Final matrix: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")
    }

    return(expr_matrix)
  }
}


#' Set dimension names for H5AD matrix
#' @keywords internal
.set_h5ad_dimnames <- function(expr_matrix, obs, var, verbose) {
  # Set row names (genes)
  if (!is.null(var) && nrow(var) == nrow(expr_matrix)) {
    if ("_index" %in% colnames(var)) {
      rownames(expr_matrix) <- as.character(var$`_index`)
      if (verbose) message("[read_h5ad]   Set gene names from var$`_index`")
    } else if (!is.null(rownames(var))) {
      rownames(expr_matrix) <- rownames(var)
      if (verbose) message("[read_h5ad]   Set gene names from var rownames")
    }
  }

  # Set column names (cells)
  if (!is.null(obs) && nrow(obs) == ncol(expr_matrix)) {
    if ("_index" %in% colnames(obs)) {
      colnames(expr_matrix) <- as.character(obs$`_index`)
      if (verbose) message("[read_h5ad]   Set cell names from obs$`_index`")
    } else if (!is.null(rownames(obs))) {
      colnames(expr_matrix) <- rownames(obs)
      if (verbose) message("[read_h5ad]   Set cell names from obs rownames")
    }
  }
}


#' Read expression matrix from H5AD file
#' @keywords internal
.read_h5ad_matrix <- function(h5file, matrix_path, verbose = TRUE) {

  matrix_group <- h5file[[matrix_path]]

  # Check if sparse matrix (CSR/CSC format)
  if (inherits(matrix_group, "H5Group")) {
    if (matrix_group$exists("data") &&
        matrix_group$exists("indices") &&
        matrix_group$exists("indptr")) {

      if (verbose) message("[read_h5ad]   Matrix format: Sparse (CSR)")

      # Read sparse matrix components
      data <- tryCatch(matrix_group[["data"]]$read(), error = function(e) {
        stop("Failed to read sparse matrix 'data' component: ", conditionMessage(e))
      })
      indices <- tryCatch(matrix_group[["indices"]]$read(), error = function(e) {
        stop("Failed to read sparse matrix 'indices' component: ", conditionMessage(e))
      })
      indptr <- tryCatch(matrix_group[["indptr"]]$read(), error = function(e) {
        stop("Failed to read sparse matrix 'indptr' component: ", conditionMessage(e))
      })

      # Get shape
      shape <- tryCatch(matrix_group$attr_open("shape")$read(), error = function(e) {
        stop("Failed to read matrix shape attribute: ", conditionMessage(e))
      })
      nrows <- shape[1]
      ncols <- shape[2]

      if (verbose) {
        message("[read_h5ad]   Sparse matrix shape: ", nrows, " x ", ncols)
        message("[read_h5ad]   Non-zero elements: ", length(data))
      }

      # Convert from 0-based to 1-based indexing
      indices <- as.integer(indices) + 1L
      indptr <- as.integer(indptr) + 1L

      # Build row indices from indptr (CSR format)
      i_indices <- integer(length(data))
      counter <- 1
      for (row in 1:(length(indptr) - 1)) {
        start_idx <- indptr[row]
        end_idx <- indptr[row + 1] - 1
        if (end_idx >= start_idx) {
          num_elements <- end_idx - start_idx + 1
          i_indices[counter:(counter + num_elements - 1)] <- row
          counter <- counter + num_elements
        }
      }

      # Create sparse matrix
      sparse_mat <- tryCatch({
        Matrix::sparseMatrix(
          i = i_indices,
          j = indices,
          x = as.numeric(data),
          dims = c(nrows, ncols)
        )
      }, error = function(e) {
        stop("Failed to construct sparse matrix: ", conditionMessage(e))
      })

      return(sparse_mat)

    } else {
      stop("Unknown sparse matrix format. Expected data/indices/indptr components.")
    }
  } else {
    # Dense matrix
    if (verbose) message("[read_h5ad]   Matrix format: Dense")
    dense_mat <- tryCatch(matrix_group$read(), error = function(e) {
      stop("Failed to read dense matrix: ", conditionMessage(e))
    })

    if (verbose) message("[read_h5ad]   Converting to sparse format...")
    return(Matrix::Matrix(dense_mat, sparse = TRUE))
  }
}


#' Read data frame from H5AD file
#' @keywords internal
.read_h5ad_dataframe <- function(h5file, df_path, verbose = TRUE) {

  if (!h5file$exists(df_path)) {
    if (verbose) message("[read_h5ad]   Path '", df_path, "' not found, skipping")
    return(NULL)
  }

  df_group <- h5file[[df_path]]

  # Read _index (row names)
  row_names <- NULL
  if (df_group$exists("_index")) {
    row_names <- tryCatch({
      as.character(df_group[["_index"]]$read())
    }, error = function(e) {
      if (verbose) warning("[read_h5ad]   Could not read _index: ", conditionMessage(e))
      NULL
    })
  }

  # Get column names
  all_items <- tryCatch(df_group$ls()$name, error = function(e) character(0))
  col_names <- setdiff(all_items, "_index")

  if (length(col_names) == 0) {
    if (verbose) message("[read_h5ad]   No columns found in ", df_path)
    return(data.frame())
  }

  # Read columns
  df_list <- list()
  for (col_name in col_names) {
    tryCatch({
      col_group <- df_group[[col_name]]

      # Handle categorical data
      if (inherits(col_group, "H5Group") && col_group$exists("categories")) {
        categories <- as.character(col_group[["categories"]]$read())
        codes <- as.integer(col_group[["codes"]]$read())
        df_list[[col_name]] <- factor(categories[codes + 1], levels = categories)
      } else {
        col_data <- col_group$read()
        df_list[[col_name]] <- col_data
      }
    }, error = function(e) {
      if (verbose) warning("[read_h5ad]   Could not read column '", col_name, "': ", conditionMessage(e))
    })
  }

  # Create data frame
  if (length(df_list) == 0) {
    return(data.frame())
  }

  df <- tryCatch({
    as.data.frame(df_list, stringsAsFactors = FALSE)
  }, error = function(e) {
    if (verbose) warning("[read_h5ad]   Error creating data frame: ", conditionMessage(e))
    return(data.frame())
  })

  # Set row names
  if (!is.null(row_names) && nrow(df) == length(row_names)) {
    rownames(df) <- row_names
  }

  return(df)
}


#' Read generic H5 file and convert to sparse matrix
#'
#' Reads a generic HDF5 file containing a gene expression matrix and converts
#' it to a sparse dgCMatrix suitable for use with gedi R6 object.
#'
#' @param file_path Character. Path to the H5 file.
#' @param matrix_path Character. Path to the expression matrix within the H5 file.
#' @param gene_names_path Character. Path to gene names. Set to NULL to skip.
#' @param cell_names_path Character. Path to cell names. Set to NULL to skip.
#' @param transpose Logical. Transpose to genes x cells format.
#' @param as_sparse Logical. Convert to sparse matrix.
#' @param verbose Logical. Print progress messages.
#'
#' @return Sparse matrix (dgCMatrix) with genes as rows and cells as columns.
#'
#' @examples
#' \dontrun{
#' # Read 10X Genomics H5 file
#' expr_matrix <- read_h5("filtered_feature_bc_matrix.h5")
#'
#' # Use with gedi R6 object
#' gedi_obj <- gedi$new(n_sample = 2)
#' gedi_obj$setup(Y = expr_matrix, sample_id = samples, K = 10)
#' }
#'
#' @importFrom hdf5r H5File
#' @importFrom Matrix sparseMatrix t Matrix
#' @export
read_h5 <- function(file_path,
                    matrix_path = "matrix",
                    gene_names_path = "features/name",
                    cell_names_path = "barcodes",
                    transpose = TRUE,
                    as_sparse = TRUE,
                    verbose = TRUE) {

  if (verbose) message("[read_h5] Starting H5 file reading process...")

  # Check dependencies
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("[read_h5] ERROR: Package 'hdf5r' required.\n",
         "  Install with: install.packages('hdf5r')\n",
         "  Or use: gedi::install_h5_dependencies()")
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("[read_h5] ERROR: Package 'Matrix' required.\n",
         "  Install with: install.packages('Matrix')")
  }

  # Validate file
  if (!file.exists(file_path)) {
    stop("[read_h5] ERROR: File not found: ", file_path)
  }

  file_size_mb <- file.info(file_path)$size / 1024^2
  if (verbose) message("[read_h5] File size: ", round(file_size_mb, 2), " MB")

  # Open H5 file
  if (verbose) message("[read_h5] Opening H5 file...")
  h5file <- tryCatch({
    hdf5r::H5File$new(file_path, mode = "r")
  }, error = function(e) {
    stop("[read_h5] ERROR: Failed to open H5 file.\n",
         "  File: ", file_path, "\n",
         "  Error: ", conditionMessage(e))
  })

  on.exit({
    if (!is.null(h5file)) {
      h5file$close()
      if (verbose) message("[read_h5] H5 file closed successfully")
    }
  }, add = TRUE)

  # Check matrix path
  if (!h5file$exists(matrix_path)) {
    available_paths <- tryCatch(h5file$ls(recursive = TRUE)$name, error = function(e) character(0))
    stop("[read_h5] ERROR: Matrix path '", matrix_path, "' not found.\n",
         "  Available paths (first 20): ", paste(head(available_paths, 20), collapse = ", "), "\n",
         "  Use list_h5_structure('", file_path, "') to explore the file.")
  }

  # Read matrix
  if (verbose) message("[read_h5] Reading expression matrix from: ", matrix_path)
  expr_matrix <- tryCatch({
    .read_h5_matrix(h5file, matrix_path, as_sparse, verbose)
  }, error = function(e) {
    stop("[read_h5] ERROR: Failed to read matrix.\n",
         "  Matrix path: ", matrix_path, "\n",
         "  Error: ", conditionMessage(e))
  })

  if (verbose) {
    message("[read_h5] Matrix dimensions (before transpose): ", nrow(expr_matrix), " x ", ncol(expr_matrix))
    message("[read_h5] Matrix class: ", class(expr_matrix)[1])
    if (inherits(expr_matrix, "sparseMatrix")) {
      message("[read_h5] Sparsity: ",
              round(100 * (1 - Matrix::nnzero(expr_matrix) / (nrow(expr_matrix) * ncol(expr_matrix))), 2), "%")
    }
  }

  # Transpose
  if (transpose) {
    if (verbose) message("[read_h5] Transposing to genes x cells format...")
    expr_matrix <- if (as_sparse) Matrix::t(expr_matrix) else t(expr_matrix)
    if (verbose) message("[read_h5] New dimensions: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")
  }

  # Read gene names
  if (!is.null(gene_names_path) && h5file$exists(gene_names_path)) {
    if (verbose) message("[read_h5] Reading gene names from: ", gene_names_path)
    gene_names <- tryCatch({
      as.character(h5file[[gene_names_path]]$read())
    }, error = function(e) {
      if (verbose) warning("[read_h5] WARNING: Could not read gene names: ", conditionMessage(e))
      NULL
    })

    if (!is.null(gene_names)) {
      if (length(gene_names) == nrow(expr_matrix)) {
        rownames(expr_matrix) <- gene_names
        if (verbose) message("[read_h5] Set ", length(gene_names), " gene names")
      } else {
        warning("[read_h5] WARNING: Gene names count (", length(gene_names),
                ") != matrix rows (", nrow(expr_matrix), ")")
      }
    }
  } else if (!is.null(gene_names_path) && verbose) {
    message("[read_h5] Gene names path '", gene_names_path, "' not found, skipping")
  }

  # Read cell names
  if (!is.null(cell_names_path) && h5file$exists(cell_names_path)) {
    if (verbose) message("[read_h5] Reading cell names from: ", cell_names_path)
    cell_names <- tryCatch({
      as.character(h5file[[cell_names_path]]$read())
    }, error = function(e) {
      if (verbose) warning("[read_h5] WARNING: Could not read cell names: ", conditionMessage(e))
      NULL
    })

    if (!is.null(cell_names)) {
      if (length(cell_names) == ncol(expr_matrix)) {
        colnames(expr_matrix) <- cell_names
        if (verbose) message("[read_h5] Set ", length(cell_names), " cell names")
      } else {
        warning("[read_h5] WARNING: Cell names count (", length(cell_names),
                ") != matrix columns (", ncol(expr_matrix), ")")
      }
    }
  } else if (!is.null(cell_names_path) && verbose) {
    message("[read_h5] Cell names path '", cell_names_path, "' not found, skipping")
  }

  if (verbose) {
    message("[read_h5] SUCCESS: H5 file read complete")
    message("[read_h5] Final matrix: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")
  }

  return(expr_matrix)
}


#' Read matrix from H5 file
#' @keywords internal
.read_h5_matrix <- function(h5file, matrix_path, as_sparse = TRUE, verbose = TRUE) {

  matrix_group <- h5file[[matrix_path]]

  # Check for 10X sparse matrix format
  if (inherits(matrix_group, "H5Group")) {
    if (matrix_group$exists("data") &&
        matrix_group$exists("indices") &&
        matrix_group$exists("indptr")) {

      if (verbose) message("[read_h5]   Matrix format: 10X Sparse (CSC)")

      # Read components
      data <- tryCatch(matrix_group[["data"]]$read(), error = function(e) {
        stop("Failed to read 'data': ", conditionMessage(e))
      })
      indices <- tryCatch(matrix_group[["indices"]]$read(), error = function(e) {
        stop("Failed to read 'indices': ", conditionMessage(e))
      })
      indptr <- tryCatch(matrix_group[["indptr"]]$read(), error = function(e) {
        stop("Failed to read 'indptr': ", conditionMessage(e))
      })
      shape <- tryCatch(matrix_group[["shape"]]$read(), error = function(e) {
        stop("Failed to read 'shape': ", conditionMessage(e))
      })

      nrows <- shape[1]
      ncols <- shape[2]

      if (verbose) {
        message("[read_h5]   Sparse matrix shape: ", nrows, " x ", ncols)
        message("[read_h5]   Non-zero elements: ", length(data))
      }

      # Convert 0-based to 1-based indexing
      indices <- as.integer(indices) + 1L
      indptr <- as.integer(indptr) + 1L

      # Build column indices from indptr (CSC format)
      j_indices <- integer(length(data))
      counter <- 1
      for (col in 1:(length(indptr) - 1)) {
        start_idx <- indptr[col]
        end_idx <- indptr[col + 1] - 1
        if (end_idx >= start_idx) {
          num_elements <- end_idx - start_idx + 1
          j_indices[counter:(counter + num_elements - 1)] <- col
          counter <- counter + num_elements
        }
      }

      # Create sparse matrix
      sparse_mat <- tryCatch({
        Matrix::sparseMatrix(
          i = indices,
          j = j_indices,
          x = as.numeric(data),
          dims = c(nrows, ncols)
        )
      }, error = function(e) {
        stop("Failed to construct sparse matrix: ", conditionMessage(e))
      })

      return(sparse_mat)

    } else {
      stop("Unknown H5 matrix format. Expected data/indices/indptr for sparse matrix.")
    }
  } else {
    # Dense matrix
    if (verbose) message("[read_h5]   Matrix format: Dense")
    dense_mat <- tryCatch(matrix_group$read(), error = function(e) {
      stop("Failed to read dense matrix: ", conditionMessage(e))
    })

    if (as_sparse) {
      if (verbose) message("[read_h5]   Converting to sparse format...")
      return(Matrix::Matrix(dense_mat, sparse = TRUE))
    } else {
      return(dense_mat)
    }
  }
}


#' List structure of H5 or H5AD file
#'
#' Helper function to explore H5/H5AD file structure.
#'
#' @param file_path Character. Path to the H5 or H5AD file.
#' @param recursive Logical. List all nested groups. Default TRUE.
#'
#' @return data.frame with file structure information
#'
#' @examples
#' \dontrun{
#' list_h5_structure("data.h5ad")
#' list_h5_structure("filtered_feature_bc_matrix.h5")
#' }
#'
#' @importFrom hdf5r H5File
#' @export
list_h5_structure <- function(file_path, recursive = TRUE) {

  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("[list_h5_structure] ERROR: Package 'hdf5r' required.\n",
         "  Install with: install.packages('hdf5r')")
  }

  if (!file.exists(file_path)) {
    stop("[list_h5_structure] ERROR: File not found: ", file_path)
  }

  h5file <- tryCatch({
    hdf5r::H5File$new(file_path, mode = "r")
  }, error = function(e) {
    stop("[list_h5_structure] ERROR: Failed to open file.\n",
         "  File: ", file_path, "\n",
         "  Error: ", conditionMessage(e))
  })

  on.exit(h5file$close(), add = TRUE)

  structure_info <- tryCatch({
    h5file$ls(recursive = recursive)
  }, error = function(e) {
    stop("[list_h5_structure] ERROR: Failed to list file contents.\n",
         "  Error: ", conditionMessage(e))
  })

  cat("\n")
  cat("Structure of:", basename(file_path), "\n")
  cat(paste(rep("=", nchar(basename(file_path)) + 14), collapse = ""), "\n")
  cat("Full path:", file_path, "\n")
  cat(paste(rep("-", nchar(file_path) + 11), collapse = ""), "\n\n")

  print(structure_info)
  cat("\n")

  return(invisible(structure_info))
}
