# ==============================================================================
# Seurat Integration Functions
# Convert between Seurat and GEDI objects
# ==============================================================================

#' Convert Seurat Object to GEDI Model
#'
#' @description
#' Extracts count data from a Seurat object and creates a GEDI model.
#' Automatically validates that the data contains raw counts (not normalized).
#' Handles both Seurat v4 and v5, including split layers in v5.
#'
#' @param seurat_object Seurat object
#' @param assay Character, which assay to use (default: "RNA")
#' @param slot Character, which slot/layer to extract (default: "counts").
#'   For Seurat v5 with split layers (e.g., counts.CTRL, counts.STIM),
#'   this will automatically detect and combine all matching layers.
#' @param sample_column Character, column name in meta.data for sample labels
#'   (default: "orig.ident")
#' @param subset_samples Character vector, subset to specific samples (default: NULL = all)
#' @param K Integer, number of latent factors (default: 10)
#' @param mode Character, normalization mode: "Bl2" or "Bsphere" (default: "Bl2")
#' @param C Gene-level prior matrix (genes × pathways) (default: NULL)
#' @param H Sample-level covariate matrix (covariates × samples) (default: NULL)
#' @param validate_counts Logical, whether to validate data appears to be counts
#'   (default: TRUE)
#' @param ... Additional arguments passed to CreateGEDIObject()
#'
#' @return GEDI R6 object
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Basic usage
#' gedi_model <- seurat_to_gedi(
#'   seurat_object = pbmc,
#'   sample_column = "orig.ident",
#'   K = 15
#' )
#'
#' # Train the model
#' gedi_model$train(iterations = 50)
#'
#' # With gene-level priors
#' gedi_model <- seurat_to_gedi(
#'   seurat_object = pbmc,
#'   sample_column = "sample_id",
#'   K = 15,
#'   C = pathway_matrix
#' )
#'
#' # Seurat v5 with split layers (counts.CTRL, counts.STIM, etc.)
#' # Automatically detects and combines split layers
#' gedi_model <- seurat_to_gedi(
#'   seurat_object = ifnb,  # Has counts.CTRL and counts.STIM
#'   assay = "RNA",
#'   slot = "counts",  # Will find and join counts.CTRL + counts.STIM
#'   sample_column = "stim",
#'   K = 15
#' )
#' }
#'
#' @export
seurat_to_gedi <- function(seurat_object,
                            assay = "RNA",
                            slot = "counts",
                            sample_column = "orig.ident",
                            subset_samples = NULL,
                            K = 10,
                            mode = "Bl2",
                            C = NULL,
                            H = NULL,
                            validate_counts = TRUE,
                            ...) {

  # Check if Seurat is available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Install with: install.packages('Seurat')",
         call. = FALSE)
  }

  # Validate inputs
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object", call. = FALSE)
  }

  if (!assay %in% names(seurat_object@assays)) {
    stop("Assay '", assay, "' not found. Available: ",
         paste(names(seurat_object@assays), collapse = ", "), call. = FALSE)
  }

  # Extract count matrix (handle both Seurat v4 and v5)
  # Seurat v5 can have split layers like counts.CTRL, counts.STIM
  assay_obj <- seurat_object@assays[[assay]]

  # Detect if we have multiple layers for the requested slot (Seurat v5 split layers)
  if (methods::is(assay_obj, "Assay5")) {
    # Seurat v5 - check for split layers
    layer_names <- names(assay_obj@layers)
    slot_pattern <- paste0("^", slot, "\\.")  # e.g., "^counts\\."

    # Find all layers matching the slot pattern
    matching_layers <- grep(slot_pattern, layer_names, value = TRUE)

    if (length(matching_layers) > 1) {
      # Multiple split layers found - need to join them
      message("Detected ", length(matching_layers), " split layers: ",
              paste(matching_layers, collapse = ", "))
      message("Joining layers with do.call(cbind)...")

      # Extract each layer and combine
      layer_matrices <- lapply(matching_layers, function(layer_name) {
        Seurat::GetAssayData(seurat_object, assay = assay, layer = layer_name)
      })

      # Combine all layers column-wise
      M <- do.call(cbind, layer_matrices)

    } else if (length(matching_layers) == 1) {
      # Single layer with pattern name (e.g., just "counts.sample1")
      M <- Seurat::GetAssayData(seurat_object, assay = assay, layer = matching_layers[1])

    } else {
      # No pattern match - try direct slot name
      M <- tryCatch(
        Seurat::GetAssayData(seurat_object, assay = assay, layer = slot),
        error = function(e) {
          stop("No layers found matching '", slot, "' in Seurat v5 assay. ",
               "Available layers: ", paste(layer_names, collapse = ", "), call. = FALSE)
        }
      )
    }

    if (slot == "data") {
      warning("Using 'data' layer instead of 'counts'. ",
              "GEDI expects raw counts for proper modeling.", call. = FALSE)
    }

  } else {
    # Seurat v4 or earlier - use traditional slot access
    if (slot == "counts") {
      M <- tryCatch(
        Seurat::GetAssayData(seurat_object, assay = assay, layer = "counts"),
        error = function(e) {
          Seurat::GetAssayData(seurat_object, assay = assay, slot = "counts")
        }
      )
    } else if (slot == "data") {
      M <- tryCatch(
        Seurat::GetAssayData(seurat_object, assay = assay, layer = "data"),
        error = function(e) {
          Seurat::GetAssayData(seurat_object, assay = assay, slot = "data")
        }
      )
      warning("Using 'data' slot instead of 'counts'. ",
              "GEDI expects raw counts for proper modeling.", call. = FALSE)
    } else {
      stop("slot must be 'counts' or 'data'", call. = FALSE)
    }
  }

  # Validate that data appears to be counts (uses deterministic check to avoid RNG)
  if (validate_counts && slot == "counts") {
    # Use deterministic sampling: take evenly spaced indices instead of random sample
    n_total <- length(M@x)  # Number of non-zero elements in sparse matrix
    if (n_total > 0) {
      # Sample up to 10000 elements evenly spaced
      n_sample <- min(10000, n_total)
      sample_indices <- as.integer(seq(1, n_total, length.out = n_sample))
      sample_values <- M@x[sample_indices]

      # Check if values are integers or very close to integers
      if (length(sample_values) > 0) {
        decimal_parts <- abs(sample_values - round(sample_values))
        non_integer_frac <- mean(decimal_parts > 1e-6)

        if (non_integer_frac > 0.01) {
          warning("Data does not appear to be raw counts (contains many non-integer values). ",
                  "GEDI expects raw count data. Set validate_counts=FALSE to skip this check.",
                  call. = FALSE)
        }
      }
    }
  }

  # Extract sample labels
  if (!sample_column %in% colnames(seurat_object@meta.data)) {
    stop("Column '", sample_column, "' not found in meta.data. Available: ",
         paste(head(colnames(seurat_object@meta.data), 10), collapse = ", "),
         call. = FALSE)
  }

  Samples <- as.character(seurat_object@meta.data[[sample_column]])

  # Subset to specific samples if requested
  if (!is.null(subset_samples)) {
    keep_cells <- Samples %in% subset_samples
    if (sum(keep_cells) == 0) {
      stop("No cells found for samples: ", paste(subset_samples, collapse = ", "),
           call. = FALSE)
    }

    M <- M[, keep_cells]
    Samples <- Samples[keep_cells]
    message("Subset to ", sum(keep_cells), " cells from ",
            length(subset_samples), " samples")
  }

  # Extract metadata (exclude sample_column to avoid duplication)
  colData <- seurat_object@meta.data
  if (!is.null(subset_samples)) {
    colData <- colData[keep_cells, , drop = FALSE]
  }

  # Remove sample column from colData (already in Samples)
  colData[[sample_column]] <- NULL

  # Create GEDI object
  message("Creating GEDI object from Seurat data...")
  message("  Genes: ", nrow(M))
  message("  Cells: ", ncol(M))
  message("  Samples: ", length(unique(Samples)))

  gedi_model <- CreateGEDIObject(
    Samples = Samples,
    M = M,
    colData = colData,
    K = K,
    mode = mode,
    C = C,
    H = H,
    ...
  )

  message("GEDI object created successfully!")

  return(gedi_model)
}


#' Convert GEDI Model to Seurat Object
#'
#' @description
#' Creates a Seurat object from a trained GEDI model, including imputed data,
#' projections, and embeddings.
#'
#' @param model GEDI model object (trained)
#' @param M Original count matrix (optional). If not provided, will use
#'   back-transformed imputed values as approximate counts.
#' @param project Character, project name for Seurat object (default: "GEDI")
#' @param assay Character, name for the main assay (default: "RNA")
#' @param use_imputed Logical, whether to add imputed data as separate assay
#'   (default: TRUE)
#' @param add_projections Logical, whether to add ZDB and DB projections as
#'   separate assays (default: TRUE)
#' @param add_embeddings Logical, whether to add UMAP and PCA embeddings if
#'   available (default: TRUE)
#' @param min_cells Integer, filter genes with counts in < min_cells (default: 0)
#' @param min_features Integer, filter cells with < min_features genes (default: 0)
#'
#' @return Seurat object with:
#'   \itemize{
#'     \item RNA assay: Original or back-transformed counts
#'     \item imputed assay: GEDI imputed expression (if use_imputed = TRUE)
#'     \item ZDB/DB assays: GEDI projections (if add_projections = TRUE)
#'     \item umap/pca reductions: Embeddings (if add_embeddings = TRUE and cached)
#'     \item meta.data: Sample labels and colData from GEDI model
#'   }
#'
#' @examples
#' \dontrun{
#' # Train GEDI model
#' gedi_model <- seurat_to_gedi(pbmc, K = 15)
#' gedi_model$train(iterations = 50)
#'
#' # Convert back to Seurat with all data
#' seurat_obj <- gedi_to_seurat(
#'   gedi_model,
#'   M = original_counts,
#'   use_imputed = TRUE,
#'   add_projections = TRUE,
#'   add_embeddings = TRUE
#' )
#'
#' # Now can use Seurat functions
#' library(Seurat)
#' DimPlot(seurat_obj, reduction = "umap", group.by = "sample")
#' FeaturePlot(seurat_obj, features = "CD3D", assay = "imputed")
#' }
#'
#' @export
gedi_to_seurat <- function(model,
                            M = NULL,
                            project = "GEDI",
                            assay = "RNA",
                            use_imputed = TRUE,
                            add_projections = TRUE,
                            add_embeddings = TRUE,
                            min_cells = 0,
                            min_features = 0) {

  # Check if Seurat is available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Install with: install.packages('Seurat')",
         call. = FALSE)
  }

  # Validate model is trained
  if (is.null(model$params)) {
    stop("Model not trained. Run model$train() first.", call. = FALSE)
  }

  message("Converting GEDI model to Seurat object...")

  # Get metadata
  gene_ids <- model$metadata$geneIDs
  cell_ids <- model$metadata$cellIDs
  samples <- model$metadata$Samples

  # Prepare count matrix
  if (is.null(M)) {
    message("  No original counts provided, using back-transformed imputed values...")
    # Get imputed log-expression and back-transform to counts (approximate)
    Y_imputed <- model$imputed$Y()
    M <- Matrix::Matrix(expm1(Y_imputed), sparse = TRUE)
  } else {
    message("  Using provided count matrix...")
    # Ensure M has correct dimensions and ordering
    if (!identical(dim(M), c(length(gene_ids), length(cell_ids)))) {
      stop("M dimensions don't match model (expected ", length(gene_ids),
           " × ", length(cell_ids), ")", call. = FALSE)
    }
  }

  # Set dimnames
  rownames(M) <- gene_ids
  colnames(M) <- cell_ids

  # Create Seurat object
  message("  Creating Seurat object...")
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = M,
    project = project,
    assay = assay,
    min.cells = min_cells,
    min.features = min_features
  )

  # Add sample information to meta.data
  seurat_obj@meta.data$Sample <- samples

  # Add colData if available
  if (!is.null(model$metadata$colData)) {
    colData <- model$metadata$colData
    # Match to cells in Seurat object (in case filtering occurred)
    common_cells <- intersect(cell_ids, colnames(seurat_obj))
    if (length(common_cells) > 0) {
      colData_matched <- colData[match(common_cells, cell_ids), , drop = FALSE]
      rownames(colData_matched) <- common_cells

      # Add each column to meta.data
      for (col in colnames(colData_matched)) {
        seurat_obj@meta.data[[col]] <- colData_matched[[col]]
      }
    }
  }

  # Add imputed data as separate assay
  if (use_imputed) {
    message("  Adding imputed assay...")
    Y_imputed <- model$imputed$Y()
    rownames(Y_imputed) <- gene_ids
    colnames(Y_imputed) <- cell_ids

    # Filter to cells in Seurat object
    common_cells <- intersect(cell_ids, colnames(seurat_obj))
    Y_imputed <- Y_imputed[, common_cells, drop = FALSE]

    # Create imputed assay (use log-transformed data as 'data' slot)
    imputed_assay <- Seurat::CreateAssayObject(
      data = Matrix::Matrix(Y_imputed, sparse = TRUE)
    )
    seurat_obj[["imputed"]] <- imputed_assay
  }

  # Add projections as separate assays
  if (add_projections) {
    message("  Adding projection assays...")

    # ZDB projection (J genes × N cells)
    ZDB <- model$projections$ZDB
    rownames(ZDB) <- gene_ids
    colnames(ZDB) <- cell_ids
    common_cells <- intersect(cell_ids, colnames(seurat_obj))
    ZDB <- ZDB[, common_cells, drop = FALSE]

    zdb_assay <- Seurat::CreateAssayObject(
      data = Matrix::Matrix(ZDB, sparse = TRUE)
    )
    seurat_obj[["ZDB"]] <- zdb_assay

    # DB projection (K factors × N cells)
    DB <- model$projections$DB
    K <- nrow(DB)
    rownames(DB) <- paste0("LV", 1:K)
    colnames(DB) <- cell_ids
    DB <- DB[, common_cells, drop = FALSE]

    db_assay <- Seurat::CreateAssayObject(
      data = Matrix::Matrix(DB, sparse = TRUE)
    )
    seurat_obj[["DB"]] <- db_assay
  }

  # Add embeddings
  if (add_embeddings) {
    cache_status <- model$cache_status()

    # Add UMAP if cached
    if (cache_status["umap"]) {
      message("  Adding UMAP embedding...")
      umap_coords <- model$embeddings$umap()
      common_cells <- intersect(cell_ids, colnames(seurat_obj))
      umap_coords <- umap_coords[match(common_cells, cell_ids), , drop = FALSE]
      rownames(umap_coords) <- common_cells
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

      seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
        embeddings = umap_coords,
        key = "UMAP_",
        assay = assay
      )
    }

    # Add PCA if cached
    if (cache_status["pca"]) {
      message("  Adding PCA embedding...")
      pca_coords <- model$embeddings$pca
      common_cells <- intersect(cell_ids, colnames(seurat_obj))
      pca_coords <- pca_coords[match(common_cells, cell_ids), , drop = FALSE]
      rownames(pca_coords) <- common_cells
      colnames(pca_coords) <- paste0("PC_", 1:ncol(pca_coords))

      seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
        embeddings = pca_coords,
        key = "PC_",
        assay = assay
      )
    }
  }

  message("Seurat object created successfully!")
  message("  Assays: ", paste(names(seurat_obj@assays), collapse = ", "))
  if (length(seurat_obj@reductions) > 0) {
    message("  Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))
  }

  return(seurat_obj)
}
