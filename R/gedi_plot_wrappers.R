# Improved Plotting API with Smart Caching and Model-First Interface
# These wrappers provide a cleaner user experience while maintaining backwards compatibility

#' Get Embedding Coordinates with Smart Caching
#'
#' @param model GEDI model object
#' @param embedding_type Character: "umap", "pca", or a custom Nx2 matrix
#' @param dims Integer vector of length 2 for which dimensions to use (default c(1,2))
#' @param verbose Logical, print messages about computation
#' @return Nx2 matrix of embedding coordinates
#' @keywords internal
.get_embedding <- function(model, embedding_type, dims = c(1, 2), verbose = TRUE) {

  # If it's already a matrix, return it
  if (is.matrix(embedding_type) || is.data.frame(embedding_type)) {
    emb <- as.matrix(embedding_type)
    if (ncol(emb) < 2) {
      stop("Custom embedding must have at least 2 columns", call. = FALSE)
    }
    return(emb[, dims])
  }

  # Handle string embedding types
  if (!is.character(embedding_type)) {
    stop("embedding must be 'umap', 'pca', or a matrix", call. = FALSE)
  }

  embedding_type <- tolower(embedding_type)

  if (embedding_type == "umap") {
    # UMAP can be accessed as field or method
    # As field: uses cached default UMAP
    # As method with no args: also uses cached default UMAP
    umap_emb <- model$embeddings$umap
    return(umap_emb[, dims])

  } else if (embedding_type == "pca") {
    # PCA returns cached value if available
    return(model$embeddings$pca[, dims])

  } else {
    stop("embedding must be 'umap', 'pca', or a matrix", call. = FALSE)
  }
}


#' Get Color Vector from Model
#'
#' @param model GEDI model object
#' @param color_by Character: "sample", metadata column name, or gene name/index
#' @param projection Character: "zdb" or "db" for gene expression
#' @return Vector of length N with color values
#' @keywords internal
.get_color_vector <- function(model, color_by, projection = "zdb") {

  if (is.null(color_by)) {
    return(NULL)
  }

  # If it's already a vector, return it
  if (is.vector(color_by) && length(color_by) > 1) {
    N <- ncol(model$params$Bi[[1]])  # Number of cells
    if (length(color_by) != N) {
      stop("color_by vector must have length N (number of cells)", call. = FALSE)
    }
    return(color_by)
  }

  # Handle string color_by
  if (!is.character(color_by)) {
    stop("color_by must be a character or vector", call. = FALSE)
  }

  # Check if it's "sample"
  if (tolower(color_by) == "sample") {
    return(model$metadata$Samples)
  }

  # Check if it's a metadata column
  if (!is.null(model$metadata$colData)) {
    if (color_by %in% colnames(model$metadata$colData)) {
      return(model$metadata$colData[[color_by]])
    }
  }

  # Check if it's a gene name
  geneIDs <- model$metadata$geneIDs
  if (color_by %in% geneIDs) {
    gene_idx <- match(color_by, geneIDs)

    # Compute projection on-demand for this gene only
    Z <- model$params$Z
    feature_weights <- Z[gene_idx, ]  # K x 1

    if (projection == "zdb") {
      expr_values <- compute_feature_projection(
        feature_weights = feature_weights,
        D = model$params$D,
        Bi_list = model$params$Bi,
        verbose = 0
      )
    } else if (projection == "db") {
      DB <- model$projections$DB  # K x N
      expr_values <- as.vector(t(DB) %*% feature_weights)  # N x 1
    } else {
      stop("projection must be 'zdb' or 'db'", call. = FALSE)
    }

    return(expr_values)
  }

  # If we get here, color_by wasn't found
  stop("color_by '", color_by, "' not found in samples, metadata, or genes", call. = FALSE)
}


#' Plot Embedding with Improved API
#'
#' Simplified interface for plotting embeddings with automatic caching.
#' Model is the first argument, and color_by handles metadata/genes automatically.
#'
#' @param model GEDI model object (or embedding matrix for backwards compatibility)
#' @param embedding Character ("umap", "pca") or Nx2 matrix
#' @param color_by Character: "sample", metadata column, gene name, or NULL
#' @param color Vector for manual coloring (overrides color_by)
#' @param projection Character: "zdb" or "db" for gene expression projection
#' @param color_limits Numeric vector c(low, high) or NULL for auto
#' @param palette Character vector of colors for continuous scale
#' @param randomize Logical, randomize point order
#' @param point_size Numeric, size of points
#' @param alpha Numeric, transparency (0-1)
#' @param raster Logical, use rasterization for large datasets
#' @param xlab Character, x-axis label
#' @param ylab Character, y-axis label
#' @param title Character, plot title
#' @param legend_title Character, legend title
#' @param verbose Logical, print computation messages
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Simple usage with smart caching
#' plot_embedding(model, embedding = "umap", color_by = "sample")
#' plot_embedding(model, embedding = "pca", color_by = "CD3D")
#'
#' # Backwards compatible - still works with matrix
#' umap_coords <- model$embeddings$umap()
#' plot_embedding(umap_coords, color = my_values)
#' }
#'
#' @export
plot_embedding <- function(model,
                          embedding = NULL,
                          color_by = NULL,
                          color = NULL,
                          projection = "zdb",
                          color_limits = NULL,
                          palette = c("blue", "lightgrey", "red"),
                          randomize = TRUE,
                          point_size = 0.3,
                          alpha = 0.9,
                          raster = FALSE,
                          xlab = "Dim 1",
                          ylab = "Dim 2",
                          title = NULL,
                          legend_title = NULL,
                          verbose = TRUE) {

  # Backwards compatibility: if model is a matrix, use old API
  if (is.matrix(model) || is.data.frame(model)) {
    # Old API: plot_embedding(embedding_matrix, color = values)
    embedding_mat <- as.matrix(model)

    if (ncol(embedding_mat) != 2) {
      stop("embedding must have exactly 2 columns", call. = FALSE)
    }

    N <- nrow(embedding_mat)

    # Handle color
    if (is.null(color)) {
      color <- rep(1, N)
      use_color <- FALSE
    } else {
      if (length(color) != N) {
        stop("color must have same length as number of rows in embedding", call. = FALSE)
      }
      use_color <- TRUE
    }

    # Create data frame
    df <- data.frame(
      Dim1 = embedding_mat[, 1],
      Dim2 = embedding_mat[, 2],
      Color = color
    )

    # Randomize once before processing
    if (randomize) {
      rand_idx <- sample(N)
      df <- df[rand_idx, ]
    }

    # Detect if color is numeric or discrete
    is_numeric <- is.numeric(color) && length(unique(color)) > 10

    # Start plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2))

    # Add points
    if (raster && N > 100000) {
      if (!use_color) {
        p <- p + ggplot2::geom_point(size = point_size, alpha = alpha)
      } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(color = Color),
                                     size = point_size, alpha = alpha)
      }
      warning("Rasterization not implemented yet. Using regular points.",
              call. = FALSE)
    } else {
      if (!use_color) {
        p <- p + ggplot2::geom_point(size = point_size, alpha = alpha)
      } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(color = Color),
                                     size = point_size, alpha = alpha)
      }
    }

    # Add color scale
    if (use_color) {
      if (is_numeric) {
        # Continuous color scale
        if (is.null(color_limits)) {
          color_limits <- compute_color_limits(color, symmetric = TRUE, quantile = 0.99)
        }

        if (is.null(legend_title)) {
          legend_title <- "Value"
        }

        p <- p + scale_color_gedi_diverging(
          limits = color_limits,
          name = legend_title
        )
      } else {
        # Discrete color scale
        if (is.null(legend_title)) {
          legend_title <- "Group"
        }

        p <- p + scale_color_gedi_discrete(name = legend_title) +
          ggplot2::guides(color = ggplot2::guide_legend(
            override.aes = list(size = 3, alpha = 1)
          ))
      }
    }

    # Add labels and theme
    p <- p +
      ggplot2::labs(x = xlab, y = ylab, title = title) +
      theme_gedi()

    return(p)
  }

  # New API: plot_embedding(model, embedding = "umap", color_by = "sample")
  if (!inherits(model, "GEDI")) {
    stop("First argument must be a GEDI model or embedding matrix", call. = FALSE)
  }

  if (is.null(embedding)) {
    stop("embedding argument is required (e.g., 'umap', 'pca')", call. = FALSE)
  }

  # Get embedding coordinates with smart caching
  embedding_mat <- .get_embedding(model, embedding, dims = c(1, 2), verbose = verbose)

  # Get color vector if color_by is specified (overrides manual color)
  if (!is.null(color_by) && is.null(color)) {
    color <- .get_color_vector(model, color_by, projection = projection)

    # Set legend title if not specified
    if (is.null(legend_title)) {
      if (is.character(color_by)) {
        legend_title <- color_by
      }
    }
  }

  # Set axis labels based on embedding type
  if (xlab == "Dim 1" && ylab == "Dim 2") {
    emb_name <- if (is.character(embedding)) toupper(embedding) else "Embedding"
    xlab <- paste(emb_name, "1")
    ylab <- paste(emb_name, "2")
  }

  # Call the base function with extracted data
  .plot_embedding_base(
    embedding = embedding_mat,
    color = color,
    color_limits = color_limits,
    palette = palette,
    randomize = randomize,
    point_size = point_size,
    alpha = alpha,
    raster = raster,
    xlab = xlab,
    ylab = ylab,
    title = title,
    legend_title = legend_title
  )
}
