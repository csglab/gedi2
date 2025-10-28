# R/gedi_plot_embedding.R

#' Plot 2D Embedding (Base Function)
#'
#' Base function for plotting 2D embeddings with customizable coloring.
#' Supports both continuous and discrete color variables.
#'
#' @param embedding Matrix (N × 2) with x and y coordinates for each cell
#' @param color Vector of length N for coloring points, or NULL for uniform color
#' @param color_limits Numeric vector c(low, high) for color scale limits,
#'   or NULL to auto-compute from data (uses 99th percentile for symmetric limits)
#' @param palette Character vector of colors for continuous scale (length 3 for diverging)
#' @param randomize Logical, whether to randomize point order before plotting
#' @param point_size Numeric, size of points
#' @param alpha Numeric, transparency of points (0-1)
#' @param raster Logical, use rasterization for large datasets (>100k points)
#' @param xlab Character, x-axis label
#' @param ylab Character, y-axis label
#' @param title Character, plot title
#' @param legend_title Character, legend title (auto-detected if NULL)
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' umap_coords <- model$embeddings$umap()
#' feature_values <- model$projections$ZDB[1, ]  # First gene
#' .plot_embedding_base(umap_coords, color = feature_values)
#' }
#'
#' @keywords internal
.plot_embedding_base <- function(embedding,
                          color = NULL,
                          color_limits = NULL,
                          palette = c("blue", "lightgrey", "red"),
                          randomize = TRUE,
                          point_size = 0.3,
                          alpha = 0.9,
                          raster = FALSE,
                          xlab = "Dim 1",
                          ylab = "Dim 2",
                          title = NULL,
                          legend_title = NULL) {
  
  # Validate inputs
  if (!is.matrix(embedding) && !is.data.frame(embedding)) {
    stop("embedding must be a matrix or data.frame", call. = FALSE)
  }
  
  if (ncol(embedding) != 2) {
    stop("embedding must have exactly 2 columns", call. = FALSE)
  }
  
  N <- nrow(embedding)
  
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
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Color = color
  )
  
  # Randomize ONCE before any processing
  if (randomize) {
    rand_idx <- sample(N)
    df <- df[rand_idx, ]
  }
  
  # Detect if color is numeric or discrete
  is_numeric <- is.numeric(color) && length(unique(color)) > 10
  
  # Start plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2))
  
  # Add points (raster or regular)
  if (raster && N > 100000) {
    if (!use_color) {
      p <- p + ggplot2::geom_point(size = point_size, alpha = alpha)
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = Color), 
                                   size = point_size, alpha = alpha)
    }
    message("Note: Rasterization not implemented yet. Using regular points.")
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


#' Plot Multiple Features on Embedding
#'
#' Efficiently plots multiple gene features on a 2D embedding using faceting.
#' Computes projections on-the-fly without storing full ZDB matrix.
#'
#' @param model GEDI model object
#' @param features Character vector of gene names or integer indices
#' @param embedding Character specifying embedding type ("umap", "pca") or
#'   a custom N × 2 matrix
#' @param projection Character, type of projection to compute ("zdb" or "db")
#' @param color_limits Character ("global" for shared scale, "individual" for
#'   per-facet scale) or numeric vector c(low, high)
#' @param ncol Integer, number of columns in facet layout
#' @param randomize Logical, whether to randomize point order
#' @param point_size Numeric, size of points
#' @param alpha Numeric, transparency of points
#' @param title Character, plot title
#'
#' @return ggplot2 object with faceted features
#'
#' @examples
#' \dontrun{
#' plot_features(model, c("CD3D", "CD79A", "LYZ"), embedding = "umap")
#' plot_features(model, c(1, 2, 3), color_limits = "individual")
#' }
#'
#' @export
plot_features <- function(model,
                         features,
                         embedding = "umap",
                         projection = "zdb",
                         color_limits = "global",
                         ncol = NULL,
                         randomize = TRUE,
                         point_size = 0.2,
                         alpha = 0.9,
                         title = NULL) {
  
  # Validate model
  if (is.null(model$params)) {
    stop("Model not trained. Run $train() first.", call. = FALSE)
  }
  
  # Get embedding coordinates
  if (is.character(embedding)) {
    if (embedding == "umap") {
      emb_mat <- model$embeddings$umap()
    } else if (embedding == "pca") {
      emb_mat <- model$embeddings$pca[, 1:2]
    } else {
      stop("embedding must be 'umap', 'pca', or a matrix", call. = FALSE)
    }
  } else if (is.matrix(embedding) || is.data.frame(embedding)) {
    emb_mat <- as.matrix(embedding)
    if (ncol(emb_mat) < 2) {
      stop("Custom embedding must have at least 2 columns", call. = FALSE)
    }
    emb_mat <- emb_mat[, 1:2]
  } else {
    stop("embedding must be a character or matrix", call. = FALSE)
  }
  
  N <- nrow(emb_mat)
  
  # Convert feature names to indices
  geneIDs <- model$metadata$geneIDs
  if (is.character(features)) {
    feature_idx <- match(features, geneIDs)
    if (any(is.na(feature_idx))) {
      missing <- features[is.na(feature_idx)]
      stop("Features not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    feature_names <- features
  } else {
    feature_idx <- as.integer(features)
    if (any(feature_idx < 1 | feature_idx > length(geneIDs))) {
      stop("Feature indices out of range", call. = FALSE)
    }
    feature_names <- geneIDs[feature_idx]
  }
  
  F <- length(feature_idx)
  
  # Extract feature weights from Z
  Z <- model$params$Z
  feature_weights <- Z[feature_idx, , drop = FALSE]  # F × K
  feature_weights <- t(feature_weights)  # K × F for C++
  
  # Compute projections using C++
  if (projection == "zdb") {
    projections <- compute_multi_feature_projection(
      feature_weights = feature_weights,
      D = model$params$D,
      Bi_list = model$params$Bi,
      verbose = 0
    )  # Returns N × F matrix
  } else if (projection == "db") {
    # For DB: just use factor weights directly
    DB <- model$projections$DB  # K × N
    projections <- t(DB) %*% feature_weights  # N × F
  } else {
    stop("projection must be 'zdb' or 'db'", call. = FALSE)
  }
  
  # Build long-format data frame
  df_list <- vector("list", F)
  for (f in 1:F) {
    df_list[[f]] <- data.frame(
      Dim1 = emb_mat[, 1],
      Dim2 = emb_mat[, 2],
      Value = projections[, f],
      Feature = feature_names[f]
    )
  }
  df <- do.call(rbind, df_list)
  df$Feature <- factor(df$Feature, levels = feature_names)
  
  # Randomize ONCE before plotting
  if (randomize) {
    df <- df[sample(nrow(df)), ]
  }
  
  # Compute color limits
  if (is.character(color_limits)) {
    if (color_limits == "global") {
      lim <- compute_color_limits(df$Value, symmetric = TRUE, quantile = 0.99)
      use_free_scale <- FALSE
    } else if (color_limits == "individual") {
      lim <- NULL
      use_free_scale <- TRUE
    } else {
      stop("color_limits must be 'global', 'individual', or numeric vector", call. = FALSE)
    }
  } else {
    lim <- color_limits
    use_free_scale <- FALSE
  }
  
  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2, color = Value)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::facet_wrap(~ Feature, ncol = ncol)
  
  # Add color scale
  if (use_free_scale) {
    # Per-facet limits - need to compute per feature
    p <- p + scale_color_gedi_diverging(name = "Expression")
  } else {
    p <- p + scale_color_gedi_diverging(limits = lim, name = "Expression")
  }
  
  # Add labels and theme
  emb_name <- if (is.character(embedding)) toupper(embedding) else "Embedding"
  p <- p +
    ggplot2::labs(
      x = paste(emb_name, "1"),
      y = paste(emb_name, "2"),
      title = if (is.null(title)) "Feature Expression" else title
    ) +
    theme_gedi()
  
  return(p)
}


#' Plot Two-Feature Comparison
#'
#' Compares two features by computing their difference or correlation in
#' the projected space. Mathematically grounded for GEDI's log-space representation.
#'
#' @param model GEDI model object
#' @param gene1 Character or integer, first gene name or index
#' @param gene2 Character or integer, second gene name or index
#' @param comparison Character, type of comparison ("difference" or "correlation")
#' @param embedding Character specifying embedding type ("umap", "pca") or
#'   a custom N × 2 matrix
#' @param projection Character, type of projection ("zdb" or "db")
#' @param color_limits Numeric vector c(low, high) or NULL for auto-compute
#' @param randomize Logical, whether to randomize point order
#' @param point_size Numeric, size of points
#' @param alpha Numeric, transparency of points
#' @param title Character, plot title
#'
#' @return ggplot2 object
#'
#' @details
#' For comparison = "difference":
#' Computes (Z[gene1,] - Z[gene2,]) * D * B, equivalent to ZDB[gene1,] - ZDB[gene2,].
#' In log-space, this represents log(gene1/gene2) in the original count space.
#' Positive values indicate gene1 > gene2, negative indicates gene2 > gene1.
#'
#' @examples
#' \dontrun{
#' # Compare CD3D vs CD79A expression
#' plot_feature_ratio(model, "CD3D", "CD79A", comparison = "difference")
#' }
#'
#' @export
plot_feature_ratio <- function(model,
                               gene1,
                               gene2,
                               comparison = "difference",
                               embedding = "umap",
                               projection = "zdb",
                               color_limits = NULL,
                               randomize = TRUE,
                               point_size = 0.3,
                               alpha = 0.9,
                               title = NULL) {
  
  # Validate model
  if (is.null(model$params)) {
    stop("Model not trained. Run $train() first.", call. = FALSE)
  }
  
  # Validate comparison type
  if (!comparison %in% c("difference", "correlation")) {
    stop("comparison must be 'difference' or 'correlation'", call. = FALSE)
  }
  
  # Get embedding coordinates
  if (is.character(embedding)) {
    if (embedding == "umap") {
      emb_mat <- model$embeddings$umap()
    } else if (embedding == "pca") {
      emb_mat <- model$embeddings$pca[, 1:2]
    } else {
      stop("embedding must be 'umap', 'pca', or a matrix", call. = FALSE)
    }
  } else {
    emb_mat <- as.matrix(embedding[, 1:2])
  }
  
  # Convert gene names to indices
  geneIDs <- model$metadata$geneIDs
  
  if (is.character(gene1)) {
    gene1_idx <- match(gene1, geneIDs)
    if (is.na(gene1_idx)) stop("gene1 not found: ", gene1, call. = FALSE)
    gene1_name <- gene1
  } else {
    gene1_idx <- as.integer(gene1)
    gene1_name <- geneIDs[gene1_idx]
  }
  
  if (is.character(gene2)) {
    gene2_idx <- match(gene2, geneIDs)
    if (is.na(gene2_idx)) stop("gene2 not found: ", gene2, call. = FALSE)
    gene2_name <- gene2
  } else {
    gene2_idx <- as.integer(gene2)
    gene2_name <- geneIDs[gene2_idx]
  }
  
  # Compute comparison
  if (comparison == "difference") {
    # Compute (Z[gene1,] - Z[gene2,]) * D * B
    Z <- model$params$Z
    weight_diff <- Z[gene1_idx, ] - Z[gene2_idx, ]  # K × 1
    
    ratio_values <- compute_feature_projection(
      feature_weights = weight_diff,
      D = model$params$D,
      Bi_list = model$params$Bi,
      verbose = 0
    )
    
    legend_title <- paste0(gene1_name, " - ", gene2_name)
    
  } else if (comparison == "correlation") {
    # Future: 2D scatter or correlation coefficient per cell
    stop("correlation comparison not yet implemented", call. = FALSE)
  }
  
  # Compute color limits
  if (is.null(color_limits)) {
    color_limits <- compute_color_limits(ratio_values, symmetric = TRUE, quantile = 0.99)
  }
  
  # Create plot using plot_embedding
  emb_name <- if (is.character(embedding)) toupper(embedding) else "Embedding"
  
  if (is.null(title)) {
    title <- paste("Feature Comparison:", gene1_name, "vs", gene2_name)
  }
  
  p <- .plot_embedding_base(
    embedding = emb_mat,
    color = ratio_values,
    color_limits = color_limits,
    randomize = randomize,
    point_size = point_size,
    alpha = alpha,
    xlab = paste(emb_name, "1"),
    ylab = paste(emb_name, "2"),
    title = title,
    legend_title = legend_title
  )
  
  return(p)
}


# R/gedi_plot_dynamics.R

#' Plot Vector Field from Dynamics Analysis
#'
#' Visualizes vector fields showing cell state transitions. Uses binned
#' aggregation for cleaner visualization without overplotting.
#'
#' @param dynamics_svd Result from model$dynamics$vector_field() or similar
#' @param color Vector of length N for coloring arrows, or NULL
#' @param alpha Vector of length N or scalar for arrow transparency
#' @param n_bins Integer, number of bins per dimension for aggregation
#' @param min_per_bin Integer, minimum observations required per bin
#' @param randomize Logical, whether to randomize data order
#' @param arrow_size Numeric, size of arrow lines
#' @param arrow_length Numeric, length of arrow heads (in cm)
#' @param arrow_color Character, color for arrows (if color is NULL)
#' @param dims Integer vector of length 2, which dimensions to plot
#' @param xlab Character, x-axis label
#' @param ylab Character, y-axis label
#' @param title Character, plot title
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Compute vector field
#' vf <- model$dynamics$vector_field(
#'   start.cond = c(0, 0),
#'   end.cond = c(1, 0)
#' )
#' 
#' # Plot with default settings
#' plot_vector_field(vf)
#' 
#' # Color by feature
#' feature_vals <- model$projections$ZDB[1, ]  # First gene
#' plot_vector_field(vf, color = feature_vals)
#' }
#'
#' @export
plot_vector_field <- function(dynamics_svd,
                              color = NULL,
                              alpha = 1,
                              n_bins = 50,
                              min_per_bin = 10,
                              randomize = TRUE,
                              arrow_size = 0.5,
                              arrow_length = 0.15,
                              arrow_color = "black",
                              dims = c(1, 2),
                              xlab = NULL,
                              ylab = NULL,
                              title = NULL) {
  
  # Validate input
  if (!inherits(dynamics_svd, "gedi_dynamics_svd")) {
    stop("dynamics_svd must be output from model$dynamics functions", call. = FALSE)
  }
  
  if (length(dims) != 2 || !all(dims %in% 1:ncol(dynamics_svd$v))) {
    stop("dims must be a vector of 2 valid dimension indices", call. = FALSE)
  }
  
  # Extract start and end points
  v_matrix <- dynamics_svd$v
  indices <- dynamics_svd$indices
  
  N <- length(indices$embedding)
  
  # Get start points (original embedding)
  start_idx <- indices$embedding
  Dim1 <- v_matrix[start_idx, dims[1]]
  Dim2 <- v_matrix[start_idx, dims[2]]
  
  # Get end points (vector field)
  if ("vector_field" %in% names(indices)) {
    # Standard vector field: start at 1:N, end at (N+1):(2N)
    end_idx <- indices$vector_field[indices$vector_field > N]
    To1 <- v_matrix[end_idx, dims[1]]
    To2 <- v_matrix[end_idx, dims[2]]
  } else {
    stop("dynamics_svd does not contain vector_field information", call. = FALSE)
  }
  
  # Handle color
  if (is.null(color)) {
    color <- rep(1, N)
    use_color <- FALSE
  } else {
    if (length(color) != N) {
      stop("color must have same length as number of cells", call. = FALSE)
    }
    use_color <- TRUE
  }
  
  # Handle alpha
  if (length(alpha) == 1) {
    alpha <- rep(alpha, N)
  } else if (length(alpha) != N) {
    stop("alpha must be scalar or have length N", call. = FALSE)
  }
  
  # Call C++ aggregation function
  agg_df <- aggregate_vectors(
    Dim1 = Dim1,
    Dim2 = Dim2,
    To1 = To1,
    To2 = To2,
    color = color,
    alpha = alpha,
    n_bins = n_bins,
    min_per_bin = min_per_bin
  )
  
  # Randomize if requested
  if (randomize && nrow(agg_df) > 0) {
    agg_df <- agg_df[sample(nrow(agg_df)), ]
  }
  
  # Detect if color is numeric
  is_numeric <- use_color && is.numeric(color) && length(unique(color)) > 10
  
  # Create base plot
  if (!use_color || length(unique(agg_df$Color)) == 1) {
    # Uniform color
    p <- ggplot2::ggplot(agg_df, ggplot2::aes(x = Dim1, y = Dim2,
                                               xend = To1, yend = To2,
                                               alpha = Alpha)) +
      ggplot2::geom_segment(
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow_length, "cm")),
        color = arrow_color,
        linewidth = arrow_size
      ) +
      ggplot2::scale_alpha_identity()
    
  } else if (is_numeric) {
    # Continuous color
    lim <- compute_color_limits(agg_df$Color, symmetric = TRUE, quantile = 0.99)
    
    p <- ggplot2::ggplot(agg_df, ggplot2::aes(x = Dim1, y = Dim2,
                                               xend = To1, yend = To2,
                                               color = Color, alpha = Alpha)) +
      ggplot2::geom_segment(
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow_length, "cm")),
        linewidth = arrow_size
      ) +
      scale_color_gedi_diverging(limits = lim, name = "Value") +
      ggplot2::scale_alpha_identity()
    
  } else {
    # Discrete color
    p <- ggplot2::ggplot(agg_df, ggplot2::aes(x = Dim1, y = Dim2,
                                               xend = To1, yend = To2,
                                               color = factor(Color), alpha = Alpha)) +
      ggplot2::geom_segment(
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow_length, "cm")),
        linewidth = arrow_size
      ) +
      scale_color_gedi_discrete(name = "Group") +
      ggplot2::scale_alpha_identity()
  }
  
  # Add labels
  if (is.null(xlab)) xlab <- paste("LV", dims[1])
  if (is.null(ylab)) ylab <- paste("LV", dims[2])
  if (is.null(title)) {
    if (!is.null(dynamics_svd$metadata$method)) {
      title <- paste("Vector Field:", dynamics_svd$metadata$method)
    } else {
      title <- "Vector Field"
    }
  }
  
  p <- p +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    theme_gedi()
  
  return(p)
}



# R/gedi_plot_diagnostics.R

#' Plot Dispersion Analysis
#'
#' Visualizes the relationship between expected and observed variance
#' for count data. Useful for assessing model fit quality.
#'
#' @param dispersion_df Data frame from compute_dispersion() with columns:
#'   Expected_Var, Observed_Var, Sample, n, bin
#' @param show_identity Logical, whether to show y=x identity line
#' @param point_size Numeric, size of points
#' @param alpha Numeric, transparency of points
#' @param title Character, plot title
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' disp <- model$imputed$dispersion(M)
#' plot_dispersion(disp)
#' }
#'
#' @export
plot_dispersion <- function(dispersion_df,
                           show_identity = TRUE,
                           point_size = 0.1,
                           alpha = 0.5,
                           title = "Dispersion Analysis") {
  
  # Validate input
  required_cols <- c("Expected_Var", "Observed_Var", "Sample")
  if (!all(required_cols %in% colnames(dispersion_df))) {
    stop("dispersion_df must contain columns: ", 
         paste(required_cols, collapse = ", "), call. = FALSE)
  }
  
  # Create plot
  p <- ggplot2::ggplot(dispersion_df, 
                       ggplot2::aes(x = Expected_Var, y = Observed_Var, color = Sample)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()
  
  # Add identity line
  if (show_identity) {
    p <- p + ggplot2::geom_abline(intercept = 0, slope = 1, 
                                   linetype = "dashed", color = "gray30")
  }
  
  # Add labels and theme
  p <- p +
    scale_color_gedi_discrete(name = "Sample") +
    ggplot2::labs(
      x = "Expected Variance",
      y = "Observed Variance",
      title = title
    ) +
    theme_gedi()
  
  return(p)
}


#' Plot Training Convergence
#'
#' Visualizes convergence of model parameters during training.
#' Supports multiple layout styles for different use cases.
#'
#' @param model GEDI model object
#' @param layout Character, layout style:
#'   \itemize{
#'     \item "faceted": Single plot with facet_wrap (default, best for reports)
#'     \item "separate": List of individual plots (for interactive exploration)
#'     \item "compact": Two-panel plot (global vs sample-specific)
#'   }
#' @param params Character vector, which parameters to include. Options:
#'   "Z", "A", "o", "Bi", "Qi", "oi", "si", "Rk", "Ro", "sigma2".
#'   NULL means all available.
#' @param log_scale Logical, use log10 scale for y-axis
#' @param smooth Logical, add smooth trend line
#' @param title Character, plot title
#'
#' @return ggplot2 object (for "faceted" or "compact") or list of ggplot2 objects (for "separate")
#'
#' @examples
#' \dontrun{
#' # Faceted view (default)
#' plot_convergence(model)
#' 
#' # Separate plots for detailed inspection
#' plots <- plot_convergence(model, layout = "separate")
#' plots$sigma2
#' 
#' # Compact two-panel view
#' plot_convergence(model, layout = "compact")
#' }
#'
#' @export
plot_convergence <- function(model,
                            layout = c("faceted", "separate", "compact"),
                            params = NULL,
                            log_scale = TRUE,
                            smooth = FALSE,
                            title = "Training Convergence") {
  
  layout <- match.arg(layout)
  
  # Validate model
  if (is.null(model$tracking)) {
    stop("No tracking data available. Model may not have been trained with tracking.", 
         call. = FALSE)
  }
  
  tracking <- model$tracking
  aux <- model$aux
  
  # Determine which parameters are available
  available_params <- c()
  if (!all(is.na(tracking$dZ))) available_params <- c(available_params, "Z")
  if (!all(is.na(tracking$dA))) available_params <- c(available_params, "A")
  if (!all(is.na(tracking$do))) available_params <- c(available_params, "o")
  if (!all(is.na(tracking$dRo))) available_params <- c(available_params, "Ro")
  if (!all(is.na(tracking$sigma2))) available_params <- c(available_params, "sigma2")
  
  # Check sample-specific parameters
  has_Bi <- !all(sapply(tracking$dBi, function(x) all(is.na(x))))
  has_Qi <- !all(sapply(tracking$dQi, function(x) all(is.na(x))))
  has_oi <- !all(sapply(tracking$doi, function(x) all(is.na(x))))
  has_si <- !all(sapply(tracking$dsi, function(x) all(is.na(x))))
  has_Rk <- !all(sapply(tracking$dRk, function(x) all(is.na(x))))
  
  if (has_Bi) available_params <- c(available_params, "Bi")
  if (has_Qi) available_params <- c(available_params, "Qi")
  if (has_oi) available_params <- c(available_params, "oi")
  if (has_si) available_params <- c(available_params, "si")
  if (has_Rk) available_params <- c(available_params, "Rk")
  
  # Filter params if specified
  if (!is.null(params)) {
    params <- intersect(params, available_params)
    if (length(params) == 0) {
      stop("None of the requested parameters are available in tracking data", 
           call. = FALSE)
    }
  } else {
    params <- available_params
  }
  
  # Build data frames for each parameter type
  plot_list <- list()
  
  # Global parameters (Z, A, o, Ro, sigma2)
  if (any(c("Z", "A", "o", "Ro") %in% params)) {
    global_df <- data.frame()
    l <- length(tracking$dZ)
    
    if ("Z" %in% params) {
      global_df <- rbind(global_df, data.frame(
        Iteration = 1:l,
        RMSD = tracking$dZ,
        Parameter = "Z"
      ))
    }
    if ("A" %in% params) {
      global_df <- rbind(global_df, data.frame(
        Iteration = 1:l,
        RMSD = tracking$dA,
        Parameter = "A"
      ))
    }
    if ("o" %in% params) {
      global_df <- rbind(global_df, data.frame(
        Iteration = 1:l,
        RMSD = tracking$do,
        Parameter = "o"
      ))
    }
    if ("Ro" %in% params) {
      global_df <- rbind(global_df, data.frame(
        Iteration = 1:l,
        RMSD = tracking$dRo,
        Parameter = "Ro"
      ))
    }
    
    # Remove NA rows
    global_df <- global_df[!is.na(global_df$RMSD), ]
    
    if (nrow(global_df) > 0) {
      p_global <- ggplot2::ggplot(global_df, 
                                   ggplot2::aes(x = Iteration, y = RMSD, color = Parameter)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 1)
      
      if (log_scale) p_global <- p_global + ggplot2::scale_y_log10()
      if (smooth) p_global <- p_global + ggplot2::geom_smooth(se = FALSE, linewidth = 0.5)
      
      p_global <- p_global +
        scale_color_gedi_discrete(name = "Parameter") +
        ggplot2::labs(x = "Iteration", y = "RMSD", 
                      title = "Global Parameters") +
        theme_gedi()
      
      plot_list$global <- p_global
    }
  }
  
  # Sample-specific parameters (Bi, Qi, oi, si)
  sample_params <- intersect(c("Bi", "Qi", "oi", "si"), params)
  if (length(sample_params) > 0) {
    
    for (param in sample_params) {
      param_df <- data.frame()
      l <- length(tracking[[paste0("d", param)]][[1]])
      
      for (i in 1:aux$numSamples) {
        param_df <- rbind(param_df, data.frame(
          Iteration = 1:l,
          RMSD = tracking[[paste0("d", param)]][[i]],
          Sample = model$metadata$sampleIDs[i]
        ))
      }
      
      param_df <- param_df[!is.na(param_df$RMSD), ]
      
      if (nrow(param_df) > 0) {
        p <- ggplot2::ggplot(param_df, 
                            ggplot2::aes(x = Iteration, y = RMSD, color = Sample)) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 1)
        
        if (log_scale) p <- p + ggplot2::scale_y_log10()
        if (smooth) p <- p + ggplot2::geom_smooth(se = FALSE, linewidth = 0.5)
        
        p <- p +
          scale_color_gedi_discrete(name = "Sample") +
          ggplot2::labs(x = "Iteration", y = "RMSD",
                       title = paste0("Parameter: ", param)) +
          theme_gedi()
        
        plot_list[[param]] <- p
      }
    }
  }
  
  # Rk parameters (if L > 0)
  if ("Rk" %in% params && has_Rk) {
    Rk_df <- data.frame()
    l <- length(tracking$dRk[[1]])
    
    for (k in 1:aux$K) {
      Rk_df <- rbind(Rk_df, data.frame(
        Iteration = 1:l,
        RMSD = tracking$dRk[[k]],
        Factor = paste0("LV", k)
      ))
    }
    
    # Add Ro
    Rk_df <- rbind(Rk_df, data.frame(
      Iteration = 1:l,
      RMSD = tracking$dRo,
      Factor = "o"
    ))
    
    Rk_df <- Rk_df[!is.na(Rk_df$RMSD), ]
    
    if (nrow(Rk_df) > 0) {
      p_Rk <- ggplot2::ggplot(Rk_df, 
                              ggplot2::aes(x = Iteration, y = RMSD, color = Factor)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 1)
      
      if (log_scale) p_Rk <- p_Rk + ggplot2::scale_y_log10()
      if (smooth) p_Rk <- p_Rk + ggplot2::geom_smooth(se = FALSE, linewidth = 0.5)
      
      p_Rk <- p_Rk +
        scale_color_gedi_discrete(name = "Factor") +
        ggplot2::labs(x = "Iteration", y = "RMSD",
                     title = "Covariate Effects (Rk)") +
        theme_gedi()
      
      plot_list$Rk <- p_Rk
    }
  }
  
  # sigma2
  if ("sigma2" %in% params && !all(is.na(tracking$sigma2))) {
    sigma2_df <- data.frame(
      Iteration = 1:length(tracking$sigma2),
      sigma2 = tracking$sigma2
    )
    
    p_sigma2 <- ggplot2::ggplot(sigma2_df, 
                                ggplot2::aes(x = Iteration, y = sigma2)) +
      ggplot2::geom_line(color = "#377EB8") +
      ggplot2::geom_point(size = 1, color = "#377EB8")
    
    if (log_scale) p_sigma2 <- p_sigma2 + ggplot2::scale_y_log10()
    if (smooth) p_sigma2 <- p_sigma2 + ggplot2::geom_smooth(se = FALSE, linewidth = 0.5)
    
    p_sigma2 <- p_sigma2 +
      ggplot2::labs(x = "Iteration", y = expression(sigma^2),
                   title = "Model Variance") +
      theme_gedi()
    
    plot_list$sigma2 <- p_sigma2
  }
  
  # Return based on layout
  if (layout == "separate") {
    return(plot_list)
    
  } else if (layout == "faceted") {
    # Combine all into single faceted plot
    combined_df <- data.frame()
    
    # Add global params
    if (!is.null(plot_list$global)) {
      global_data <- plot_list$global$data
      global_data$Type <- "Global"
      combined_df <- rbind(combined_df, global_data)
    }
    
    # Add sample-specific params
    for (param in sample_params) {
      if (!is.null(plot_list[[param]])) {
        param_data <- plot_list[[param]]$data
        param_data$Parameter <- param
        param_data$Type <- param
        names(param_data)[names(param_data) == "Sample"] <- "Group"
        combined_df <- rbind(combined_df, param_data)
      }
    }
    
    # Add Rk
    if (!is.null(plot_list$Rk)) {
      Rk_data <- plot_list$Rk$data
      Rk_data$Parameter <- "Rk"
      Rk_data$Type <- "Rk"
      names(Rk_data)[names(Rk_data) == "Factor"] <- "Group"
      combined_df <- rbind(combined_df, Rk_data)
    }
    
    # Add sigma2
    if (!is.null(plot_list$sigma2)) {
      sigma2_data <- plot_list$sigma2$data
      sigma2_data$RMSD <- sigma2_data$sigma2
      sigma2_data$Parameter <- "sigma2"
      sigma2_data$Type <- "sigma2"
      sigma2_data$Group <- "sigma2"
      combined_df <- rbind(combined_df, sigma2_data[, c("Iteration", "RMSD", "Parameter", "Type", "Group")])
    }
    
    if (nrow(combined_df) == 0) {
      stop("No valid tracking data to plot", call. = FALSE)
    }
    
    p <- ggplot2::ggplot(combined_df, 
                        ggplot2::aes(x = Iteration, y = RMSD, color = Group)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 0.5) +
      ggplot2::facet_wrap(~ Type, scales = "free_y")
    
    if (!log_scale) {
      p <- p + ggplot2::scale_y_continuous()
    } else {
      p <- p + ggplot2::scale_y_log10()
    }
    
    p <- p +
      scale_color_gedi_discrete(name = "Group") +
      ggplot2::labs(x = "Iteration", y = "RMSD", title = title) +
      theme_gedi()
    
    return(p)
    
  } else if (layout == "compact") {
    # Two-panel: global + sample-specific
    stop("Compact layout not yet implemented. Use 'faceted' or 'separate'.", 
         call. = FALSE)
  }
}



