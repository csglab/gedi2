# R/gedi_plot_utils.R

#' GEDI Plot Theme
#' @keywords internal
theme_gedi <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      # Legend
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = base_size * 0.9, face = "plain", hjust = 0.5),
      legend.text = ggplot2::element_text(size = base_size * 0.8),
      legend.key.height = ggplot2::unit(0.4, "cm"),
      legend.key.width = ggplot2::unit(3, "cm"),
      legend.margin = ggplot2::margin(t = 5, b = 5),
      
      # Title
      plot.title = ggplot2::element_text(size = base_size * 1.2, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.95, hjust = 0),
      
      # Facets
      strip.text = ggplot2::element_text(size = base_size * 0.95, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "gray90", color = "black", linewidth = 0.3)
    )
}


#' Diverging Color Scale for GEDI Plots
#' @keywords internal
scale_color_gedi_diverging <- function(limits = NULL, name = "Value", ...) {
  ggplot2::scale_color_gradientn(
    colors = c("blue", "lightgrey", "red"),
    values = scales::rescale(c(-1, 0, 1)),
    limits = limits,
    oob = scales::squish,
    name = name,
    guide = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = ggplot2::unit(4, "cm"),
      barheight = ggplot2::unit(0.4, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.5,
      ticks.colour = "black",
      ticks.linewidth = 0.5
    ),
    ...
  )
}


#' Fill Scale for GEDI Plots
#' @keywords internal
scale_fill_gedi_diverging <- function(limits = NULL, name = "Value", ...) {
  ggplot2::scale_fill_gradientn(
    colors = c("blue", "lightgrey", "red"),
    values = scales::rescale(c(-1, 0, 1)),
    limits = limits,
    oob = scales::squish,
    name = name,
    guide = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = ggplot2::unit(4, "cm"),
      barheight = ggplot2::unit(0.4, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.5,
      ticks.colour = "black",
      ticks.linewidth = 0.5
    ),
    ...
  )
}


#' Discrete Color Palette for GEDI
#' @keywords internal
scale_color_gedi_discrete <- function(name = "Group", ...) {
  ggplot2::scale_color_manual(
    values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
               "#FFFF33", "#A65628", "#F781BF", "#999999"),
    name = name,
    ...
  )
}


#' Compute Color Limits from Data
#' @keywords internal
compute_color_limits <- function(values, symmetric = TRUE, quantile = 0.99) {
  if (symmetric) {
    lim <- stats::quantile(abs(values), quantile, na.rm = TRUE)
    c(-lim, lim)
  } else {
    c(min(values, na.rm = TRUE), max(values, na.rm = TRUE))
  }
}