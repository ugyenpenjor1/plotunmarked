
#' Plot occupancy interaction surface from an `unmarked` model
#'
#' Generates an interaction surface plot of occupancy probability for two interacting site-level covariates
#' from a fitted `unmarked` occupancy model (class `unmarkedFitOccu`).
#' This function creates a prediction grid over the range of two covariates, makes predictions in chunks
#' to improve performance, and visualises the interaction using a color gradient.
#'
#' @param model A fitted `unmarkedFitOccu` model object created using `unmarked::occu()`.
#' @param cov1 A character string indicating the first covariate involved in the interaction (x-axis).
#' @param cov2 A character string indicating the second covariate involved in the interaction (color scale).
#' @param grid_length Integer. Number of points to generate along each covariate axis (default is 400).
#' @param chunk_size Integer. Number of rows to predict per chunk (default is 50,000).
#' @param palette Character. Name of the color palette to use for the gradient. Must be compatible with `hcl.colors()`.
#'        Default is `"cividis"`. If `NULL`, uses `"cividis"` as fallback.
#' @param xlab Optional character string to override the x-axis label. Defaults to the name of `cov1`.
#'
#' @return A `ggplot2` plot object showing predicted occupancy across a grid of `cov1` and `cov2` values.
#'
#' @import ggplot2
#' @importFrom unmarked predict
#' @importFrom grDevices hcl.colors
#' @importFrom progress progress_bar
#'
#' @examples
#' if (requireNamespace("unmarked", quietly = TRUE)) {
#'   # Simulate data
#'   set.seed(123)
#'   n_sites <- 50
#'   n_visits <- 3
#'   site_covs <- data.frame(
#'     elev = scale(rnorm(n_sites)),
#'     forest = scale(rnorm(n_sites))
#'   )
#'   obs_covs <- list(effort_s = matrix(rnorm(n_sites * n_visits), n_sites, n_visits))
#'   y <- matrix(rbinom(n_sites * n_visits, 1, 0.5), n_sites, n_visits)
#'
#'   umf <- unmarked::unmarkedFrameOccu(y = y, siteCovs = site_covs, obsCovs = obs_covs)
#'   mod <- unmarked::occu(~ effort_s ~ elev * forest, data = umf)
#'
#'   # Plot interaction surface
#'   plot <- plot_interaction_surface1(
#'     model = mod,
#'     cov1 = "elev",
#'     cov2 = "forest",
#'     grid_length = 100,
#'     chunk_size = 5000,
#'     xlab = "elev"
#'   )
#'   print(plot)
#' }
#'
#' @export

plot_occupancy_interaction <- function(
    model,
    cov1,
    cov2,
    grid_length = 400,
    chunk_size = 50000,
    palette = NULL,  # palette is optional; uses "cividis" by default
    xlab = NULL
) {

  # Extract site-level covariates
  site_covs <- model@data@siteCovs

  # Validate covariate names
  if (!(cov1 %in% names(site_covs)) | !(cov2 %in% names(site_covs))) {
    stop("Both covariates must be present in model.")
  }

  # Get covariate ranges
  cov1_range <- range(site_covs[[cov1]], na.rm = TRUE)
  cov2_range <- range(site_covs[[cov2]], na.rm = TRUE)

  # Create prediction grid
  grid <- expand.grid(
    cov1_vals = seq(cov1_range[1], cov1_range[2], length.out = grid_length),
    cov2_vals = seq(cov2_range[1], cov2_range[2], length.out = grid_length)
  )
  names(grid) <- c(cov1, cov2)

  # Initialise prediction column
  grid$psi <- NA_real_
  n_rows <- nrow(grid)

  # Chunking
  chunks <- split(grid, ceiling(seq_len(n_rows) / chunk_size))
  pb <- progress::progress_bar$new(
    format = "Predicting [:bar] :percent ETA: :eta",
    total = length(chunks), clear = FALSE, width = 60
  )

  # All site-level covariate names from the original model
  all_covs <- names(site_covs)

  for (i in seq_along(chunks)) {
    chunk_data <- chunks[[i]]

    # Find any covariates that are missing from the chunk of the prediction grid
    # Fill in missing covariates with numeric NA
    missing_covs <- setdiff(all_covs, names(chunk_data))
    for (var in missing_covs) {
      chunk_data[[var]] <- NA_real_   # a missing numeric NA
    }

    # Match column order to model
    chunk_data <- chunk_data[all_covs]

    # Predict state (occupancy)
    preds <- unmarked::predict(model, type = "state", newdata = chunk_data, appendData = FALSE)$Predicted
    grid$psi[as.numeric(rownames(chunk_data))] <- preds

    pb$tick()
  }

  # Set up color palette
  col_palette <- if (is.null(palette)) {
    hcl.colors(256, "cividis")
  } else {
    hcl.colors(256, palette)
  }

  # Use provided or default x-axis label
  x_label <- if (is.null(xlab)) cov1 else xlab

  # Plot
  p <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = .data[[cov1]],
      y = psi,
      group = .data[[cov2]],
      colour = .data[[cov2]]
    )
  ) +
    ggplot2::geom_line(linewidth = 0.5, alpha = 0.8) +
    ggplot2::scale_colour_gradientn(
      colours = col_palette,
      name = cov2
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Probability of occupancy"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 15),
      axis.title = ggplot2::element_text(size = 17, vjust = 0.8),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 15),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    )

  return(p)
}
