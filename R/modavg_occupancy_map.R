#' Model-averaged occupancy map prediction and visualisation
#'
#' Generates spatial predictions of occupancy probabilities using model averaging from a set of fitted `unmarked` models. The function aggregates input covariate rasters for efficiency, computes model-averaged predictions, and optionally plots and saves the resulting map with associated uncertainty.
#'
#' @param d2_mod_list A candidate model set formatted as a list of fitted `unmarkedFit` objects, as required by `AICcmodavg::modavgPred()`. Must include model selection information.
#' @param pred_cov_stack A `SpatRaster` (from the `terra` package) or a `RasterStack`/`RasterBrick` (from the `raster` package) containing named covariate layers used in the candidate models.
#' @param agg_factor Integer. Aggregation factor for reducing raster resolution (default is `4`). Helps reduce computational cost for large rasters.
#' @param plot_ci Logical. If `TRUE`, includes 95% confidence interval bounds (`lower`, `upper`) in the plot. Default is `FALSE`.
#' @param save_plot Logical. If `TRUE`, saves the plot to a file using the specified filename and format. Default is `FALSE`.
#' @param plot_filename Character. Base name of the output file (excluding extension). Default is `"occupancy_map"`.
#' @param plot_format Character. File format for saving the plot (e.g., `"png"`, `"jpeg"`). Default is `"png"`.
#'
#' @return A (named) list containing:
#' \describe{
#'   \item{`pred_df`}{A data frame with x/y coordinates, predicted occupancy (`Predicted`), standard error (`SE`), and confidence interval bounds (`lower`, `upper`).}
#'   \item{`pred_raster`}{A `SpatRaster` of predicted occupancy probabilities.}
#'   \item{`mod_avg_pred`}{A list of numeric vectors: `mod.avg.pred`, `uncond.se`, `lower.CL`, `upper.CL` for each raster cell.}
#' }
#' The raster resolution of the output matches the aggregated resolution specified by `agg_factor`.
#'
#' @details This function performs model averaging over multiple candidate models using the `AICcmodavg` package and spatial covariate data. Predictions are performed on an aggregated version of the input raster stack to improve speed. A chunking mechanism is used to process raster cells in batches and a progress bar is shown during computation.
#'
#' If `plot_ci = TRUE`, the output map will include lower and upper 95% confidence bounds. All maps are created using `ggplot2`.
#'
#' @importFrom terra rast aggregate
#' @importFrom raster stack brick
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_viridis_c coord_equal theme_minimal labs facet_wrap theme element_rect element_text ggsave
#' @importFrom progress progress_bar
#' @importFrom AICcmodavg modavgPred
#' @export
#'
#' @examples
#' \dontrun{
#' library(unmarked)
#' library(terra)
#' library(AICcmodavg)
#'
#' # Create synthetic data and fit two simple models
#' umf <- unmarkedFrameOccu(y = matrix(sample(0:1, 100, TRUE), ncol = 5),
#'                          siteCovs = data.frame(elev = rnorm(20), forest = rnorm(20)))
#' m1 <- occu(~1 ~ elev, data = umf)
#' m2 <- occu(~1 ~ forest, data = umf)
#' mod_list <- list(m1 = m1, m2 = m2)
#'
#' # Prepare raster covariates
#' elev_rast <- rast(matrix(rnorm(100), 10, 10))
#' forest_rast <- rast(matrix(rnorm(100), 10, 10))
#' names(elev_rast) <- "elev"
#' names(forest_rast) <- "forest"
#' cov_stack <- c(elev_rast, forest_rast)
#'
#' # Using raster package
#' elev_rast <- raster(matrix(rnorm(100), 10, 10))
#' forest_rast <- raster(matrix(rnorm(100), 10, 10))
#' names(elev_rast) <- "elev"
#' names(forest_rast) <- "forest"
#' cov_stack <- stack(elev_rast, forest_rast)
#'
#' # Generate occupancy map using model averaging
#' result <- modavg_occupancy_map(d2_mod_list = mod_list, pred_cov_stack = cov_stack)
#' }


modavg_occupancy_map <- function(
    d2_mod_list,
    pred_cov_stack,
    agg_factor = 4,
    plot_ci = FALSE,
    save_plot = FALSE,
    plot_filename = "occupancy_map",
    plot_format = "png"
) {
  # Convert raster::Raster* objects to terra::SpatRaster
  if (inherits(pred_cov_stack, "Raster")) {
    pred_cov_stack <- terra::rast(pred_cov_stack)
  } else if (!inherits(pred_cov_stack, "SpatRaster")) {
    stop("pred_cov_stack must be a 'SpatRaster' (terra) or 'RasterStack/Brick' (raster).")
  }

  if (is.null(names(pred_cov_stack))) {
    stop("pred_cov_stack must have layer names matching model covariates.")
  }

  pred_cov_reduced <- terra::aggregate(
    pred_cov_stack,
    fact = agg_factor,
    fun = mean
  )
  newdata <- as.data.frame(pred_cov_reduced, xy = TRUE)

  # Chunking for computation efficiency
  n <- nrow(newdata)
  chunk_size <- if (n <= 1000) n else if (n <= 5000) 500 else 1000

  # Progress bar
  pb <- progress::progress_bar$new(
    format = "Preparing your map [:bar] :percent ETA: :eta",
    total = ceiling(n / chunk_size), clear = FALSE, width = 60
  )

  # List to hold predicted values
  mod_avg_preds <- list(
    mod.avg.pred = numeric(n),
    uncond.se = numeric(n),
    lower.CL = numeric(n),
    upper.CL = numeric(n)
  )

  idx_start <- seq(1, n, by = chunk_size)
  idx_end <- pmin(idx_start + chunk_size - 1, n)

  # Do model averaging
  for (i in seq_along(idx_start)) {
    rows <- idx_start[i]:idx_end[i]
    chunk_newdata <- newdata[rows, , drop = FALSE]

    chunk_pred <- AICcmodavg::modavgPred(
      cand.set = d2_mod_list,
      parm.type = "psi",
      newdata = chunk_newdata
    )

    mod_avg_preds$mod.avg.pred[rows] <- chunk_pred$mod.avg.pred
    mod_avg_preds$uncond.se[rows] <- chunk_pred$uncond.se
    mod_avg_preds$lower.CL[rows] <- chunk_pred$lower.CL
    mod_avg_preds$upper.CL[rows] <- chunk_pred$upper.CL

    pb$tick()
  }

  pred_df <- data.frame(
    Predicted = mod_avg_preds$mod.avg.pred,
    SE = mod_avg_preds$uncond.se,
    lower = mod_avg_preds$lower.CL,
    upper = mod_avg_preds$upper.CL,
    newdata
  )

  # Clamp values to valid probability range [0, 1]
  pred_df$Predicted <- pmin(pmax(pred_df$Predicted, 0), 1)
  pred_df$SE <- pmax(pred_df$SE, 0)
  pred_df$lower <- pmin(pmax(pred_df$lower, 0), 1)
  pred_df$upper <- pmin(pmax(pred_df$upper, 0), 1)

  # Convert to long format for ggplot2
  plot_df <- rbind(
    data.frame(
      x = pred_df$x,
      y = pred_df$y,
      value = pred_df$Predicted,
      type = "Predicted occupancy (mean)"
    ),
    data.frame(
      x = pred_df$x,
      y = pred_df$y,
      value = pred_df$SE,
      type = "Standard error"
    )
  )

  if (plot_ci) {
    plot_df <- rbind(
      plot_df,
      data.frame(
        x = pred_df$x,
        y = pred_df$y,
        value = pred_df$lower,
        type = "95% CI - lower"
      ),
      data.frame(
        x = pred_df$x,
        y = pred_df$y,
        value = pred_df$upper,
        type = "95% CI - upper"
      )
    )
  }

  # Reorder facet panels
  plot_df$type <- factor(
    plot_df$type,
    levels = c(
      "Predicted occupancy (mean)",
      "Standard error",
      "95% CI - lower",
      "95% CI - upper"
    )
  )

  # Raster map
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(
      name = "Occupancy\nprobability",
      option = "viridis",
      na.value = "transparent",
      limits = c(0, 1)
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::labs(fill = NULL, x = NULL, y = NULL) +
    ggplot2::facet_wrap(~ type, ncol = 2) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"),
      axis.text.y = ggplot2::element_text(angle = 90, vjust = 0, hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )

  print(p)

  # Save if requested
  if (save_plot) {
    ggplot2::ggsave(
      filename = paste0(plot_filename, ".", tolower(plot_format)),
      plot = p,
      width = 10,
      height = if (plot_ci) 8 else 5,
      dpi = 300,
      bg = "white",
      device = plot_format
    )
  }

  # Create predicted raster for output
  xyz_mat <- cbind(x = pred_df$x, y = pred_df$y, z = pred_df$Predicted)
  pred_raster <- terra::rast(xyz_mat, type = "xyz")

  return(
    invisible(
      list(
        pred_df = pred_df,
        pred_raster = pred_raster,
        mod_avg_pred = mod_avg_preds
      )
    )
  )
}
