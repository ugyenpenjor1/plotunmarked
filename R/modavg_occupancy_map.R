
#' Create and plot model-averaged occupancy prediction map
#'
#' Generates a spatial map of model-averaged occupancy predictions using a list of fitted models,
#' prediction covariate raster layers, and optionally plots confidence intervals. The prediction raster
#' can be aggregated to reduce resolution, and the resulting map can be saved to a file.
#'
#' @param d2_mod_list A list of candidate occupancy models (usually output from \code{unmarked} or
#'   compatible modeling functions) to be used for model-averaged prediction.
#' @param pred_cov_stack A \code{SpatRaster} (from the \code{terra} package) containing covariate layers
#'   used for prediction. Must have layer names matching the covariates in the models.
#' @param agg_factor Integer. Factor by which to aggregate (downsample) the raster resolution to speed up
#'   processing. Default is 4.
#' @param plot_ci Logical. If \code{TRUE}, 95% confidence interval maps (lower and upper bounds) will be plotted
#'   alongside predicted occupancy and standard error. Default is \code{FALSE}.
#' @param save_plot Logical. If \code{TRUE}, the generated plot will be saved to disk. Default is \code{FALSE}.
#' @param plot_filename Character. Filename (without extension) to save the plot. Default is \code{"occupancy_plot"}.
#' @param plot_format Character. File format to save the plot (e.g., \code{"jpeg"}, \code{"png"}, \code{"pdf"}).
#'   Default is \code{"jpeg"}.
#'
#' @return A list (invisibly) containing:
#' \item{pred_df}{A data.frame with predicted occupancy, standard errors, confidence intervals, and coordinates.}
#' \item{pred_raster}{A \code{SpatRaster} of predicted occupancy values.}
#' \item{mod_avg_pred}{A list of numeric vectors containing model-averaged predictions and uncertainty measures.}
#'
#' @importFrom terra aggregate rast
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_viridis_c coord_equal theme_minimal labs facet_wrap theme element_rect element_text element_blank ggsave
#' @importFrom progress progress_bar
#' @importFrom AICcmodavg modavgPred
#'
#' @examples
#' \dontrun{
#' # Assuming d2_mod_list is your list of models and pred_cov_stack is your raster stack:
#' result <- modavg_occupancy_map(d2_mod_list, pred_cov_stack, agg_factor = 2, plot_ci = TRUE)
#' }
#'
#' @export

modavg_occupancy_map <- function(
    d2_mod_list,
    pred_cov_stack,
    agg_factor = 4,
    plot_ci = FALSE,
    save_plot = FALSE,
    plot_filename = "occupancy_plot",
    plot_format = "jpeg"
) {
  if (is.null(names(pred_cov_stack))) stop("pred_cov_stack must have layer names matching model covariates.")

  pred_cov_reduced <- terra::aggregate(pred_cov_stack, fact = agg_factor, fun = mean)
  newdata <- as.data.frame(pred_cov_reduced, xy = TRUE)

  n <- nrow(newdata)
  chunk_size <- if (n <= 1000) n else if (n <= 5000) 500 else 1000

  pb <- progress::progress_bar$new(
    format = "Preparing your map [:bar] :percent ETA: :eta",
    total = ceiling(n / chunk_size), clear = FALSE, width = 60
  )

  mod_avg_preds <- list(
    mod.avg.pred = numeric(n),
    uncond.se = numeric(n),
    lower.CL = numeric(n),
    upper.CL = numeric(n)
  )

  idx_start <- seq(1, n, by = chunk_size)
  idx_end <- pmin(idx_start + chunk_size - 1, n)

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
  pred_df$SE <- pmax(pred_df$SE, 0)  # SE can't be negative but can exceed 1 if wildly wrong
  pred_df$lower <- pmin(pmax(pred_df$lower, 0), 1)
  pred_df$upper <- pmin(pmax(pred_df$upper, 0), 1)

  # ---- Convert to long format for ggplot2 ----
  plot_df <- rbind(
    data.frame(x = pred_df$x, y = pred_df$y, value = pred_df$Predicted, type = "Predicted occupancy (mean)"),
    data.frame(x = pred_df$x, y = pred_df$y, value = pred_df$SE, type = "Standard error")
  )

  if (plot_ci) {
    plot_df <- rbind(
      plot_df,
      data.frame(x = pred_df$x, y = pred_df$y, value = pred_df$lower, type = "95% CI - lower"),
      data.frame(x = pred_df$x, y = pred_df$y, value = pred_df$upper, type = "95% CI - upper")
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

  # ---- ggplot2 raster map ----
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(
      name = "Occupancy\nprobability",
      option = "magma",
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

  return(invisible(list(
    pred_df = pred_df,
    pred_raster = pred_raster,
    mod_avg_pred = mod_avg_preds
  )))
}
