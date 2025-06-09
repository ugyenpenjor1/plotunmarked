
#' Model-averaged occupancy prediction and mapping
#'
#' Generates model-averaged occupancy predictions across a spatial landscape using
#' a list of fitted occupancy models (`unmarked` objects) and a raster stack of covariates.
#' The function supports polynomial and interaction terms used in occupancy formulas,
#' performs raster aggregation for computational efficiency, and optionally plots
#' the predictions with uncertainty estimates.
#'
#' @param d2_mod_list A list of fitted occupancy models (from `unmarked::occu()`),
#'   created with the same dataset and covariate names as used in `pred_cov_stack`.
#' @param pred_cov_stack A `terra::SpatRaster` or `raster::RasterStack`/`RasterBrick`
#'   of scaled covariates used for prediction. Layer names must match model covariates.
#' @param agg_factor Integer. Factor by which to spatially aggregate the prediction raster
#'   for efficiency (default is 4). Use 1 to skip aggregation.
#' @param plot_ci Logical. If `TRUE`, plots 95% confidence intervals for predictions.
#'   Default is `FALSE`.
#' @param save_plot Logical. If `TRUE`, saves the generated plot to disk. Default is `FALSE`.
#' @param plot_filename Character. File name (without extension) to save the plot if `save_plot = TRUE`.
#'   Default is `"occupancy_map"`.
#' @param plot_format Character. Format of the saved plot (e.g., `"png"`, `"pdf"`). Default is `"png"`.
#'
#' @return An invisible list containing:
#' \describe{
#'   \item{`pred_df`}{A `data.frame` with x, y coordinates and predicted occupancy values with uncertainty.}
#'   \item{`pred_raster`}{A `terra::SpatRaster` with the predicted occupancy surface.}
#'   \item{`mod_avg_pred`}{A list with numeric vectors of model-averaged predictions and standard errors.}
#' }
#'
#' @details
#' The function parses each model formula to identify any polynomial (e.g., `I(elev^2)`)
#' or interaction (e.g., `elev:forest`) terms and automatically generates those terms
#' in the prediction dataset. This ensures that `modavgPred()` from `AICcmodavg` can
#' correctly evaluate predictions even when complex terms are used.
#'
#' @note Covariates in `pred_cov_stack` must already be scaled consistently with those
#' used in model fitting.
#'
#' @importFrom terra rast aggregate
#' @importFrom AICcmodavg modavgPred
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_viridis_c coord_equal theme_minimal facet_wrap
#' @importFrom progress progress_bar
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume d2_mod_list is a list of fitted occu models
#' # and pred_cov_stack is a SpatRaster of scaled covariates
#' result <- modavg_occupancy_map2(
#'   d2_mod_list = d2_mod_list,
#'   pred_cov_stack = terra_list1,
#'   agg_factor = 4,
#'   plot_ci = TRUE,
#'   save_plot = TRUE,
#'   plot_filename = "map_output",
#'   plot_format = "png"
#' )
#' }

modavg_occupancy_map2 <- function(
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

  # Helper function to add polynomial and interaction terms
  build_model_matrix_terms <- function(formula_list, base_data) {
    all_terms <- list()

    for (mod in formula_list) {
      form <- formula(mod@formula)
      psi_rhs <- as.formula(paste("~", as.character(form)[3]))
      mm <- model.matrix(psi_rhs, base_data)
      all_terms[[length(all_terms) + 1]] <- colnames(mm)
    }

    final_terms <- unique(unlist(all_terms))
    final_terms <- setdiff(final_terms, "(Intercept)")

    for (term in final_terms) {
      if (!(term %in% names(base_data))) {
        # Safely evaluate and create the missing term
        try({
          base_data[[term]] <- eval(parse(text = term), envir = base_data)
        }, silent = TRUE)
      }
    }
    return(base_data)
  }

  # Generate any polynomial or interaction terms required
  newdata <- build_model_matrix_terms(d2_mod_list, newdata)

  # Chunking for computation efficiency
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

  pred_df$Predicted <- pmin(pmax(pred_df$Predicted, 0), 1)
  pred_df$SE <- pmax(pred_df$SE, 0)
  pred_df$lower <- pmin(pmax(pred_df$lower, 0), 1)
  pred_df$upper <- pmin(pmax(pred_df$upper, 0), 1)

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

  plot_df$type <- factor(
    plot_df$type,
    levels = c(
      "Predicted occupancy (mean)",
      "Standard error",
      "95% CI - lower",
      "95% CI - upper"
    )
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    # ggplot2::scale_fill_viridis_c(
    #   name = "Occupancy\nprobability",
    #   option = "inferno",
    #   na.value = "transparent",
    #   limits = c(0, 1)
    # ) +
    ggplot2::scale_fill_gradientn(
      name = "Occupancy\nprobability",
      colours = hcl.colors(256, "Batlow")#,
      #limits = c(0, 1)
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

  # Create raster output
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
