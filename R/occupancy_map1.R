#' Generate occupancy prediction map from a single unmarked model
#'
#' This function predicts occupancy probabilities from a single fitted \code{unmarked}
#' occupancy model over a spatial covariate stack.
#' It aggregates input covariates, computes predictions in chunks for efficiency,
#' and plots the results with options for confidence intervals and saving the plot.
#'
#' The function outputs a list containing prediction data, a raster of predicted occupancy,
#' and the raw predicted values with standard errors and confidence limits.
#'
#' @param fitted_model An \code{unmarkedFit} model object fitted with \code{unmarked}.
#' @param pred_cov_stack A \code{terra::SpatRaster} or \code{raster::RasterLayer} containing covariate layers named to match those used in the model. Raster inputs from the \code{raster} package are automatically converted.
#' @param agg_factor Integer factor to aggregate (resample) the covariate raster to coarser resolution for faster prediction. Default is 4.
#' @param plot_ci Logical indicating whether to include 95% confidence intervals in the plot. Default is \code{FALSE}.
#' @param save_plot Logical indicating whether to save the occupancy plot to disk. Default is \code{FALSE}.
#' @param plot_filename Character string filename (without extension) for saving the plot. Default is \code{"occupancy_map"}.
#' @param plot_format Character string for the file format (e.g., \code{"png"}, \code{"pdf"}). Default is \code{"png"}.
#' @param mask_raster Optional \code{terra::SpatRaster} used to mask or clip the output raster. Default is \code{NULL}.
#'
#' @return A (invisible) list with components:
#' \item{pred_df}{A data frame containing predicted occupancy, standard error, confidence intervals, and covariate values.}
#' \item{pred_raster}{A \code{terra::SpatRaster} raster of predicted occupancy probabilities.}
#' \item{pred_values}{A data frame of raw prediction results including standard errors and confidence limits.}
#'
#' @details
#' The function automatically constructs polynomial and interaction terms required by the model formula
#' and predicts occupancy probabilities chunk-wise to improve performance on large datasets.
#'
#' @examples
#' \dontrun{
#' library(unmarked)
#' library(terra)
#'
#' # Assume 'umf' is an unmarkedFrame and 'fitted_model' is an unmarkedFit occupancy model
#' # Assume 'covariate_stack' is a SpatRaster (also accepts RasterLayer object) with layers named as in the model
#'
#' result <- occupancy_map1(
#'   fitted_model = fitted_model,
#'   pred_cov_stack = covariate_stack,
#'   agg_factor = 4,
#'   plot_ci = TRUE,
#'   save_plot = TRUE,
#'   plot_filename = "occupancy_map",
#'   plot_format = "png"
#' )
#'
#' # View predictions dataframe
#' head(result$pred_df)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn coord_equal theme_minimal labs facet_wrap theme element_rect element_text element_blank ggsave
#' @importFrom terra rast aggregate mask crs
#' @importFrom progress progress_bar
#' @export

occupancy_map1 <- function(
    fitted_model,
    pred_cov_stack,
    agg_factor = 4,
    plot_ci = FALSE,
    save_plot = FALSE,
    plot_filename = "occupancy_map",
    plot_format = "png",
    mask_raster = NULL
) {
  # Check model class
  if (!inherits(fitted_model, "unmarkedFit")) {
    stop("fitted_model must be an 'unmarkedFit' object.")
  }

  # Convert to terra::SpatRaster if needed
  if (inherits(pred_cov_stack, "Raster")) {
    pred_cov_stack <- terra::rast(pred_cov_stack)
  } else if (!inherits(pred_cov_stack, "SpatRaster")) {
    stop("pred_cov_stack must be a 'SpatRaster' or 'RasterStack/Brick'.")
  }

  if (is.null(names(pred_cov_stack))) {
    stop("pred_cov_stack must have layer names matching model covariates.")
  }

  pred_cov_reduced <- terra::aggregate(pred_cov_stack, fact = agg_factor, fun = mean)
  newdata <- as.data.frame(pred_cov_reduced, xy = TRUE)

  # Build interaction/polynomial terms
  build_model_matrix_terms <- function(formula_obj, base_data) {
    rhs <- as.formula(paste("~", as.character(formula(formula_obj))[3]))
    mm <- model.matrix(rhs, base_data)
    needed_terms <- setdiff(colnames(mm), "(Intercept)")

    for (term in needed_terms) {
      if (!(term %in% names(base_data))) {
        try({
          base_data[[term]] <- eval(parse(text = term), envir = base_data)
        }, silent = TRUE)
      }
    }

    # Warn if still missing
    missing_terms <- needed_terms[!needed_terms %in% names(base_data)]
    if (length(missing_terms) > 0) {
      warning("The following terms could not be created: ", paste(missing_terms, collapse = ", "))
    }

    return(base_data)
  }

  newdata <- build_model_matrix_terms(fitted_model, newdata)

  # Predict in chunks
  n <- nrow(newdata)
  chunk_size <- if (n <= 1000) n else if (n <= 5000) 500 else 1000

  pb <- progress::progress_bar$new(
    format = "Preparing your map [:bar] :percent ETA: :eta",
    total = ceiling(n / chunk_size), clear = FALSE, width = 60
  )

  pred_values <- data.frame(
    Predicted = numeric(n),
    SE = numeric(n),
    lower = numeric(n),
    upper = numeric(n)
  )

  idx_start <- seq(1, n, by = chunk_size)
  idx_end <- pmin(idx_start + chunk_size - 1, n)

  for (i in seq_along(idx_start)) {
    rows <- idx_start[i]:idx_end[i]
    chunk_newdata <- newdata[rows, , drop = FALSE]

    chunk_pred <- tryCatch(
      unmarked::predict(
        fitted_model,
        type = "state",
        newdata = chunk_newdata,
        appendData = FALSE
      ),
      error = function(e) {
        warning("Prediction failed for chunk ", i, ": ", e$message)
        data.frame(Predicted = rep(NA, nrow(chunk_newdata)),
                   SE = NA, lower = NA, upper = NA)
      }
    )

    pred_values[rows, ] <- chunk_pred
    pb$tick()
  }

  # Clamp values
  pred_values$Predicted <- pmin(pmax(pred_values$Predicted, 0), 1)
  pred_values$SE <- pmax(pred_values$SE, 0)
  pred_values$lower <- pmin(pmax(pred_values$lower, 0), 1)
  pred_values$upper <- pmin(pmax(pred_values$upper, 0), 1)

  pred_df <- cbind(pred_values, newdata)

  # Prepare for plotting
  res_x <- min(diff(sort(unique(pred_df$x))))
  res_y <- min(diff(sort(unique(pred_df$y))))

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

  plot_df$type <- factor(
    plot_df$type,
    levels = c("Predicted occupancy (mean)", "Standard error", "95% CI - lower", "95% CI - upper")
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile(width = res_x, height = res_y) +
    ggplot2::scale_fill_gradientn(
      name = "Occupancy\nprobability",
      colours = hcl.colors(256, "Viridis"),
      na.value = "transparent"
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
  pred_raster <- terra::rast(xyz_mat, type = "xyz", crs = terra::crs(pred_cov_stack))

  # Apply optional mask
  if (!is.null(mask_raster)) {
    pred_raster <- terra::mask(pred_raster, mask_raster)
  }

  return(invisible(list(
    pred_df = pred_df,
    pred_raster = pred_raster,
    pred_values = pred_values
  )))
}
