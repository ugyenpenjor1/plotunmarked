#' Predict and visualise occupancy from an unmarked model
#'
#' Computes spatial predictions of occupancy probabilities from a fitted `unmarkedFit` model and raster covariates. Optionally plots the map and saves it to a file.
#'
#' @param model A fitted model of class `unmarkedFit`, typically from the `unmarked` package.
#' @param covariate_rasters A named list of raster covariates (`SpatRaster` or `RasterLayer`), where names match the occupancy covariates used in the model. The list must include all covariates used in the model, excluding the intercept.
#' @param return_class Character. Specifies the class of raster returned. Options are `"same"` (default, returns the same class as input), `"terra"` (returns `SpatRaster`), or `"raster"` (returns `RasterLayer`).
#' @param plot_map Logical. If `TRUE` (default), plots the predicted occupancy map and uncertainty layers.
#' @param plot_ci Logical. If `TRUE`, plots confidence interval bounds on the occupancy scale.
#' @param ci_level Numeric. Confidence level for intervals (default is `0.95` for 95% CI).
#' @param save_plot Logical. If `TRUE`, saves the plot to a file using the specified format and filename.
#' @param plot_filename Character. Base filename to save the plot (without file extension). Default is `"occupancy_map"`.
#' @param plot_format Character. File format to use when saving the plot (e.g., `"jpeg"`, `"png"`). Default is `"png"`.
#'
#' @return A (named) list of raster objects containing:
#' \describe{
#'   \item{`mean`}{Raster of predicted occupancy probabilities.}
#'   \item{`se`}{Raster of standard errors on the probability scale.}
#'   \item{`lower`, `upper`}{(If `plot_ci = TRUE`) Rasters of lower and upper bounds of the confidence interval.}
#' }
#' All returned raster layers match the format specified by `return_class`.
#'
#' @details The function uses model coefficients and their variance-covariance matrix to generate occupancy predictions across a spatial surface defined by the input covariate rasters. The intercept term is handled internally by generating a constant raster. Variance propagation is performed via the delta method.
#'
#' If `plot_map = TRUE`, the function creates a `ggplot2`-based map visualising the predicted occupancy, standard error, and optionally, the confidence interval bounds. Plots are saved using `ggsave()` if `save_plot = TRUE`.
#'
#' A progress bar is shown during uncertainty calculation if the `pbapply` package is installed.
#'
#' @importFrom terra rast app values global
#' @importFrom raster raster cellStats
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_viridis_c coord_equal theme_minimal labs facet_wrap theme element_rect element_text ggsave
#' @importFrom stats qnorm
#' @importFrom unmarked coef vcov
#' @importFrom pbapply pbapply
#' @export
#'
#' @examples
#' \dontrun{
#' library(unmarked)
#' library(terra)
#'
#' # Fit a simple model (example assumes you have data in unmarkedFrame)
#' umf <- unmarkedFrameOccu(y = matrix(sample(0:1, 100, TRUE), ncol = 5),
#'                          siteCovs = data.frame(elev = rnorm(20), forest = rnorm(20)))
#' mod <- occu(~1 ~ elev + forest, data = umf)
#'
#' # Prepare raster covariates
#' elev_rast <- rast(matrix(rnorm(100), 10, 10))
#' forest_rast <- rast(matrix(rnorm(100), 10, 10))
#' names(elev_rast) <- "elev"
#' names(forest_rast) <- "forest"
#'
#' covariate_rasters <- list(elev = elev_rast, forest = forest_rast)
#'
#' # Generate and plot map
#' result <- occupancy_map(model = mod, covariate_rasters = covariate_rasters)
#' }


occupancy_map <- function(
    model,
    covariate_rasters,
    return_class = c("same", "terra", "raster"),
    plot_map = TRUE,
    plot_ci = FALSE,
    ci_level = 0.95,
    save_plot = FALSE,
    plot_filename = "occupancy_map",
    plot_format = "png"
) {
  return_class <- match.arg(return_class)

  if (!inherits(model, "unmarkedFit")) {
    stop("`model` must be an 'unmarkedFit' object.")
  }

  if (!is.list(covariate_rasters) || is.null(names(covariate_rasters))) {
    stop("`covariate_rasters` must be a named list of rasters with names matching model covariates.")
  }

  # # Extract coefficients and VCV matrix
  beta <- unmarked::coef(model)
  vcov_mat <- unmarked::vcov(model)

  psi_coefs <- beta[grep("^psi\\(", names(beta))]
  model_covariate_names <- gsub("^psi\\((.*)\\)$", "\\1", names(psi_coefs))
  names(psi_coefs) <- model_covariate_names

  # Check covariate raster names match
  # Ignore intercept ("Int") when checking for covariate rasters
  missing_covs <- setdiff(setdiff(model_covariate_names, "Int"), names(covariate_rasters))
  if (length(missing_covs) > 0) {
    stop("Missing raster(s) for the following model covariates: ", paste(missing_covs, collapse = ", "))
  } else {
    message("Raster list names match model covariates: ", paste(model_covariate_names, collapse = ", "))
  }

  coef_names <- names(psi_coefs)
  beta_vec <- psi_coefs

  # Extract relevant VCV matrix
  beta_idx <- grep("^psi\\(", rownames(vcov_mat))
  V_sub <- vcov_mat[beta_idx, beta_idx]
  rownames(V_sub) <- colnames(V_sub) <- model_covariate_names

  # Create raster stack of covariates (suppor both terra and raster)
  X_stack <- terra::rast(lapply(coef_names, function(cov) {
    if (cov == "Int") {
      # Intercept: create constant raster
      base_raster <- covariate_rasters[[1]]
      if (inherits(base_raster, "Raster")) base_raster <- terra::rast(base_raster)
      return(base_raster * 0 + 1)
    } else {
      rast <- covariate_rasters[[cov]]
      if (!inherits(rast, c("SpatRaster", "RasterLayer"))) {
        stop(sprintf("Raster for '%s' must be SpatRaster or RasterLayer.", cov))
      }
      if (inherits(rast, "Raster")) {
        rast <- terra::rast(rast)
      }
      return(rast)
    }
  }))
  names(X_stack) <- coef_names

  # Check raster consistency
  for (cov in setdiff(coef_names, "Int")) {
    if (!(cov %in% names(covariate_rasters))) {
      stop("Missing raster for covariate: ", cov)
    }
    r <- covariate_rasters[[cov]]
    if (inherits(r, "SpatRaster")) {
      if (terra::global(r, fun = function(x) sum(!is.na(x)))[1, 1] == 0) {
        stop(sprintf("Raster for '%s' has no non-NA values.", cov))
      }
    } else {
      if (raster::cellStats(r, stat = function(x, ...) sum(!is.na(x))) == 0) {
        stop(sprintf("Raster for '%s' has no non-NA values.", cov))
      }
    }
  }

  # Compute logit prediction
  logitPsi <- terra::app(X_stack, fun = function(x_row) {
    sum(x_row * beta_vec)
  })

  message("Preparing your map...")

  if (requireNamespace("pbapply", quietly = TRUE)) {
    vals <- terra::values(X_stack)
    pbapply <- getNamespace("pbapply")
    var_vals <- pbapply::pbapply(vals, 1, function(x) {
      x_mat <- matrix(x, ncol = 1)
      as.numeric(t(x_mat) %*% V_sub %*% x_mat)
    })
  } else {
    warning("Install 'pbapply' package to see progress bar.")
    vals <- terra::values(X_stack)
    var_vals <- apply(vals, 1, function(x) {
      x_mat <- matrix(x, ncol = 1)
      as.numeric(t(x_mat) %*% V_sub %*% x_mat)
    })
  }

  # Reconstruct variance raster
  var_logitPsi <- logitPsi
  terra::values(var_logitPsi) <- var_vals

  # Compute mean and SE on probability scale
  psi_mean <- 1 / (1 + exp(-logitPsi))
  psi_se <- (exp(-logitPsi) / (1 + exp(-logitPsi))^2) * sqrt(var_logitPsi)

  # CI bounds
  if (plot_ci) {
    z_val <- qnorm(1 - (1 - ci_level) / 2)
    lower_logit <- logitPsi - z_val * sqrt(var_logitPsi)
    upper_logit <- logitPsi + z_val * sqrt(var_logitPsi)

    psi_lower <- 1 / (1 + exp(-lower_logit))
    psi_upper <- 1 / (1 + exp(-upper_logit))
  }

  # Convert raster class if needed
  convert_rast <- function(r) {
    if (return_class == "raster" && inherits(r, "SpatRaster")) {
      raster::raster(r)
    } else if (return_class == "terra" && inherits(r, "RasterLayer")) {
      terra::rast(r)
    } else {
      r
    }
  }

  psi_mean <- convert_rast(psi_mean)
  psi_se <- convert_rast(psi_se)
  if (plot_ci) {
    psi_lower <- convert_rast(psi_lower)
    psi_upper <- convert_rast(psi_upper)
  }

  # Plotting
  if (plot_map) {
    # Convert to terra::SpatRaster for plotting only
    plot_mean <- if (inherits(psi_mean, "SpatRaster")) psi_mean else terra::rast(psi_mean)
    plot_se <- if (inherits(psi_se, "SpatRaster")) psi_se else terra::rast(psi_se)

    # Convert rasters to data frames
    df_mean <- as.data.frame(psi_mean, xy = TRUE, na.rm = TRUE)
    names(df_mean)[3] <- "value"
    df_mean$type <- "Predicted occupancy (mean)"

    df_se <- as.data.frame(psi_se, xy = TRUE, na.rm = TRUE)
    names(df_se)[3] <- "value"
    df_se$type <- "Standard error"

    plot_df <- rbind(df_mean, df_se)

    if (plot_ci) {
      # Convert to terra::SpatRaster for plotting only
      plot_lower <- if (inherits(psi_lower, "SpatRaster")) psi_lower else terra::rast(psi_lower)
      plot_upper <- if (inherits(psi_upper, "SpatRaster")) psi_upper else terra::rast(psi_upper)

      df_lower <- as.data.frame(psi_lower, xy = TRUE, na.rm = TRUE)
      names(df_lower)[3] <- "value"
      df_lower$type <- paste0(ci_level * 100, "% CI - lower")

      df_upper <- as.data.frame(psi_upper, xy = TRUE, na.rm = TRUE)
      names(df_upper)[3] <- "value"
      df_upper$type <- paste0(ci_level * 100, "% CI - upper")

      plot_df <- rbind(plot_df, df_lower, df_upper)
    }

    # Reorder facet levels
    plot_df$type <- factor(
      plot_df$type,
      levels = c(
        "Predicted occupancy (mean)",
        "Standard error",
        paste0(ci_level * 100, "% CI - lower"),
        paste0(ci_level * 100, "% CI - upper")
      )
    )

    # Create ggplot
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
        axis.text.y = ggplot2::element_text(angle = 90, vjust = 0, hjust = 0.5)
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
  }

  # Output
  result <- list(mean = psi_mean, se = psi_se)
  if (plot_ci) {
    result$lower <- psi_lower
    result$upper <- psi_upper
  }

  return(invisible(result))
}
