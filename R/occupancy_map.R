#' Generate occupancy probability maps from an unmarked model
#'
#' This function creates spatial maps of predicted occupancy probabilities
#' based on an `unmarkedFit` model and a set of covariate rasters. It can
#' return the maps as either terra or raster objects, optionally plot the maps,
#' show confidence intervals, and save the plot to a file.
#'
#' @param model An `unmarkedFit` model object (from the `unmarked` package).
#' @param covariate_rasters A **named list** of raster objects (either `SpatRaster` or `RasterLayer`)
#'   matching the covariates used in the occupancy model. The names should correspond
#'   to the covariate names in the model.
#' @param return_class Character specifying the raster class of output maps.
#'   One of `"same"` (default, match input raster class), `"terra"`, or `"raster"`.
#' @param plot_map Logical; if `TRUE` (default), a map plot is displayed.
#' @param plot_ci Logical; if `TRUE`, confidence interval maps are plotted along with the mean and SE.
#' @param ci_level Numeric; confidence level for intervals (default 0.95).
#' @param save_plot Logical; if `TRUE`, the plot will be saved to disk.
#' @param plot_filename Character; base filename for saving the plot (default `"occupancy_plot"`).
#' @param plot_format Character; format to save the plot (e.g., `"jpeg"`, `"png"`; default `"jpeg"`).
#'
#' @return A **named list** of rasters:
#' \item{mean}{Raster layer of predicted occupancy probabilities.}
#' \item{se}{Raster layer of standard errors.}
#' \item{lower}{(If `plot_ci = TRUE`) Raster layer of lower confidence interval bound.}
#' \item{upper}{(If `plot_ci = TRUE`) Raster layer of upper confidence interval bound.}
#'
#' The result is returned invisibly.
#'
#' @details
#' The function extracts coefficients and variance-covariance matrix from the
#' unmarked occupancy model, combines them with spatial covariate rasters to compute
#' predicted occupancy probabilities on the map. It supports plotting using ggplot2
#' with options to display uncertainty and save plots.
#'
#' @importFrom stats qnorm
#' @importFrom terra rast app values global
#' @importFrom raster raster cellStats
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_viridis_c coord_equal theme_minimal labs facet_wrap theme element_rect element_text ggsave
#' @importFrom pbapply pbapply
#' @importFrom stats setNames coef sd vcov
#'
#' @examples
#' \dontrun{
#' # Assuming `fit` is an unmarkedFit model and `covs` a named list of rasters
#' result <- occupancy_map(fit, covs, plot_map = TRUE, plot_ci = TRUE, save_plot = TRUE)
#' }
#'
#' @export

occupancy_map <- function(
    model,
    covariate_rasters,
    return_class = c("same", "terra", "raster"),
    plot_map = TRUE,
    plot_ci = FALSE,
    ci_level = 0.95,
    save_plot = FALSE,
    plot_filename = "occupancy_plot",
    plot_format = "jpeg"
) {
  return_class <- match.arg(return_class)

  if (!inherits(model, "unmarkedFit")) {
    stop("`model` must be an 'unmarkedFit' object.")
  }

  if (!is.list(covariate_rasters) || is.null(names(covariate_rasters))) {
    stop("`covariate_rasters` must be a named list of rasters with names matching model covariates.")
  }

  # Extract coefficients and VCV matrix
  beta <- coef(model)
  vcov_mat <- vcov(model)

  psi_coefs <- beta[grep("^psi\\(", names(beta))]
  names_clean <- gsub("^psi\\((.*)\\)$", "\\1", names(psi_coefs))
  names(psi_coefs) <- names_clean

  coef_names <- names(psi_coefs)
  beta_vec <- psi_coefs

  # Extract relevant VCV matrix
  beta_idx <- grep("^psi\\(", rownames(vcov_mat))
  V_sub <- vcov_mat[beta_idx, beta_idx]
  rownames(V_sub) <- colnames(V_sub) <- names_clean

  # Create raster stack of covariates
  X_stack <- terra::rast(lapply(coef_names, function(cov) {
    if (cov == "Int") {
      base_raster <- covariate_rasters[[1]]
      return(base_raster * 0 + 1)
    } else {
      rast <- covariate_rasters[[cov]]
      if (!inherits(rast, c("SpatRaster", "RasterLayer"))) {
        stop(sprintf("Raster for '%s' must be SpatRaster or RasterLayer.", cov))
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
      if (raster::cellStats(r, stat = function(x) sum(!is.na(x))) == 0) {
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
    # Convert rasters to data frames
    df_mean <- as.data.frame(psi_mean, xy = TRUE, na.rm = TRUE)
    names(df_mean)[3] <- "value"
    df_mean$type <- "Predicted occupancy (mean)"

    df_se <- as.data.frame(psi_se, xy = TRUE, na.rm = TRUE)
    names(df_se)[3] <- "value"
    df_se$type <- "Standard error"

    plot_df <- rbind(df_mean, df_se)

    if (plot_ci) {
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
        option = "magma",
        na.value = "transparent",
        imits = c(0, 1)
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
