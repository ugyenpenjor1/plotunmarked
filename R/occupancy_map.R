
#' Predict and map occupancy probabilities from a fitted unmarked model
#'
#' This function computes predicted occupancy probabilities across a spatial landscape
#' using raster covariates and a fitted `unmarkedFit` model. The result is a raster
#' layer showing spatial variation in occupancy probability.
#'
#' @param model A fitted `unmarkedFit` object (e.g., from `occu()`).
#' @param covariate_rasters A named list of raster covariates (as `SpatRaster` or `RasterLayer`),
#'   where each name matches a covariate used in the model.
#' @param return_class Desired class of the returned raster: `"same"` (default), `"terra"`, or `"raster"`.
#' @param plot_map Logical; if `TRUE`, plots the resulting occupancy map (default is `TRUE`).
#'
#' @return A raster layer (`SpatRaster` or `RasterLayer`) representing predicted occupancy probabilities.
#'
#' @importFrom terra global rast
#' @importFrom raster cellStats raster
#' @importFrom methods new
#'
#' @examples
#' # Example usage (requires fitted model and covariate rasters)
#' # result <- occupancy_map(model, covariate_rasters)
#'
#' @export

occupancy_map <- function(model, covariate_rasters, return_class = c("same", "terra", "raster"), plot_map = TRUE) {
  return_class <- match.arg(return_class)

  # Validate model
  if (!inherits(model, "unmarkedFit")) {
    stop("`model` must be an 'unmarkedFit' object.")
  }

  # Validate raster list
  if (!is.list(covariate_rasters) || is.null(names(covariate_rasters))) {
    stop("`covariate_rasters` must be a named list of rasters with names matching model covariates.")
  }

  # Extract and clean occupancy coefficients
  beta <- coef(model)
  psi_coefs <- beta[grep("^psi\\(", names(beta))]
  names_clean <- gsub("^psi\\((.*)\\)$", "\\1", names(psi_coefs))
  names(psi_coefs) <- names_clean

  intercept <- if ("Int" %in% names(psi_coefs)) psi_coefs["Int"] else 0
  covars <- setdiff(names(psi_coefs), "Int")

  # Check missing covariates
  missing_covars <- setdiff(covars, names(covariate_rasters))
  if (length(missing_covars) > 0) {
    stop("Missing raster(s) for covariate(s): ", paste(missing_covars, collapse = ", "))
  }

  # Check raster values
  for (name in covars) {
    rast <- covariate_rasters[[name]]
    if (!inherits(rast, c("SpatRaster", "RasterLayer"))) {
      stop(sprintf("Raster for '%s' must be SpatRaster or RasterLayer.", name))
    }

    # Check if raster has values
    if (inherits(rast, "SpatRaster")) {
      has_vals <- terra::global(rast, fun = function(x) sum(!is.na(x)))[1, 1]
    } else {
      has_vals <- raster::cellStats(rast, stat = function(x) sum(!is.na(x)))
    }

    if (is.na(has_vals) || has_vals == 0) {
      stop(sprintf("Raster for '%s' has no non-NA values.", name))
    }
  }

  # Start with intercept
  base_raster <- covariate_rasters[[covars[1]]]
  logitPsi <- base_raster * 0 + intercept

  # Apply covariate rasters
  for (cov in covars) {
    logitPsi <- logitPsi + psi_coefs[cov] * covariate_rasters[[cov]]
  }

  # Convert to probability
  psi <- exp(logitPsi) / (1 + exp(logitPsi))

  # Convert format if needed
  if (return_class == "raster" && inherits(psi, "SpatRaster")) {
    psi <- raster::raster(psi)
  } else if (return_class == "terra" && inherits(psi, "RasterLayer")) {
    psi <- terra::rast(psi)
  }

  if (plot_map) {
    plot(psi, main = "Predicted occupancy probability")
  }

  return(psi)
}
