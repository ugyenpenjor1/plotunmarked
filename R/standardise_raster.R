
#' Standardise a raster Layer based on Covariate Data
#'
#' This function standardises a raster layer by subtracting the mean and dividing by the standard deviation
#' of a specified covariate from a provided data frame. It supports raster objects from either the
#' `raster` or `terra` packages and can return the result as either class.
#'
#' @param raster_layer A raster object to standardise. Must be a `SpatRaster` (from `terra`) or `RasterLayer` (from `raster`).
#' @param cov_data A data frame containing the covariate data used to compute mean and standard deviation.
#' @param cov_name A character string specifying the name of the covariate column in `cov_data`.
#' @param return_class A character string indicating the class of the returned raster: `"same"` (default; same as input),
#' `"terra"` to return a `SpatRaster`, or `"raster"` to return a `RasterLayer`.
#'
#' @return A standardized raster object of the class specified by `return_class`.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(system.file("ex/logo.tif", package="terra"))
#' cov_df <- data.frame(elevation = runif(100, 100, 500))
#' std_r <- standardise_raster(r, cov_df, "elevation")
#' plot(std_r)
#' }
#'
#' @export

standardise_raster <- function(raster_layer, cov_data, cov_name, return_class = c("same", "terra", "raster")) {
  # Match return class
  return_class <- match.arg(return_class)

  # Validate raster class
  if (!inherits(raster_layer, c("SpatRaster", "RasterLayer"))) {
    stop("raster_layer must be either a SpatRaster (terra) or RasterLayer (raster) object")
  }

  # Validate covariate data
  if (!is.data.frame(cov_data)) {
    stop("cov_data must be a data.frame")
  }
  if (!(cov_name %in% names(cov_data))) {
    stop(paste("Covariate", cov_name, "not found in cov_data"))
  }

  # Compute standardisation values
  mean_val <- mean(cov_data[[cov_name]], na.rm = TRUE)
  sd_val <- sd(cov_data[[cov_name]], na.rm = TRUE)
  if (sd_val == 0) stop("Standard deviation is zero. Cannot standardize.")

  # Standardize based on type
  if (inherits(raster_layer, "RasterLayer")) {
    # raster package
    standardised <- (raster_layer - mean_val) / sd_val
    if (return_class == "terra") {
      standardised <- terra::rast(standardised)
    }
  } else {
    # terra package
    standardised <- (raster_layer - mean_val) / sd_val
    if (return_class == "raster") {
      standardised <- raster::raster(standardised)
    }
  }

  # Return result
  return(standardised)
}
