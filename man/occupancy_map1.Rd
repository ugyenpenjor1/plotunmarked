% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occupancy_map1.R
\name{occupancy_map1}
\alias{occupancy_map1}
\title{Generate occupancy prediction map from a single unmarked model}
\usage{
occupancy_map1(
  fitted_model,
  pred_cov_stack,
  agg_factor = 4,
  plot_ci = FALSE,
  save_plot = FALSE,
  plot_filename = "occupancy_map",
  plot_format = "png",
  mask_raster = NULL
)
}
\arguments{
\item{fitted_model}{An \code{unmarkedFit} model object fitted with \code{unmarked}.}

\item{pred_cov_stack}{A \code{terra::SpatRaster} or \code{raster::RasterLayer} containing covariate layers named to match those used in the model. Raster inputs from the \code{raster} package are automatically converted.}

\item{agg_factor}{Integer factor to aggregate (resample) the covariate raster to coarser resolution for faster prediction. Default is 4.}

\item{plot_ci}{Logical indicating whether to include 95\% confidence intervals in the plot. Default is \code{FALSE}.}

\item{save_plot}{Logical indicating whether to save the occupancy plot to disk. Default is \code{FALSE}.}

\item{plot_filename}{Character string filename (without extension) for saving the plot. Default is \code{"occupancy_map"}.}

\item{plot_format}{Character string for the file format (e.g., \code{"png"}, \code{"pdf"}). Default is \code{"png"}.}

\item{mask_raster}{Optional \code{terra::SpatRaster} used to mask or clip the output raster. Default is \code{NULL}.}
}
\value{
A (invisible) list with components:
\item{pred_df}{A data frame containing predicted occupancy, standard error, confidence intervals, and covariate values.}
\item{pred_raster}{A \code{terra::SpatRaster} raster of predicted occupancy probabilities.}
\item{pred_values}{A data frame of raw prediction results including standard errors and confidence limits.}
}
\description{
This function predicts occupancy probabilities from a single fitted \code{unmarked}
occupancy model over a spatial covariate stack.
It aggregates input covariates, computes predictions in chunks for efficiency,
and plots the results with options for confidence intervals and saving the plot.
}
\details{
The function outputs a list containing prediction data, a raster of predicted occupancy,
and the raw predicted values with standard errors and confidence limits.

The function automatically constructs polynomial and interaction terms required by the model formula
and predicts occupancy probabilities chunk-wise to improve performance on large datasets.
}
\examples{
\dontrun{
library(unmarked)
library(terra)

# Assume 'umf' is an unmarkedFrame and 'fitted_model' is an unmarkedFit occupancy model
# Assume 'covariate_stack' is a SpatRaster (also accepts RasterLayer object) with layers named as in the model

result <- occupancy_map1(
  fitted_model = fitted_model,
  pred_cov_stack = covariate_stack,
  agg_factor = 4,
  plot_ci = TRUE,
  save_plot = TRUE,
  plot_filename = "occupancy_map",
  plot_format = "png"
)

# View predictions dataframe
head(result$pred_df)
}

}
