% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occupancy_map.R
\name{occupancy_map}
\alias{occupancy_map}
\title{Predict and visualise occupancy (map) from an unmarked model}
\usage{
occupancy_map(
  model,
  covariate_rasters,
  return_class = c("same", "terra", "raster"),
  plot_map = TRUE,
  plot_ci = FALSE,
  ci_level = 0.95,
  save_plot = FALSE,
  plot_filename = "occupancy_map",
  plot_format = "png"
)
}
\arguments{
\item{model}{A fitted model of class \code{unmarkedFit}, typically from the \code{unmarked} package.}

\item{covariate_rasters}{A named list of raster covariates (\code{SpatRaster} or \code{RasterLayer}),
where names match the occupancy covariates used in the model. The list must
include all covariates used in the model, excluding the intercept.}

\item{return_class}{Character. Specifies the class of raster returned.
Options are \code{"same"} (default, returns the same class as input),
\code{"terra"} (returns \code{SpatRaster}), or \code{"raster"} (returns \code{RasterLayer}).}

\item{plot_map}{Logical. If \code{TRUE} (default), plots the predicted occupancy map and uncertainty layers.}

\item{plot_ci}{Logical. If \code{TRUE}, plots confidence interval bounds on the occupancy scale.}

\item{ci_level}{Numeric. Confidence level for intervals (default is \code{0.95} for 95\% CI).}

\item{save_plot}{Logical. If \code{TRUE}, saves the plot to a file using the specified format and filename.}

\item{plot_filename}{Character. Base filename to save the plot (without file extension). Default is \code{"occupancy_map"}.}

\item{plot_format}{Character. File format to use when saving the plot (e.g., \code{"jpeg"}, \code{"png"}). Default is \code{"png"}.}
}
\value{
A (named) list of raster objects containing:
\describe{
\item{\code{mean}}{Raster of predicted occupancy probabilities.}
\item{\code{se}}{Raster of standard errors on the probability scale.}
\item{\code{lower}, \code{upper}}{(If \code{plot_ci = TRUE}) Rasters of lower and upper bounds of the confidence interval.}
}
All returned raster layers match the format specified by \code{return_class}.
}
\description{
Computes and visualises spatial predictions of occupancy probabilities from
a fitted \code{unmarkedFit} model and raster covariates. Supports interaction
terms and returns mean predictions, standard errors and confidence intervals.
Optionally plots the map and saves it to a file.
}
\details{
This function supports main effects and interaction terms of the
form \code{var1:var2}. The function uses model coefficients and their
variance-covariance matrix to generate occupancy predictions across a spatial
surface defined by the input covariate rasters. The intercept term is handled
internally by generating a constant raster. Variance propagation is performed
via the delta method.

Internally, the function computes the linear predictor on the logit scale,
applies the inverse logit to obtain occupancy probabilities, and uses the
delta method to estimate the standard error. Confidence intervals are
computed using the normal approximation if \code{plot_ci = TRUE}.
}
\examples{
\dontrun{
library(unmarked)
library(terra) # also supports library(raster)

# Fit a simple model (example assumes you have data in unmarkedFrame)
umf <- unmarkedFrameOccu(y = matrix(sample(0:1, 100, TRUE), ncol = 5),
                         siteCovs = data.frame(elev = rnorm(20), forest = rnorm(20)))
mod <- occu(~1 ~ elev + forest, data = umf)

# Prepare raster covariates
elev_rast <- rast(matrix(rnorm(100), 10, 10))
forest_rast <- rast(matrix(rnorm(100), 10, 10))
names(elev_rast) <- "elev"
names(forest_rast) <- "forest"

covariate_rasters <- list(elev = elev_rast, forest = forest_rast)

# Generate and plot map
result <- occupancy_map(model = mod, covariate_rasters = covariate_rasters)
}
}
