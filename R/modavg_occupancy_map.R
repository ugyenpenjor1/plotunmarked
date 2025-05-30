
#' Model-averaged occupancy prediction map
#'
#' Generates a spatial map of model-averaged occupancy probabilities using a
#' list of fitted occupancy models and a raster stack of covariates.
#'
#' @param d2_mod_list A candidate model set, typically created using
#'   \code{AICcmodavg::aictab} or as a list of fitted unmarked models.
#' @param pred_cov_stack A \code{terra::SpatRaster} with named layers
#'   corresponding to the covariates used in the models.
#' @param agg_factor An integer specifying the aggregation factor to reduce
#'   raster resolution before prediction. Default is 4.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{pred_df}}{A data frame containing predictions and coordinates.}
#'   \item{\code{pred_raster}}{A \code{terra::SpatRaster} object of predicted occupancy.}
#'   \item{\code{mod_avg_pred}}{A list containing predicted values, standard errors,
#'   and confidence intervals.}
#' }
#'
#' @details The function splits the raster prediction task into manageable chunks,
#' making it memory-efficient for large areas. It uses
#' \code{AICcmodavg::modavgPred()} for model-averaged predictions. A progress bar
#' from \code{progress} package is displayed during execution.
#'
#' @importFrom terra aggregate rast plot
#' @importFrom AICcmodavg modavgPred
#' @importFrom wesanderson wes_palette
#' @importFrom progress progress_bar
#'
#' @examples
#' \dontrun{
#' # Assuming 'models' is a list of fitted unmarked models,
#' # and 'cov_stack' is a SpatRaster of covariates
#' result <- modavg_occupancy_map(models, cov_stack, agg_factor = 4)
#' }
#'
#' @export

modavg_occupancy_map <- function(
    d2_mod_list,
    pred_cov_stack,
    agg_factor = 4
) {
  if (is.null(names(pred_cov_stack))) stop("pred_cov_stack must have layer names matching model covariates.")

  pred_cov_reduced <- terra::aggregate(pred_cov_stack, fact = agg_factor, fun = mean)
  newdata <- as.data.frame(pred_cov_reduced, xy = TRUE)

  n <- nrow(newdata)

  # Automatically set chunk size
  if (n <= 1000) {
    chunk_size <- n      # no chunking for small data
  } else if (n <= 5000) {
    chunk_size <- 500
  } else {
    chunk_size <- 1000
  }

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

  xyz_mat <- cbind(
    x = pred_df$x,
    y = pred_df$y,
    z = pred_df$Predicted
  )

  pred_raster <- terra::rast(xyz_mat, type = "xyz")

  terra::plot(
    pred_raster,
    col = wesanderson::wes_palette("Zissou1", 100, type = "continuous"),
    main = "Model-averaged occupancy prediction"
  )

  return(list(
    pred_df = pred_df,
    pred_raster = pred_raster,
    mod_avg_pred = mod_avg_preds
  ))
}
