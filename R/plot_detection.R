#' Plot detection probability
#'
#' This function plots predicted detection probabilities from an `unmarkedFit` model
#' based on a single covariate. It handles both site-level and observation-level
#' covariates, generating predictions while keeping other covariates fixed at specified
#' or default values.
#'
#' @param model An object of class `unmarkedFit`, typically from `unmarked::occu()`.
#' @param covariate A character string indicating the name of the covariate to plot.
#' @param fixed_vals Optional named list of fixed values for other covariates.
#' @param xlab Optional custom x-axis label.
#'
#' @return A `ggplot2` plot of predicted detection probability.
#' @examples
#' # Example usage (assumes a fitted model and the 'unmarked' package):
#' # plot_detection(model, "wind")
#'
#' @export

plot_detection <- function(model, covariate, fixed_vals = NULL, xlab = NULL) {
  if (!inherits(model, "unmarkedFit")) stop("Model must be of class 'unmarkedFit'")
  if (!is.character(covariate) || length(covariate) != 1) stop("covariate must be a single character string")

  site_covs <- model@data@siteCovs
  obs_covs <- model@data@obsCovs
  cov_in_site <- covariate %in% names(site_covs)
  cov_in_obs  <- !is.null(obs_covs) && covariate %in% names(obs_covs)

  if (!cov_in_site && !cov_in_obs) {
    stop(paste("Covariate", covariate, "not found in site or observation covariates"))
  }

  # If covariate is site-level
  if (cov_in_site) {
    cov_range <- range(site_covs[[covariate]], na.rm = TRUE)
    pred_vals <- seq(cov_range[1], cov_range[2], by = 0.1)
    newdata <- data.frame(pred_vals)
    names(newdata) <- covariate

    other_covs <- setdiff(names(site_covs), covariate)
    for (cov in other_covs) {
      newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
        fixed_vals[[cov]]
      } else {
        0
      }
    }

    # Remove row names to avoid warning in predict
    rownames(newdata) <- NULL

    pred <- unmarked::predict(model, type = "det", newdata = newdata, appendData = TRUE)
    pred_df <- as.data.frame(pred)
    cov_sym <- rlang::sym(covariate)

  } else {
    # covariate is in obsCovs — time-varying
    cov_range <- range(obs_covs[[covariate]], na.rm = TRUE)
    pred_vals <- seq(cov_range[1], cov_range[2], by = 0.1)

    n_sites <- nrow(model@data@y)
    n_occasions <- ncol(model@data@y)
    total_obs <- n_sites * n_occasions

    obs_df <- data.frame(matrix(NA, nrow = total_obs, ncol = 0))
    obs_df[[covariate]] <- rep(pred_vals, length.out = total_obs)

    all_obs_covs <- colnames(obs_covs)
    other_covs <- setdiff(all_obs_covs, covariate)
    for (cov in other_covs) {
      obs_df[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
        fixed_vals[[cov]]
      } else {
        0
      }
    }

    # Remove row names to avoid warning in predict
    rownames(obs_df) <- NULL

    # Create dummy site covs
    site_df <- as.data.frame(site_covs[1, , drop = FALSE])
    for (col in names(site_df)) site_df[[col]] <- 0
    site_df <- site_df[rep(1, n_sites), ]
    rownames(site_df) <- NULL

    # Create dummy response (not used in prediction)
    y_dummy <- matrix(1, nrow = n_sites, ncol = n_occasions)
    rownames(y_dummy) <- NULL  # Though this is a matrix, just to be consistent

    umf_pred <- unmarked::unmarkedFrameOccu(
      y = y_dummy,
      siteCovs = site_df,
      obsCovs = obs_df
    )

    # Remove row names from umf_pred components
    rownames(umf_pred@y) <- NULL
    rownames(umf_pred@siteCovs) <- NULL
    rownames(umf_pred@obsCovs) <- NULL

    pred <- unmarked::predict(model, type = "det", newdata = umf_pred, appendData = TRUE)
    pred_df <- as.data.frame(pred)
    pred_df[[covariate]] <- obs_df[[covariate]]
    cov_sym <- rlang::sym(covariate)
  }

  # Define x-axis label
  x_label <- if (!is.null(xlab)) xlab else paste0(covariate, " (standardised)")

  # Plot
  p <- ggplot2::ggplot(pred_df) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = !!cov_sym, ymin = lower, ymax = upper),
      fill = "#56b4e9", alpha = 0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = !!cov_sym, y = Predicted),
      colour = "#0072b2", linewidth = 1.2
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Probability of detection"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 15),
      axis.title = ggplot2::element_text(size = 17, vjust = 0.8),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = 'black')
    )

  return(p)
}
