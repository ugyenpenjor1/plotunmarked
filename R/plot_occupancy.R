#' Plot occupancy estimates against a continuous covariate
#'
#' This function generates a plot of occupancy probabilities predicted by a fitted
#' \code{unmarkedFit} model, as a function of a selected site-level covariate.
#' A ribbon shows the 95% confidence interval around the predicted values.
#'
#' @param model An object of class \code{unmarkedFit}, typically the result of fitting a model using the \code{unmarked} package.
#' @param covariate A character string giving the name of the site-level covariate to be plotted on the x-axis.
#' @param fixed_vals Optional named list of fixed values for other site-level covariates not being plotted. If not provided, they default to 0.
#' @param xlab Optional custom label for the x-axis. If not specified, defaults to the covariate name with " (standardised)" appended.
#'
#' @return A \code{ggplot} object showing the predicted occupancy probabilities with confidence intervals.
#' @export

plot_occupancy <- function(model, covariate, fixed_vals = NULL, xlab = NULL, ci_level = 0.95) {
  if (!inherits(model, "unmarkedFit")) stop("Model must be of class 'unmarkedFit'")
  if (!is.character(covariate) || length(covariate) != 1) stop("covariate must be a single character string")

  model_data <- model@data@siteCovs

  if (!(covariate %in% colnames(model_data))) {
    stop(paste("Covariate", covariate, "not found in model's site covariates"))
  }

  cov_range <- range(model_data[[covariate]], na.rm = TRUE)
  pred_vals <- seq(cov_range[1], cov_range[2], by = 0.1)
  newdata <- data.frame(pred_vals)

  # Rename to match covariate
  newdata[[covariate]] <- newdata$pred_vals
  newdata$pred_vals <- NULL

  # Fix other covariates
  other_covs <- setdiff(colnames(model_data), covariate)
  for (cov in other_covs) {
    newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
      fixed_vals[[cov]]
    } else {
      0
    }
  }

  # Ask for desired CI level
  occ_prob <- predict(model, type = "state", newdata = newdata, level = ci_level)
  pred <- as.data.frame(occ_prob)

  # Add covariate column back to pred (from newdata)
  pred[[covariate]] <- newdata[[covariate]]

  cov_sym <- rlang::sym(covariate)
  x_label <- if (!is.null(xlab)) xlab else paste0(covariate, " (standardised)")

  p <- ggplot2::ggplot(pred) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = !!cov_sym, ymin = lower, ymax = upper),
      fill = "#56b4e9", alpha = 0.4
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = !!cov_sym, y = Predicted),
      colour = "#0072b2", linewidth = 1.2
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Probability of occupancy"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 15),
      axis.title = ggplot2::element_text(size = 17, vjust = 0.8),
      legend.text = ggplot2::element_text(size = 25),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = 'black')
    )

  return(p)
}
