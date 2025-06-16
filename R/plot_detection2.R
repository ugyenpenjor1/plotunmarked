
#' Plot detection probability vs. covariate from unmarked model
#'
#' Creates a plot of predicted detection probabilities from a fitted `unmarked` model
#' as a function of a specified covariate. Supports both continuous and categorical
#' covariates and polynomial terms.
#'
#' @param model An object of class `unmarkedFit`, typically created by fitting a model
#'   using functions from the `unmarked` package such as `pcount`, `occu`, etc.
#' @param focal_cov Character string specifying the name of the focal covariate (either
#'   detection or site-level) to plot on the x-axis.
#' @param cat_cov Logical. If `TRUE`, the focal covariate is treated as categorical.
#'   Defaults to `FALSE`.
#' @param group_cov Optional character string specifying the name of a grouping
#'   categorical covariate. If `NULL` and `cat_cov = TRUE`, defaults to `focal_cov`.
#' @param xlab Optional character string for labeling the x-axis. If `NULL`, uses
#'   `focal_cov` as the label.
#'
#' @return A `ggplot2` object showing the predicted detection probability as a function
#'   of the specified covariate, with optional grouping.
#'
#' @details
#' This function uses the detection formula of the `unmarkedFit` model to
#' generate predicted detection probabilities over the range of a focal covariate.
#' When the covariate is continuous, a sequence of values is generated for plotting.
#' If the covariate is categorical, the function plots estimates for each category.
#' If a second grouping variable is supplied, separate lines or error bars are
#' generated for each group.
#'
#' If the model includes a quadratic term (e.g., `I(x^2)`) for the focal covariate,
#' this term is automatically constructed and included in prediction.
#'
#' All other covariates are held constant at their reference levels (for factors) or
#' at 0 (for numeric covariates).
#'
#' @import ggplot2
#' @importFrom unmarked predict
#'
#' @examples
#' \dontrun{
#' # Assuming `umf` is an unmarkedFrame and `fit` is a model fit with occu()
#' fit <- occu(~date+I(date^2)+observer+wind~forest+elevation, data =  umf)
#' plot_detection2(fit, focal_cov = "wind", cat_cov = FALSE)
#' plot_detection2(fit, focal_cov = "observer", cat_cov = TRUE)
#' plot_detection2(fit, focal_cov = "date", group_cov = NULL, cat_cov = FALSE)
#' }
#'
#' @export


plot_detection2 <- function(
    model,
    focal_cov,
    cat_cov = FALSE,
    group_cov = NULL,
    xlab = NULL
) {

  # Detection and site covariates
  det_covs <- model@data@obsCovs
  site_covs <- model@data@siteCovs

  # Helper function to check if covariate is detection or site
  get_cov_source <- function(cov) {
    if (cov %in% names(det_covs)) return("obsCovs")
    if (cov %in% names(site_covs)) return("siteCovs")
    stop(paste("Covariate", cov, "not found in detection or site covariates."))
  }

  # Check focal covariate exists in either obsCovs or siteCovs
  focal_source <- get_cov_source(focal_cov)

  # If group_cov not provided, default to focal_cov if categorical
  if (is.null(group_cov) & cat_cov) {
    group_cov <- focal_cov
  }

  if (!is.null(group_cov)) {
    group_source <- get_cov_source(group_cov)
  }

  # Extract vectors and levels accordingly
  focal_data <- if (focal_source == "obsCovs") det_covs[[focal_cov]] else site_covs[[focal_cov]]
  if (!is.null(group_cov)) {
    group_data <- if (group_source == "obsCovs") det_covs[[group_cov]] else site_covs[[group_cov]]
  }

  # Create range for focal if continuous
  if (!cat_cov) {
    seq_vals <- seq(
      min(focal_data, na.rm = TRUE),
      max(focal_data, na.rm = TRUE),
      length.out = 51
    )
  }

  # Polynomial term creator (for continuous focal covariate)
  make_poly_terms <- function(df, var) {
    det_form <- model@call$detformula
    if (is.null(det_form)) return(df)

    poly_term <- paste0("I(", var, "^2)")
    if (poly_term %in% attr(terms(det_form), "term.labels")) {
      df[[poly_term]] <- df[[var]]^2
    }
    return(df)
  }

  # Create newdata for prediction
  if (!cat_cov) {
    # Continuous focal covariate

    # If grouping provided, expand grid
    if (!is.null(group_cov)) {
      if (!is.factor(group_data)) {
        stop("Group covariate must be a factor.")
      }

      newdata <- expand.grid(
        focal = seq_vals,
        group = levels(group_data)
      )
      names(newdata) <- c(focal_cov, group_cov)

      # Fill other covariates in obsCovs
      for (cov in names(det_covs)) {
        if (!(cov %in% c(focal_cov, group_cov))) {
          if (is.factor(det_covs[[cov]])) {
            newdata[[cov]] <- factor(
              levels(det_covs[[cov]])[1], # baseline (first in the category)
              levels = levels(det_covs[[cov]])
            )
          } else {
            newdata[[cov]] <- 0
          }
        }
      }

      # Fill other covariates in siteCovs
      for (cov in names(site_covs)) {
        if (!(cov %in% c(focal_cov, group_cov))) {
          if (is.factor(site_covs[[cov]])) {
            newdata[[cov]] <- factor(
              levels(site_covs[[cov]])[1], # baseline (first in the category)
              levels = levels(site_covs[[cov]])
            )
          } else {
            newdata[[cov]] <- 0
          }
        }
      }

      newdata <- make_poly_terms(newdata, focal_cov)

      preds <- unmarked::predict(
        model,
        type = "det",
        newdata = newdata,
        appendData = TRUE
      )

      p <- ggplot2::ggplot(
        preds,
        ggplot2::aes(
          x = .data[[focal_cov]],
          y = .data[["Predicted"]],
          color = .data[[group_cov]],
          group = .data[[group_cov]]
        )
      ) +
        ggplot2::geom_line(linewidth = 1.2) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower, ymax = upper, fill = .data[[group_cov]]),
          alpha = 0.2,
          color = NA
        ) +
        ggplot2::theme_bw() +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(
          x = ifelse(is.null(xlab), focal_cov, xlab),
          y = "Probability of detection",
          color = group_cov,
          fill = group_cov
        ) +
        ggplot2::theme(
          axis.text = ggplot2::element_text(size = 15),
          axis.title = ggplot2::element_text(size = 17, vjust = 0.8),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = 'black')
        )

    } else {
      # No grouping
      newdata <- data.frame(seq_vals)
      names(newdata) <- focal_cov

      # Fill other detection covariates
      for (cov in names(det_covs)) {
        if (cov != focal_cov) {
          if (is.factor(det_covs[[cov]])) {
            newdata[[cov]] <- factor(
              levels(det_covs[[cov]])[1],
              levels = levels(det_covs[[cov]])
            )
          } else {
            newdata[[cov]] <- 0
          }
        }
      }
      # Fill site covariates
      for (cov in names(site_covs)) {
        if (cov != focal_cov) {
          if (is.factor(site_covs[[cov]])) {
            newdata[[cov]] <- factor(
              levels(site_covs[[cov]])[1],
              levels = levels(site_covs[[cov]])
            )
          } else {
            newdata[[cov]] <- 0
          }
        }
      }

      newdata <- make_poly_terms(newdata, focal_cov)

      preds <- unmarked::predict(
        model,
        type = "det",
        newdata = newdata,
        appendData = TRUE
      )

      p <- ggplot2::ggplot(
        preds,
        ggplot2::aes(x = .data[[focal_cov]], y = .data[["Predicted"]])
      ) +
        ggplot2::geom_line(color = "steelblue3", linewidth = 1.2) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower, ymax = upper),
          alpha = 0.25,
          fill = "dodgerblue"
        ) +
        ggplot2::theme_bw() +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(
          x = ifelse(is.null(xlab), focal_cov, xlab),
          y = "Probability of detection"
        ) +
        ggplot2::theme(
          axis.text = ggplot2::element_text(size = 15),
          axis.title = ggplot2::element_text(size = 17, vjust = 0.8),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = 'black')
        )
    }

  } else {
    # Categorical focal covariate

    if (is.null(group_cov)) {
      group_cov <- focal_cov
      group_source <- focal_source
      group_data <- focal_data
    }

    if (!is.factor(group_data)) {
      stop("Group covariate must be a factor.")
    }

    if (group_cov == focal_cov) {
      # just levels of focal covariate

      newdata <- data.frame(focal_levels = levels(focal_data))
      names(newdata) <- focal_cov

      # Fill other detection covariates
      for (cov in names(det_covs)) {
        if (cov != focal_cov) {
          if (is.factor(det_covs[[cov]])) {
            newdata[[cov]] <- factor(
              levels(det_covs[[cov]])[1],
              levels = levels(det_covs[[cov]])
            )
          } else {
            newdata[[cov]] <- 0
          }
        }
      }
      # Fill site covariates
      for (cov in names(site_covs)) {
        if (cov != focal_cov) {
          if (is.factor(site_covs[[cov]])) {
            newdata[[cov]] <- factor(
              levels(site_covs[[cov]])[1],
              levels = levels(site_covs[[cov]])
            )
          } else {
            newdata[[cov]] <- 0
          }
        }
      }

      preds <- unmarked::predict(
        model,
        type = "det",
        newdata = newdata,
        appendData = TRUE
      )

      p <- ggplot2::ggplot(
        preds,
        ggplot2::aes(
          x = .data[[focal_cov]],
          y = .data[["Predicted"]],
          group = .data[[focal_cov]],
          #color = .data[[focal_cov]]
          color = NA
        )
      ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = lower, ymax = upper),
          width = 0,
          linewidth = 1.2,
          colour = "steelblue"
        ) +
        ggplot2::geom_point(size = 4.5, colour = "maroon") +
        ggplot2::geom_line(linewidth = 1.2) +
        ggplot2::theme_bw() +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(
          x = ifelse(is.null(xlab), focal_cov, xlab),
          y = "Probability of detection",
          color = focal_cov
        ) +
        ggplot2::theme(
          legend.position = "none",
          axis.text = ggplot2::element_text(size = 15),
          axis.title = ggplot2::element_text(size = 17, vjust = 0.8),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = 'black')
        )

    } else {
      # different group covariate

      newdata <- expand.grid(
        focal_levels = levels(focal_data),
        group_levels = levels(group_data)
      )
      names(newdata) <- c(focal_cov, group_cov)

      # Fill other detection covariates
      for (cov in names(det_covs)) {
        if (!(cov %in% c(focal_cov, group_cov))) {
          if (is.factor(det_covs[[cov]])) {
            newdata[[cov]] <- factor(
              levels(det_covs[[cov]])[1],
              levels = levels(det_covs[[cov]])
            )
          } else {
            newdata[[cov]] <- 0
          }
        }
      }
      # Fill site covariates
      for (cov in names(site_covs)) {
        if (!(cov %in% c(focal_cov, group_cov))) {
          if (is.factor(site_covs[[cov]])) {
            newdata[[cov]] <- factor(
              levels(site_covs[[cov]])[1],
              levels = levels(site_covs[[cov]])
            )
          } else {
            newdata[[cov]] <- 0
          }
        }
      }

      preds <- unmarked::predict(
        model,
        type = "det",
        newdata = newdata,
        appendData = TRUE
      )

      p <- ggplot2::ggplot(
        preds,
        ggplot2::aes(
          x = .data[[focal_cov]],
          y = .data[["Predicted"]],
          #color = .data[[group_cov]],
          color = NA,
          group = .data[[group_cov]]
        )
      ) +
        ggplot2::geom_point(
          size = 3.5,
          position = ggplot2::position_dodge(width = 0.5)
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = lower, ymax = upper),
          width = 0,
          linewidth = 0.8,
          position = ggplot2::position_dodge(width = 0.5)
        ) +
        ggplot2::geom_line(
          position = ggplot2::position_dodge(width = 0.5),
          size = 1.2
        ) +
        ggplot2::theme_bw() +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(
          x = ifelse(is.null(xlab), focal_cov, xlab),
          y = "Probability of detection",
          color = group_cov
        ) +
        ggplot2::theme(
          legend.position = "none",
          axis.text = ggplot2::element_text(size = 15),
          axis.title = ggplot2::element_text(size = 17, vjust = 0.8),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = 'black')
        )
    }
  }

  return(p)
}
