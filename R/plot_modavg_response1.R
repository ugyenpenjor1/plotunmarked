#' Plot model-averaged response curve
#'
#' Generates a model-averaged response curve for a single covariate based on a
#' set of occupancy or detection models, selected using a ΔAIC threshold.
#'
#' @param model_selection_table A model selection table, typically the result of `modSel()`
#'   from the `unmarked` package. Can also be of class `unmarkedModSel`.
#' @param full_model_list A full `fitList` object containing all candidate models (not pre-filtered).
#' @param covariate Character. The name of the covariate for which the response curve will be plotted.
#' @param fixed_vals Named list. Fixed values for other covariates not being plotted.
#'   If not provided, these covariates are set to zero by default.
#' @param response_type Character. Whether to plot response for `"state"` (occupancy) or `"det"` (detection).
#'   Default is `"state"`.
#' @param xlab Character. Custom label for the x-axis. If `NULL`, the function uses the covariate name.
#' @param ci_level Numeric. Confidence interval level (default is 0.95).
#' @param delta_cutoff Numeric. The ΔAIC threshold for selecting models to include in model averaging (default is 2).
#' @param suppress_warnings Logical. If `TRUE`, warnings about unspecified covariates being set to zero are suppressed.
#'
#' @return A `ggplot2` object showing the model-averaged predicted response curve, with confidence intervals.
#'
#' @details
#' This function computes model-averaged predictions from an `unmarkedFitList` of occupancy models.
#' It supports plotting either the state (occupancy probability) or detection probability (`type = "det"`)
#' for a specified covariate. All other covariates are either held constant at user-specified values (`fixed_vals`)
#' or at zero.
#'
#' For detection probability curves, the covariate can be either site-level or observation-level.
#' If it is observation-level, the function constructs a design matrix with consistent values across surveys.
#'
#' Models are filtered based on ΔAIC cutoff, and model-averaged predictions are weighted by relative AIC weights.
#'
#' @importFrom unmarked predict
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line labs ylim coord_cartesian theme_bw theme element_text element_line element_blank
#' @importFrom rlang sym
#'
#' @examples
#' \dontrun{
#' plot_modavg_response1(
#'   model_selection_table = ms,
#'   full_model_list = mod_ls,
#'   covariate = "cov1",
#'   fixed_vals = list(cov2 = 0, cov3 = 0),
#'   response_type = "state"
#' )
#' }
#'
#' @export

plot_modavg_response1 <- function(
    model_selection_table,
    full_model_list,
    covariate,
    fixed_vals = NULL,
    response_type = c("state", "det"),
    xlab = NULL,
    ci_level = 0.95,
    delta_cutoff = 2,
    suppress_warnings = FALSE
  ) {
  # Match the response type argument
  response_type <- match.arg(response_type)

  # If input is unmarkedModSel object, extract the AIC table dataframe
  if (inherits(model_selection_table, "unmarkedModSel")) {
    model_selection_table <- model_selection_table@Full
  }

  # Ensure model_selection_table has required columns 'delta' and 'model'
  if (!("delta" %in% colnames(model_selection_table)))
    stop("model_selection_table must contain a 'delta' column")
  if (!("model" %in% colnames(model_selection_table)))
    stop("model_selection_table must contain a 'model' column")

  # Filter models based on delta AIC cutoff (default 2)
  ms_sel <- model_selection_table[model_selection_table$delta <= delta_cutoff, ]
  if (nrow(ms_sel) == 0)
    stop("No models within delta AIC ≤ ", delta_cutoff, ".")

  # Extract selected model names and compute model weights from delta AIC
  model_names_sel <- ms_sel$model
  model_weights <- exp(-0.5 * ms_sel$delta)
  model_weights <- model_weights / sum(model_weights)  # normalise weights

  # Extract the fitted models from the full model list by matching model names
  model_list_fits <- full_model_list@fits
  models_to_use_fits <- model_list_fits[model_names_sel]

  # Create a new unmarkedFitList object for selected models
  models_to_use <- methods::new("unmarkedFitList", fits = models_to_use_fits)

  # Extract data object from the first model (assumed all same data)
  umf_data <- models_to_use@fits[[1]]@data
  site_covs <- umf_data@siteCovs   # Site-level covariates
  obs_covs <- umf_data@obsCovs     # Observation-level covariates

  # Check if covariate is in site or observation covariates
  cov_in_site <- covariate %in% colnames(site_covs)
  cov_in_obs <- covariate %in% colnames(obs_covs)

  if (response_type == "state") {
    # For state (occupancy) prediction, covariate must be site-level
    if (!cov_in_site) {
      stop(paste("Covariate", covariate, "not found in site covariates, required for response_type = 'state'"))
    }
    # Generate sequence of values over covariate range for prediction
    cov_range <- range(site_covs[[covariate]], na.rm = TRUE)
    pred_vals <- seq(cov_range[1], cov_range[2], length.out = 100)

    # Prepare newdata for prediction: start with covariate sequence
    newdata <- data.frame(pred_vals)
    newdata[[covariate]] <- newdata$pred_vals
    newdata$pred_vals <- NULL

    # Identify other site covariates (non-target) that need fixed values
    other_covs <- setdiff(colnames(site_covs), covariate)
    # Warn if any non-target covariates are missing fixed values
    missing_fixed <- setdiff(other_covs, names(fixed_vals))
    if (length(missing_fixed) > 0 && !suppress_warnings) {
      warning("The following covariates are not specified in 'fixed_vals' and will be set to 0: ",
              paste(missing_fixed, collapse = ", "))
    }

    # Assign fixed values or zero for other covariates
    for (cov in other_covs) {
      newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
        fixed_vals[[cov]]
      } else {
        0
      }
    }

  } else if (response_type == "det") {
    # For detection probability, covariate may be site- or obs-level
    if (!cov_in_site && !cov_in_obs) {
      stop(paste("Covariate", covariate, "not found in site or observation covariates"))
    }

    # Extract covariate values accordingly
    cov_data <- if (cov_in_obs) obs_covs[[covariate]] else site_covs[[covariate]]
    cov_range <- range(cov_data, na.rm = TRUE)
    pred_vals <- seq(cov_range[1], cov_range[2], length.out = 100)

    if (cov_in_obs) {
      # For observation covariates, replicate covariate values across surveys
      n_surveys <- ncol(umf_data@y)
      newdata <- do.call(rbind, lapply(pred_vals, function(val) {
        # Create a single-row dataframe where all surveys have the same covariate value
        data.frame(matrix(val, nrow = 1, ncol = n_surveys))
      }))
      colnames(newdata) <- rep(covariate, ncol(newdata))

      # Add fixed site covariates with user-supplied values or zero
      for (cov in colnames(site_covs)) {
        newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
          fixed_vals[[cov]]
        } else {
          0
        }
      }

    } else {
      # If covariate is site-level for detection, build newdata similarly to state
      newdata <- data.frame(pred_vals)
      newdata[[covariate]] <- newdata$pred_vals
      newdata$pred_vals <- NULL

      other_covs <- setdiff(colnames(site_covs), covariate)
      missing_fixed <- setdiff(other_covs, names(fixed_vals))
      if (length(missing_fixed) > 0 && !suppress_warnings) {
        warning("The following covariates are not specified in 'fixed_vals' and will be set to 0: ",
                paste(missing_fixed, collapse = ", "))
      }

      for (cov in other_covs) {
        newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
          fixed_vals[[cov]]
        } else {
          0
        }
      }
    }
  }

  # Generate model-averaged predictions for each covariate value
  pred_out <- unmarked::predict(
    models_to_use,
    type = response_type,
    newdata = newdata,
    level = ci_level
  )

  # Convert predictions to data frame for plotting
  pred_df <- as.data.frame(pred_out)

  # Add covariate values to prediction dataframe
  pred_df[[covariate]] <- if (cov_in_obs && response_type == "det") {
    rep(pred_vals, each = 1)
  } else {
    newdata[[covariate]]
  }

  cov_sym <- rlang::sym(covariate)  # For tidy eval in ggplot2

  # Set x-axis label if none provided
  if (is.null(xlab)) xlab <- paste0(covariate, " (standardised)")

  # Create ggplot object with ribbon (CI) and line (prediction)
  p <- ggplot2::ggplot(pred_df) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = !!cov_sym, ymin = lower, ymax = upper),
      fill = "#56b4e9",
      alpha = 0.4
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = !!cov_sym, y = Predicted),
      colour = "#0072b2",
      linewidth = 1.2
    ) +
    ggplot2::labs(
      x = xlab,
      y = ifelse(
        response_type == "state",
        "Probability of occupancy",
        "Probability of detection"
      )
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
      axis.line = ggplot2::element_line(colour = "black")
    )

  # Return ggplot object
  return(p)
}
