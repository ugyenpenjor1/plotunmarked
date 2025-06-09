
#' Plot model-averaged response from multiple occupancy models (competing models with some delta threshold)
#'
#' This function creates a ggplot showing model-averaged predictions of occupancy or detection probabilities
#' based on a set of models selected from an AIC model selection table. The function supports prediction across
#' a single covariate while holding others constant.
#'
#' @param model_selection_table An `unmarkedModSel` object or a `data.frame` containing at least 'model' and 'delta' columns.
#' @param model_list An `unmarkedFitList` object containing all candidate models.
#' @param covariate Character string specifying the covariate to vary on the x-axis.
#' @param fixed_vals Optional named list specifying fixed values for other covariates.
#' @param response_type A character string: either `"state"` for occupancy or `"det"` for detection. Defaults to `"state"`.
#' @param xlab Optional x-axis label. If `NULL`, the covariate name is used.
#'
#' @return A `ggplot` object showing the model-averaged response curve with confidence ribbons.
#'
#' @examples
#' \dontrun{
#' library(unmarked)
#' # Assume model_selection_table and model_list are already created
#' plot_modavg_response(model_selection_table, model_list, "elevation")
#' }
#'
#' @export

plot_modavg_response <- function(
    model_selection_table,
    model_list,
    covariate,
    fixed_vals = NULL,
    response_type = c("state", "det"),
    xlab = NULL,
    ci_level = 0.95
) {

  response_type <- match.arg(response_type)

  # Handle model_selection_table from AICcmodavg
  if (inherits(model_selection_table, "unmarkedModSel")) {
    model_selection_table <- model_selection_table@Full
  }

  # Basic checks
  if (!("delta" %in% colnames(model_selection_table)))
    stop("model_selection_table must contain a 'delta' column")
  if (!("model" %in% colnames(model_selection_table)))
    stop("model_selection_table must contain a 'model' column")

  # Select models within delta AIC ≤ 2
  ms_sel <- model_selection_table[model_selection_table$delta <= 2, ]
  if (nrow(ms_sel) == 0)
    stop("No models within delta AIC ≤ 2.")

  model_names_sel <- ms_sel$model
  model_list_fits <- model_list@fits
  models_to_use_fits <- model_list_fits[model_names_sel]
  models_to_use <- methods::new("unmarkedFitList", fits = models_to_use_fits)

  # Access covariates
  umf_data <- models_to_use@fits[[1]]@data
  site_covs <- umf_data@siteCovs
  obs_covs <- umf_data@obsCovs

  # Automatically detect if the covariate is time-varying (obs-level) or site-level
  cov_in_site <- covariate %in% colnames(site_covs)
  cov_in_obs  <- covariate %in% colnames(obs_covs)

  if (response_type == "state") {
    if (!cov_in_site) {
      stop(paste("Covariate", covariate, "not found in site covariates, required for response_type = 'state'"))
    }

    cov_range <- range(site_covs[[covariate]], na.rm = TRUE)
    pred_vals <- seq(cov_range[1], cov_range[2], length.out = 100)
    newdata <- data.frame(pred_vals)
    newdata[[covariate]] <- newdata$pred_vals
    newdata$pred_vals <- NULL

    other_covs <- setdiff(colnames(site_covs), covariate)
    for (cov in other_covs) {
      newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
        fixed_vals[[cov]]
      } else {
        0
      }
    }

  } else if (response_type == "det") {
    if (!cov_in_site && !cov_in_obs) {
      stop(paste("Covariate", covariate, "not found in site or observation covariates"))
    }

    cov_data <- if (cov_in_obs) obs_covs[[covariate]] else site_covs[[covariate]]
    cov_range <- range(cov_data, na.rm = TRUE)
    pred_vals <- seq(cov_range[1], cov_range[2], length.out = 100)

    if (cov_in_obs) {
      # Time-varying covariate: one row per survey
      n_surveys <- ncol(umf_data@y)

      # Construct newdata for time-varying detection covariates
      newdata <- do.call(rbind, lapply(pred_vals, function(val) {
        data.frame(matrix(val, nrow = 1, ncol = n_surveys))
      }))
      colnames(newdata) <- rep(covariate, ncol(newdata))

      # Add dummy site-level covariates (set to 0 or fixed_vals)
      for (cov in colnames(site_covs)) {
        newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
          fixed_vals[[cov]]
        } else {
          0
        }
      }

    } else {
      # Site-level detection covariate
      newdata <- data.frame(pred_vals)
      newdata[[covariate]] <- newdata$pred_vals
      newdata$pred_vals <- NULL

      other_covs <- setdiff(colnames(site_covs), covariate)
      for (cov in other_covs) {
        newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
          fixed_vals[[cov]]
        } else {
          0
        }
      }
    }
  }

  # Predict from model-averaged fit
  pred_out <- unmarked::predict(
    models_to_use,
    type = response_type,
    newdata = newdata,
    level = ci_level
  )
  pred_df <- as.data.frame(pred_out)

  # Assign the covariate to the plotting data frame
  pred_df[[covariate]] <- if (cov_in_obs && response_type == "det") {
    rep(pred_vals, each = 1)  # 1 row per pred_val (since 1 site with full survey columns)
  } else {
    newdata[[covariate]]
  }

  # Plot
  cov_sym <- rlang::sym(covariate)
  if (is.null(xlab))
    xlab <- paste0(covariate, " (standardised)")

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
        response_type == "state", "Probability of occupancy", "Probability of detection")
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

  return(p)
}
