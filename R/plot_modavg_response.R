
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
    model_selection_table,    # unmarkedModSel object OR data.frame of model selection results
    model_list,               # unmarkedFitList object containing all candidate models
    covariate,                # covariate name to vary on x-axis
    fixed_vals = NULL,        # named list of fixed values for other covariates
    response_type = c("state", "det"),  # which response to predict: occupancy ("state") or detection ("det")
    xlab = NULL               # optional x-axis label; if NULL defaults to covariate name
) {

  response_type <- match.arg(response_type)

  # If user passed full unmarkedModSel object, extract @Full slot (data.frame)
  if (inherits(model_selection_table, "unmarkedModSel")) {
    model_selection_table <- model_selection_table@Full
  }

  # Check required columns
  if (!("delta" %in% colnames(model_selection_table))) stop("model_selection_table must contain a 'delta' column")
  if (!("model" %in% colnames(model_selection_table))) stop("model_selection_table must contain a 'model' column with model names")

  # Select models within delta AIC 2 (or only 1 if less than 2 models qualify)
  ms_sel <- model_selection_table[model_selection_table$delta <= 2, ]
  if (nrow(ms_sel) == 0) stop("No models within delta AIC 2")

  model_names_sel <- ms_sel$model

  # Extract list of models from unmarkedFitList (slot 'fits')
  model_list_fits <- model_list@fits

  # Subset the fits list to selected models
  models_to_use_fits <- model_list_fits[model_names_sel]

  # Wrap back into unmarkedFitList
  models_to_use <- methods::new("unmarkedFitList", fits = models_to_use_fits)

  # Extract site-level covariates from first model (assuming all have same covariates)
  site_covs <- models_to_use@fits[[1]]@data@siteCovs
  if (!(covariate %in% colnames(site_covs))) {
    stop(paste("Covariate", covariate, "not found in site covariates"))
  }

  # Generate sequence for covariate
  cov_range <- range(site_covs[[covariate]], na.rm = TRUE)
  pred_vals <- seq(cov_range[1], cov_range[2], by = 0.1) # step size for covariate sequence

  # Build newdata for prediction
  newdata <- data.frame(pred_vals)
  names(newdata) <- covariate

  # Fill other covariates with fixed values or zero
  other_covs <- setdiff(colnames(site_covs), covariate)
  for (cov in other_covs) {
    newdata[[cov]] <- if (!is.null(fixed_vals) && cov %in% names(fixed_vals)) {
      fixed_vals[[cov]]
    } else {
      0
    }
  }

  # Predict model-averaged response using unmarked::predict on unmarkedFitList
  pred_out <- unmarked::predict(models_to_use, type = response_type, newdata = newdata, appendData = TRUE)
  pred_df <- as.data.frame(pred_out)

  # Prepare covariate symbol for ggplot
  cov_sym <- rlang::sym(covariate)

  # Default x-axis label if none provided
  if (is.null(xlab)) xlab <- paste0(covariate, " (standardised)")

  # Plot
  p <- ggplot2::ggplot(pred_df) +
    ggplot2::geom_ribbon(ggplot2::aes(x = !!cov_sym, ymin = lower, ymax = upper), fill = "#56b4e9", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x = !!cov_sym, y = Predicted), colour = "#0072b2", size = 1.2) +
    ggplot2::labs(
      x = xlab,
      y = ifelse(
        response_type == "state", "Probability of occupancy", "Probability of detection"
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

  return(p)
}
