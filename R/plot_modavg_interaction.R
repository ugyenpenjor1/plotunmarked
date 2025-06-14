#' Plot model-averaged occupancy interaction effects
#'
#' This function creates a plot of model-averaged occupancy probabilities
#' across a range of two interacting covariates using a set of models selected
#' based on an AIC threshold.
#'
#' @param model_selection_table A model selection table, typically generated by `modSel()` from the `unmarked` package.
#'   Can also be of class `unmarkedModSel`.
#' @param full_model_list A full `fitList` of all candidate models (not pre-filtered by AIC).
#' @param cov1 Character. The first covariate involved in the interaction (x-axis of the plot).
#' @param cov2 Character. The second covariate involved in the interaction (used to color lines).
#' @param delta_cutoff Numeric. Threshold for selecting models within a given ΔAIC (default is 2).
#' @param grid_length Integer. Number of points to use in the prediction grid for each covariate (default is 400).
#' @param chunk_size Integer. Number of rows to predict at once (useful for memory management; default is 50,000).
#' @param palette Character. Optional name of a color palette passed to `hcl.colors()`. Defaults to `"cividis"`.
#' @param xlab Character. Optional label for the x-axis. If `NULL`, `cov1` is used.
#' @param nonInt_cov Named list. Values for non-interacting covariates to be held constant during prediction.
#'   If not specified, these covariates are set to their median value from the data.
#'
#' @return A `ggplot2` object showing the model-averaged relationship between occupancy and the interaction
#'   between `cov1` and `cov2`, with lines colored by `cov2`.
#'
#' @details
#' Only models containing both `cov1` and `cov2` as site-level covariates will be used. The function generates
#' a grid of predictions across a combination of values of `cov1` and `cov2`, and computes model-averaged
#' occupancy probabilities weighted by AIC support.
#'
#' Any other site-level covariates not specified in `nonInt_cov` will be set to their median values.
#'
#' @importFrom unmarked predict
#' @importFrom ggplot2 ggplot aes geom_line scale_colour_gradientn labs theme_bw theme element_text element_line element_blank
#' @importFrom progress progress_bar
#' @importFrom grDevices hcl.colors
#'
#' @examples
#' \dontrun{
#' plot_modavg_interaction(
#'   model_selection_table = ms,
#'   full_model_list = mod_ls,
#'   cov1 = "cov1",
#'   cov2 = "cov2",
#'   delta_cutoff = 2,
#'   grid_length = 400,
#'   nonInt_cov = list(cov3 = 0, cov4 = 0)
#' )
#' }
#'
#' @export

plot_modavg_interaction <- function(model_selection_table,
                                     full_model_list,
                                     cov1,
                                     cov2,
                                     delta_cutoff = 2,
                                     grid_length = 400,
                                     chunk_size = 50000,
                                     palette = NULL,
                                     xlab = NULL,
                                     nonInt_cov = list()) {

  # Convert model_selection_table if it's from unmarkedModSel
  if (inherits(model_selection_table, "unmarkedModSel")) {
    model_selection_table <- model_selection_table@Full
  }

  if (!("delta" %in% colnames(model_selection_table)) || !("model" %in% colnames(model_selection_table))) {
    stop("model_selection_table must contain 'delta' and 'model' columns")
  }

  # Subset models by delta AIC threshold
  ms_sel <- model_selection_table[model_selection_table$delta <= delta_cutoff, ]
  if (nrow(ms_sel) == 0) stop("No models within specified delta AIC threshold.")

  model_names_sel <- ms_sel$model
  model_list_fits <- full_model_list@fits
  models_to_use_fits <- model_list_fits[model_names_sel]
  models_to_use <- methods::new("unmarkedFitList", fits = models_to_use_fits)

  # Get site covariates from first model
  umf_data <- models_to_use@fits[[1]]@data
  site_covs <- umf_data@siteCovs
  all_covs <- names(site_covs)

  if (!(cov1 %in% all_covs) || !(cov2 %in% all_covs)) {
    stop("Both covariates must be present in the site covariates.")
  }

  # Create grid of covariate combinations
  cov1_range <- range(site_covs[[cov1]], na.rm = TRUE)
  cov2_range <- range(site_covs[[cov2]], na.rm = TRUE)

  grid <- expand.grid(
    cov1_vals = seq(cov1_range[1], cov1_range[2], length.out = grid_length),
    cov2_vals = seq(cov2_range[1], cov2_range[2], length.out = grid_length)
  )
  names(grid) <- c(cov1, cov2)

  # Identify covariates not explicitly set by user
  missing_covs <- setdiff(all_covs, c(cov1, cov2, names(nonInt_cov)))
  if (length(missing_covs) > 0) {
    warning("The following covariates are not specified in nonInt_cov and will be set to their median: ",
            paste(missing_covs, collapse = ", "))
  }

  # Fill in non-interaction covariates
  for (cov in setdiff(all_covs, c(cov1, cov2))) {
    grid[[cov]] <- if (cov %in% names(nonInt_cov)) {
      nonInt_cov[[cov]]
    } else {
      median(site_covs[[cov]], na.rm = TRUE)
    }
  }

  # AIC weights
  model_weights <- exp(-0.5 * ms_sel$delta)
  model_weights <- model_weights / sum(model_weights)

  # Predict and model average
  grid$psi <- 0
  pb <- progress::progress_bar$new(
    format = "Model averaging [:bar] :percent ETA: :eta",
    total = length(models_to_use@fits),
    clear = FALSE,
    width = 60
  )

  for (i in seq_along(models_to_use@fits)) {
    mod <- models_to_use@fits[[i]]
    preds <- unmarked::predict(mod, type = "state", newdata = grid, appendData = FALSE)$Predicted
    grid$psi <- grid$psi + preds * model_weights[i]
    pb$tick()
  }

  # Plotting
  col_palette <- if (is.null(palette)) {
    hcl.colors(256, "cividis")
  } else {
    hcl.colors(256, palette)
  }

  x_label <- if (is.null(xlab)) cov1 else xlab

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data[[cov1]], y = psi,
                                          group = .data[[cov2]], colour = .data[[cov2]])) +
    ggplot2::geom_line(linewidth = 0.5, alpha = 0.8, na.rm = TRUE) +
    ggplot2::scale_colour_gradientn(colours = col_palette, name = cov2) +
    ggplot2::labs(x = x_label, y = "Probability of occupancy") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 15),
      axis.title = ggplot2::element_text(size = 17, vjust = 0.8),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 15),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    )

  return(p)
}
