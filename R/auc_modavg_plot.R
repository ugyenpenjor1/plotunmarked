
#' ROC curve and AUC for model-averaged occupancy predictions
#'
#' Computes and visualises the Receiver Operating Characteristic (ROC) curve and
#' Area Under the Curve (AUC) from a set of occupancy models using model-averaged
#' predictions.
#'
#' @param model_list A list of \code{unmarkedFit} occupancy models.
#' @param use_detection Logical; if \code{TRUE}, the binary response is derived
#'   from observed detection history. If \code{FALSE}, it uses the mean posterior
#'   occupancy probability (EBUP) across models.
#' @param detection_history Optional detection/non-detection matrix (sites x surveys).
#'   Required if \code{use_detection = TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{roc}}{An object of class \code{roc} from \pkg{pROC}.}
#'   \item{\code{auc}}{AUC value from the ROC analysis.}
#'   \item{\code{auc_ci}}{Confidence interval for AUC.}
#'   \item{\code{y_true}}{Observed binary occupancy response used in ROC analysis.}
#'   \item{\code{predicted}}{Predicted model-averaged occupancy probabilities.}
#'   \item{\code{plot}}{A \code{ggplot2} ROC curve plot.}
#' }
#'
#' @details This function uses model-averaged predictions from a list of fitted
#' occupancy models (using \code{fitList}) to evaluate the model's discriminatory
#' power. It computes the AUC with 95\% confidence intervals using \pkg{pROC}, and
#' displays the ROC curve using \pkg{ggplot2}.
#'
#' @importFrom unmarked predict ranef bup fitList
#' @importFrom pROC roc
#' @importFrom ROCR prediction performance
#' @importFrom purrr map reduce
#' @importFrom ggplot2 ggplot aes geom_line geom_abline annotate labs theme_bw
#'
#' @examples
#' \dontrun{
#' # Given a list of unmarkedFit models and a detection matrix:
#' result <- auc_modavg_plot(model_list, use_detection = TRUE, detection_history = det_hist)
#' result$plot
#' result$auc
#' }
#'
#' @export

auc_modavg_plot <- function(model_list, use_detection = TRUE, detection_history = NULL) {

  # Create model-averaged fitList
  fit_list <- fitList(fits = setNames(model_list, paste0("mod", seq_along(model_list))))

  # Extract site covariates from the first model for prediction
  site_covariates <- model_list[[1]]@data@siteCovs

  # Predict model-averaged occupancy probabilities (psi)
  pred <- unmarked::predict(fit_list, newdata = site_covariates, type = "state")
  psi <- pred$Predicted  # predicted occupancy probabilities

  # Calculate mean latent occupancy from bup(ranef(.))
  bup_list <- purrr::map(model_list, ~ unmarked::bup(unmarked::ranef(.x), stat = "mean"))
  mean_b <- purrr::reduce(bup_list, `+`) / length(bup_list)

  # Decide on y_true based on user input
  if (use_detection) {
    if (is.null(detection_history)) {
      stop("If use_detection = TRUE, you must supply detection_history.")
    }
    y_true <- apply(detection_history, 1, function(x) as.integer(any(x == 1, na.rm = TRUE)))
  } else {
    y_true <- round(mean_b)
  }

  # Compute ROC and AUC
  roc_obj <- pROC::roc(
    response = y_true,
    predictor = psi,
    smoothed = TRUE,
    ci = TRUE,
    ci.alpha = 0.95,
    plot = FALSE
  )

  # Generate AUC text
  auc_ci <- as.numeric(roc_obj$ci)
  auc_text <- sprintf("AUC: %.2f\n95%% CI: %.2f-%.2f", auc_ci[2], auc_ci[1], auc_ci[3])

  # Create data for ggplot using ROCR
  pred_rocr <- ROCR::prediction(psi, y_true)
  perf <- ROCR::performance(pred_rocr, "tpr", "fpr")
  auc_df <- data.frame(fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])

  # Plot ROC
  p <- ggplot2::ggplot(auc_df, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_line(color = "blue", linewidth = 1.5) +
    ggplot2::geom_abline(linetype = "dashed") +
    ggplot2::annotate("text", x = 0.55, y = 0.3, label = auc_text, size = 4.5, hjust = 0) +
    ggplot2::labs(
      title = "",
      x = "1 - Specificity (FPR)",
      y = "Sensitivity (TPR)"
    ) +
    ggplot2::theme_bw(base_size = 14)

  print(p)

  # Return components
  return(list(
    roc = roc_obj,
    auc = roc_obj$auc,
    auc_ci = roc_obj$ci,
    y_true = y_true,
    predicted = psi,
    plot = p
  ))
}
