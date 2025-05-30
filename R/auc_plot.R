
#' ROC curve and AUC for single occupancy model
#'
#' Computes and visualises the Receiver Operating Characteristic (ROC) curve
#' and Area Under the Curve (AUC) for a single fitted occupancy model from
#' the \pkg{unmarked} package.
#'
#' @param model An \code{unmarkedFit} object from the \pkg{unmarked} package.
#' @param use_detection Logical; if \code{TRUE}, the binary presence/absence
#'   response is based on observed detection history. If \code{FALSE}, it is
#'   based on the model's posterior mean occupancy (EBUP).
#' @param detection_history Optional matrix of detection/non-detection data
#'   (1 = detection, 0 = no detection, NA = missing). Required if
#'   \code{use_detection = TRUE}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{roc}}{An object of class \code{roc} from \pkg{pROC} containing ROC data.}
#'   \item{\code{auc}}{The AUC value.}
#'   \item{\code{auc_ci}}{The 95\% confidence interval for AUC.}
#'   \item{\code{y_true}}{Binary observed occupancy (used as response).}
#'   \item{\code{predicted}}{Predicted occupancy probabilities.}
#'   \item{\code{plot}}{A \code{ggplot2} ROC curve plot.}
#' }
#'
#' @details The function calculates predicted occupancy probabilities from
#'   \code{model}, compares them to either observed presence/absence
#'   (via \code{detection_history}) or to the model's posterior mean estimates,
#'   and computes AUC statistics using \pkg{pROC}. The ROC curve is then
#'   visualised using \pkg{ggplot2}.
#'
#' @importFrom pROC roc
#' @importFrom unmarked ranef bup predict
#' @importFrom ROCR prediction performance
#' @importFrom ggplot2 ggplot aes geom_line geom_abline annotate labs theme_bw
#'
#' @examples
#' \dontrun{
#' # Assuming 'fitted_model' is an unmarked occupancy model
#' # and 'y_mat' is the detection matrix:
#' result <- auc_single(fitted_model, use_detection = TRUE, detection_history = y_mat)
#' print(result$plot)
#' result$auc
#' }
#'
#' @export

auc_plot <- function(model, use_detection = TRUE, detection_history = NULL) {

  # Extract site-level covariates from model
  site_covariates <- model@data@siteCovs

  # Predict occupancy probabilities
  occ_pred <- unmarked::predict(model, newdata = site_covariates, type = "state", inf.rm = TRUE)
  psi <- occ_pred$Predicted

  # Get latent occupancy estimate (posterior mean from ranef)
  EBUP <- unmarked::bup(unmarked::ranef(model), stat = "mean")

  # Decide on "truth" (y_true)
  if (use_detection) {
    if (is.null(detection_history)) {
      stop("You must provide detection_history if use_detection = TRUE")
    }
    y_true <- apply(detection_history, 1, function(x) as.integer(any(x == 1, na.rm = TRUE)))
  } else {
    y_true <- round(EBUP)
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

  # Format AUC text
  auc_ci <- as.numeric(roc_obj$ci)
  auc_text <- sprintf("AUC: %.2f\n95%% CI: %.2f-%.2f", auc_ci[2], auc_ci[1], auc_ci[3])

  # Build ROC data using ROCR for plotting
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

  # Return AUC components
  return(
    invisible(
      list(
        roc = roc_obj,
        auc = roc_obj$auc,
        auc_ci = roc_obj$ci,
        y_true = y_true,
        predicted = psi,
        plot = p
      )
    )
  )
}

