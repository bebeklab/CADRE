#' Calculate AUROC and AUPRC
#'
#' This function calculates the Area Under the Receiver Operating Characteristic (AUROC)
#' and the Area Under the Precision-Recall Curve (AUPRC) for a given set of drug combinations
#' and their associated probabilities. The function also outputs prediction results.
#'
#' @param tissue_folder The folder containing the RDS files.
#' @param tissue The specific tissue type for which calculations are performed.
#' @return A list containing AUROC and AUPRC values.
#' @importFrom pROC roc auc
#' @importFrom PRROC pr.curve
#' @export
calculate_performance_metrics <- function(tissue_folder, tissue, threshold=0.99) {
  probability_list <- readRDS(paste0(tissue_folder, "RDS/probability_list_all.RDS"))
  combination_tri_score <- readRDS(paste0(tissue_folder, "RDS/tri_scores.RDS"))

  prediction0 <- merge(probability_list, combination_tri_score, by = c("drug1", "drug2"))
  length1 <- nrow(prediction0)

  if (length1 < 500) {
    synergy_threshold_list <- quantile(prediction0$score, (length1 - 5) / length1) - (1-threshold)
  } else {
    synergy_threshold_list <- quantile(prediction0$score, threshold)
  }

  synergy_threshold <- synergy_threshold_list[1]
  prediction0$score <- ifelse(prediction0$score > synergy_threshold, 1, 0)
  write.csv(prediction0, paste0(tissue_folder, "prediction.csv"))

  roc_obj <- roc(prediction0$score, prediction0$probability, levels = c(0, 1), direction = "<")
  auroc_value <- auc(roc_obj)


  pr_curve <- pr.curve(scores.class0 = prediction0$score, weights.class0 = prediction0$probability, curve = TRUE)
  auprc_value <- pr_curve$auc.integral

  cat(paste0("AUROC for tissue ", tissue, " = ", auroc_value),paste0("AUPRC for tissue ", tissue, " = ", auprc_value),sep = "\n")
  metrics_output_path <- paste0(tissue_folder, tissue, "_auc.txt")
  write.table(data.frame(AUROC = auroc_value, AUPRC = auprc_value),
              file = metrics_output_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

  return(prediction0)
}
