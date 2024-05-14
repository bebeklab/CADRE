#' Classify Sensitivity Based on Resistance Indices
#'
#' @param ri_data Data frame containing resistance indices.
#' @param tissue_folder Folder path for saving results.
#' @return Data frame with sensitivity classifications.
classifySensitivity <- function(ri_data, tissue_folder) {
  ri_drug_senres <- data.frame(matrix(ncol = ncol(ri_data) + 1, nrow = 0))
  colnames(ri_drug_senres) <- c(colnames(ri_data), "senres")

  for (drug in unique(ri_data$drugname)) {
    temp_senres <- ri_data[ri_data$drugname == drug, ]
    temp_senres$senres <- cut(temp_senres$ri_mean, breaks = quantile(temp_senres$ri_mean, probs = seq(0, 1, length.out = 4), na.rm = TRUE), include.lowest = TRUE, labels = c("resistant", "medium", "sensitive"))
    ri_drug_senres <- rbind(ri_drug_senres, temp_senres)
  }

  return(ri_drug_senres)
}
