#' Analyze Drug-Cell Interactions and calculate Tri-score
#'
#' This function calculates mean resistance indices, classifies synergy and antagonism,
#' and classifies cell line sensitivity to drug combinations, generate Tri-score.
#'
#' @param ALMANAC_sub_tissue_converted Data frame containing drug and cell line data.
#' @param tissue_folder Path to the folder for saving the output files.
#' @return NULL Invisible; this function saves several files to disk.
#' @export
analyzeDrugCellInteractions <- function(ALMANAC_sub_tissue_converted, tissue_folder,unique_cellnames,unique_drugnames) {
  # Ensure the directory exists
  dir.create(tissue_folder, recursive = TRUE, showWarnings = FALSE)
  subdirs <- c("RDS", "csv")
  lapply(subdirs, function(sub) dir.create(file.path(tissue_folder, sub), recursive = TRUE, showWarnings = FALSE))

  # Process RI Data
  ri_dataframe <- calculateMeanRI(ALMANAC_sub_tissue_converted, unique_cellnames, unique_drugnames)
  saveRDS(ri_dataframe, file.path(tissue_folder, "RDS", "ri_dataframe.RDS"))


  min_values_per_row <- apply(ALMANAC_sub_tissue_converted[,c('synergy_bliss','synergy_hsa','synergy_zip')], 1, min, na.rm = TRUE)
  # Find the minimum of these maximum values
  max_of_row_mins <- max(min_values_per_row)
  thresholds=c(floor(min(min_values_per_row)):floor(max(min_values_per_row)),Inf)

  # Process Synergy Scores
  synergy_results <- classifySynergyAntagonism(ALMANAC_sub_tissue_converted, thresholds)
  saveRDS(synergy_results, file.path(tissue_folder, "RDS", "synergy_results.RDS"))

  # Process Tri-Scores
  tri_scores <- calculateTriScores(ALMANAC_sub_tissue_converted, synergy_results$share_synergy_list, thresholds)
  saveRDS(tri_scores, file.path(tissue_folder, "RDS", "tri_scores.RDS"))

  # Sensitivity Classification
  sensitivity_classification <- classifySensitivity(ri_dataframe, tissue_folder)
  saveRDS(sensitivity_classification, file.path(tissue_folder, "RDS", "sensitivity_classification.RDS"))

  invisible(NULL)
}
