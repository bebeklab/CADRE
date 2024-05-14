#' Classify Synergy and Antagonism
#'
#' @param data Data frame containing synergy scores.
#' @param thresholds Numeric vector of thresholds for classification.
#' @return List containing synergy list, antagonism list, and summaries.
classifySynergyAntagonism <- function(data, thresholds) {
  min_values_per_row <- apply(data[,c('synergy_bliss', 'synergy_hsa', 'synergy_zip')], 1, min, na.rm = TRUE)
  max_of_row_mins <- max(min_values_per_row)

  score_sheet <- data.frame()
  share_synergy_list <- list()
  share_antagonism_list <- list()

  for (i in 1:(length(thresholds) - 1)) {
    lower_bound <- thresholds[i]
    upper_bound <- thresholds[i + 1]
    shared_synergy <- subset(data, (synergy_zip > lower_bound & synergy_hsa > lower_bound & synergy_bliss > lower_bound) & !(synergy_zip > upper_bound & synergy_hsa > upper_bound & synergy_bliss > upper_bound))
    shared_antagonism <- subset(data, (synergy_zip < -lower_bound & synergy_hsa < -lower_bound & synergy_bliss < -lower_bound) & !(synergy_zip < -upper_bound & synergy_hsa < -upper_bound & synergy_bliss < -upper_bound))

    score_sheet <- rbind(score_sheet, data.frame(synergy_count = nrow(shared_synergy), antagonism_count = nrow(shared_antagonism)))
    share_synergy_list[[i]] <- shared_synergy
    share_antagonism_list[[i]] <- shared_antagonism
  }

  return(list(score_sheet = score_sheet, share_synergy_list = share_synergy_list, share_antagonism_list = share_antagonism_list))
}
