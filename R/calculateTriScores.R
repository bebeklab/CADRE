#' Calculate Tri-Scores for Drug Combinations
#'
#' @param dataset dataset frame with drug interaction dataset.
#' @param synergy_list List of dataset frames containing synergy summaries for different thresholds.
#' @param thresholds Numeric vector of thresholds for scoring.
#' @return dataset frame with tri-scores for each drug combination.
calculateTriScores <- function(dataset, synergy_list, thresholds) {
  combination_tri_score = list()


  unique_drugs <- unique(c(dataset$drug_row_id, dataset$drug_col_id))

  for (i in 1:(length(unique_drugs) - 1)) {
    for (j in (i + 1):length(unique_drugs)) {

      num_cellline=dim(subset(dataset,
                              (dataset$drug_row_id== unique_drugs[i] & dataset$drug_col_id== unique_drugs[j])|
                                (dataset$drug_col_id== unique_drugs[i] & dataset$drug_row_id== unique_drugs[j]
                                )))[1]
      if(num_cellline!=0){
        temp_sum=NULL
        for(thresholdi in 1:(length(thresholds) - 1)){
          share_synergy_list_sub=synergy_list[[thresholdi]]
          temp <- subset(share_synergy_list_sub,
                         ((share_synergy_list_sub$drug_row_id == unique_drugs[i]) &
                            (share_synergy_list_sub$drug_col_id == unique_drugs[j])) |
                           ((share_synergy_list_sub$drug_col_id == unique_drugs[i]) &
                              (share_synergy_list_sub$drug_row_id == unique_drugs[j]))  )
          temp_sum=rbind(temp_sum,nrow(temp)*(thresholdi+min(thresholds)))
        }
        combination_tri_score[[length(combination_tri_score) + 1]] <- data.frame(
          drug1 = unique_drugs[i],
          drug2 = unique_drugs[j],
          score = sum(temp_sum) / (num_cellline / length(unique(dataset$cell_line_id))),
          stringsAsFactors = FALSE)
      }

    }
  }


  # Combine all dataset frames in the list into a single dataset frame
  combination_tri_score <- do.call(rbind, combination_tri_score)


  combination_tri_score <- order_2columns(combination_tri_score)
    return(combination_tri_score)
}
