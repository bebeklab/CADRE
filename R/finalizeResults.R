#' Convert list of interesting_pathways to data frame and save it as csv
#'
#' This function processes pathway results for different drugs and generates a CSV file summarizing the results.
#' @param pathway_results A list of data frames where each data frame corresponds to a different drug, containing columns 'pathway' and 'NES'.
#' @param pathway_gene_set A named list of pathway gene sets.
#' @param tissue_folder A string specifying the folder path for saving the output CSV file.
#' @param tissue A string specifying the tissue type, used in naming the output file.
#' @return A data frame summarizing the pathway analysis results for each drug, also saved as a CSV file in the specified location.
#' @export
finalizeResults <- function(pathway_results, pathway_gene_set, tissue_folder, tissue) {
  # Filter out empty results
  indices <- sapply(pathway_results, function(x) dim(x)[1])
  pathway_results <- pathway_results[indices != 0]

  # Prepare pathway list
  pathway_list <- data.frame(all_list = names(pathway_gene_set))

  # Process each drug
  for (drug in names(pathway_results)) {
    result_vec <- integer(length(pathway_list$all_list))
    matched_indexes <- match(pathway_list$all_list, pathway_results[[drug]]$pathway)
    non_na_indices <- which(!is.na(matched_indexes))

    result_vec[non_na_indices] <- ifelse(pathway_results[[drug]]$NES[matched_indexes[non_na_indices]] > 0, 1, -1)
    pathway_list[[drug]] <- result_vec
  }

  # Format the final data frame
  drug_pathway_df <- pathway_list[,-1]
  rownames(drug_pathway_df) <- pathway_list$all_list

  # Save to CSV
  csv_file_path <- file.path(tissue_folder, "csv", paste0("drug_pathway_", tissue, ".csv"))
  write.csv(drug_pathway_df, csv_file_path)

  return(drug_pathway_df)
}
