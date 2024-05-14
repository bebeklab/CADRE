#' Calculate Mean Resistance Index
#'
#' @param data Data frame containing drug and cell line data.
#' @param cellnames Vector of unique cell line names.
#' @param drugnames Vector of unique drug names.
#' @return Data frame with mean resistance indices.
calculateMeanRI <- function(data, cellnames, drugnames) {
  ri_dataframe <- data.frame(cellname = character(), drugname = character(), ri_mean = numeric(), stringsAsFactors = FALSE)
  for (cellname in cellnames) {
    for (drugname in drugnames) {
      ri_list <- c(data$ri_row[data$drug_row_id == drugname & data$cell_line_id == cellname],
                   data$ri_col[data$drug_col_id == drugname & data$cell_line_id == cellname])
      if (length(na.omit(ri_list)) > 0) {
        ri_mean <- mean(na.omit(ri_list))
        ri_dataframe <- rbind(ri_dataframe, data.frame(cellname = cellname, drugname = drugname, ri_mean = ri_mean))
      }
    }
  }
  return(ri_dataframe)
}
