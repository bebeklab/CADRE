#' Process Tissue Cancer Information and Convert Drug and Cell IDs
#'
#' This function processes tissue cancer information from specified rows of a data frame,
#' creates necessary directories, subsets and cleans the data, and converts drug and cell IDs.
#'
#' @param tissue_cancer_infosheet A data frame containing tissue and cancer information.
#' @param summary_df_ALMANAC A data frame summarizing ALMANAC data by tissue.
#' @param summary_df_api A data frame from API summarizing additional details.
#' @param DC_drug A data frame containing drug information including 'id' and 'stitch_name'.
#' @param DC_cell A data frame containing cell line information including 'id' and 'name'.
#' @param row_number The row number of cancer from tissue_cancer_infosheet to process.
#' @return A data frame of the processed ALMANAC tissue data with converted IDs.
#' @export
processTissueCancerInfo <- function(tissue_cancer_infosheet, summary_df_ALMANAC, summary_df_api, DC_drug, DC_cell, row_number) {
  # Subset specific row from infosheet
  tissue_cancer_infosheet_temp <- tissue_cancer_infosheet[row_number, ]
  tissue <- tissue_cancer_infosheet_temp[1, 1]
  cancer <- tissue_cancer_infosheet_temp[1, 2]

  # Create necessary folders
  tissue_folder <- paste0("./result/", tissue, "/")
  assign("tissue_folder", tissue_folder, .GlobalEnv)
  dir.create(tissue_folder, recursive = TRUE)
  dir.create(paste0(tissue_folder, "RDS"), recursive = TRUE)
  dir.create(paste0(tissue_folder, "csv"), recursive = TRUE)

  # Subset data by tissue
  ALMANAC_sub_tissue <- subset(summary_df_ALMANAC, tissue_name == tissue)
  ALMANAC_sub_tissue$cell_line_name <- gsub("[^a-zA-Z0-9]", "", toupper(ALMANAC_sub_tissue$cell_line_name))
  ALMANAC_sub_tissue2 <- subset(summary_df_api, block_id %in% ALMANAC_sub_tissue$block_id)

  # Convert drug and cell IDs
  ALMANAC_sub_tissue_converted <- ALMANAC_sub_tissue2
  drug_row_indices <- match(ALMANAC_sub_tissue_converted$drug_row_id, DC_drug$id)
  drug_col_indices <- match(ALMANAC_sub_tissue_converted$drug_col_id, DC_drug$id)
  cell_line_indices <- match(ALMANAC_sub_tissue_converted$cell_line_id, DC_cell$id)

  temp_row <- DC_drug$stitch_name[drug_row_indices]
  temp_row[is.na(temp_row)] <- ALMANAC_sub_tissue_converted$drug_row_id[is.na(temp_row)]
  temp_col <- DC_drug$stitch_name[drug_col_indices]
  temp_col[is.na(temp_col)] <- ALMANAC_sub_tissue_converted$drug_col_id[is.na(temp_col)]
  temp_cell_line <- DC_cell$name[cell_line_indices]

  ALMANAC_sub_tissue_converted$drug_row_id <- temp_row
  ALMANAC_sub_tissue_converted$drug_col_id <- temp_col
  ALMANAC_sub_tissue_converted$cell_line_id <- toupper(temp_cell_line)
  ALMANAC_sub_tissue_converted <- ALMANAC_sub_tissue_converted[ALMANAC_sub_tissue_converted$drug_row_id != ALMANAC_sub_tissue_converted$drug_col_id, ]

  return(ALMANAC_sub_tissue_converted)
}
