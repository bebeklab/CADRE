#' Load and Preprocess Drug and Cell Line Data
#'
#' This function takes paths to drug and cell line data files and optionally other data files,
#' loads them, applies necessary transformations, and assigns each dataframe to a variable
#' in the specified environment.
#'
#' @param drug_data_path Path to the drug data CSV file from DrugComb.
#' @param cell_data_path Path to the cell line data CSV file from DrugComb.
#' @param summary_df_ALMANAC_path Path to the summary data filtered by ALMANAC, CSV file from DrugComb.
#' @param summary_df_api_path Path to the summary data download by api, CSV file from DrugComb.
#' @param pathway_gene_set_path Path to the pathway gene set file, download from https://www.gsea-msigdb.org/.
#' @param tissue_cancer_infosheet_path Path to the cancer-tissue info sheet of GSE entry, CSV file.
#' @param env Environment in which to assign the dataframes (default is the global environment).
#' @importFrom dplyr mutate %>%
#' @importFrom fgsea gmtPathways
#' @return NULL Invisible; this function release variables to .GlobalEnv.
#' @export
load_and_preprocess_data <- function(drug_data_path=NULL,
                                     cell_data_path=NULL,
                                     summary_df_ALMANAC_path=NULL,
                                     summary_df_api_path=NULL,
                                     pathway_gene_set_path=NULL,
                                     tissue_cancer_infosheet_path=NULL,
                                     env = .GlobalEnv) {
  # Load and preprocess drug data
  if(!is.null(drug_data_path)){
    DC_drug <- read.csv(drug_data_path)
    DC_drug <- dplyr::mutate(DC_drug, stitch_name = ifelse(stitch_name == "" | is.na(stitch_name), dname, stitch_name))
    assign("DC_drug", DC_drug, env)}



  # Load and preprocess cell line data
  if(!is.null(cell_data_path)){
    DC_cell <- read.csv(cell_data_path)
    DC_cell$name <- gsub("[^a-zA-Z0-9]", "", toupper(DC_cell$name))
    assign("DC_cell", DC_cell, env)}

  # Process optional other data files

  if(!is.null(summary_df_ALMANAC_path)){
    summary_df_ALMANAC <- read.csv(summary_df_ALMANAC_path)
    assign("summary_df_ALMANAC", summary_df_ALMANAC, env)}

  if(!is.null(summary_df_api_path)){
    summary_df_api <- read.csv(summary_df_api_path)
    assign("summary_df_api", summary_df_api, env)}

  if(!is.null(pathway_gene_set_path)){
    pathway_gene_set <- gmtPathways(pathway_gene_set_path)
    assign("pathway_gene_set", pathway_gene_set, env)}

  if(!is.null(tissue_cancer_infosheet_path)){
    tissue_cancer_infosheet <- read.csv(tissue_cancer_infosheet_path)
    assign("tissue_cancer_infosheet", tissue_cancer_infosheet, env)}

  invisible(NULL)
}
