#' Analyze Dysregulated Pathways from GEO Data
#'
#' This function retrieves and analyzes gene expression data from GEO to identify dysregulated pathways based on differential expression analysis.
#'
#' @param geoID GEO Series ID to fetch data for analysis.
#' @param tissue_folder Directory path to store intermediate and final results.
#' @param cancer Specific tissue type for analysis.
#' @param cancer Specific cancer type for analysis.
#' @param gene_sets List of gene sets for pathway analysis.
#' @return Returns a list of dysregulated pathways for further analysis.
#' @export
analyzeDysregulatedPathways <- function(geoID,sensitivity_classification,
                                        pathway_gene_set,unique_drugnames,
                                        tissue_folder,tissue,cancer) {

  eset_selected <- retrieveGEOData(geoID, tissue_folder)
  eset_final <- collapseProbesAndValidate(eset_selected)
  groupdata_list <- normalizeGroupData(eset_final,pData(eset_final), cancer)
  group_info=groupdata_list$group_info_sub_tissue
  eset_sub_tissue=groupdata_list$eset_sub_tissue

  interesting_pathways <- performPathwayAnalysis(group_info, eset_sub_tissue,tissue_folder, env=environment(),sensitivity_classification,pathway_gene_set,unique_drugnames)
  interesting_pathways=finalizeResults(interesting_pathways, pathway_gene_set, tissue_folder, tissue)

  return(interesting_pathways)
}
