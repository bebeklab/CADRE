#' Generate and Analyze Graph from Drug Pathway Data with igraph package
#'
#' This function prepares a graph from drug pathway data.
#'
#' @param drug_pathway_df A data frame where rows correspond to probes and columns to different drugs.
#' @param tissue_folder A string specifying the directory to save output files.
#' @param tissue A string specifying the tissue type, used in naming the output files.
#' @return Generates a graph object and saves it as an RDS file
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join select filter group_by summarise
#' @importFrom igraph graph_from_data_frame simplify E
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_minimal
#' @export
generateGraph <- function(drug_pathway_df, tissue_folder, tissue) {
  # Add probe names as a column and reshape
  drug_pathway_df <- cbind("probes" = rownames(drug_pathway_df), drug_pathway_df)
  property_data <- reshape2::melt(drug_pathway_df, id.vars = "probes", variable.name = "Node", value.name = "Property")
  property_data <- property_data[property_data$Property != 0, ]

  property_data1 <- as.data.frame(cbind(paste0(property_data$probes, "_", property_data$Property), as.character(property_data$Node)))
  colnames(property_data1) <- c("property", "node")

  # Create adjacency matrix
  weight_drug_pathway_df <- dplyr::full_join(property_data1, property_data1, by = c('property' = 'property')) %>%
    dplyr::select(-property) %>%
    dplyr::filter(node.x != node.y) %>%
    dplyr::group_by(node.x, node.y) %>%
    dplyr::summarise(weight = dplyr::n(), .groups = 'drop')

  # Save the adjacency matrix
  saveRDS(weight_drug_pathway_df, file = paste0(tissue_folder, "RDS/", "weight_drug_pathway_", tissue, "_kegg.RDS"))

  # Form the graph and save it
  graph_pathway <- igraph::graph_from_data_frame(weight_drug_pathway_df, directed = FALSE)
  graph_pathway <- igraph::simplify(graph_pathway)
  igraph::E(graph_pathway)$weight <- igraph::E(graph_pathway)$weight / 2

  saveRDS(graph_pathway, file = paste0(tissue_folder, "RDS/graph_pathway.RDS"))
  return(graph_pathway)
}

