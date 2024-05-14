#' Predict drug synergy from cluster
#'
#' This function calculates synergy scores for drug combinations
#' based on clustering of in graph.
#' It outputs a list of probabilities indicating potential synergies.
#'
#' @param tissue_folder The folder containing the RDS files.
#' @param threshold_graph The threshold to use for the initial filter on graph edges.
#' @return A list of drug combinations with their synergy probabilities.
#' @importFrom dplyr %>% ntile mutate rowwise
#' @importFrom igraph V E cluster_walktrap membership delete.edges delete.vertices degree get.edgelist
#' @export
prediction_by_cluster <- function(tissue_folder, threshold_graph = 95,tri_scores,graph_pathway) {
  drug_pair_sheet=readRDS(paste0(tissue_folder,"RDS/tri_scores.RDS"))[,1:2]
  drug_pair_sheet$appear=0
  drug_pair_sheet$synergy=0
  threshold_list=seq(0, 99, by=1)
  for (threshold_edge_vector in (threshold_graph+1):100){
    threshold_edge=threshold_list[threshold_edge_vector]

    if(length(E(graph_pathway))>400){
      threshold1=max(E(graph_pathway)$weight[ntile(E(graph_pathway)$weight, 100)==threshold_edge])
    }else{
      threshold1=as.numeric(quantile(E(graph_pathway)$weight,((length(E(graph_pathway))-(100-threshold_edge)*4)/length(E(graph_pathway)))))
    }
    graph_pathway_trim=delete.edges(graph_pathway, E(graph_pathway)[E(graph_pathway)$weight <= threshold1])
    graph_pathway_trim=delete.vertices(graph_pathway_trim, which(degree(graph_pathway_trim)<=0))

    set.seed(42)
    clusters=cluster_walktrap(graph_pathway_trim)
    membership=as.data.frame(membership(clusters))
    cluster_results=data.frame(id=rownames(membership),cluster=membership[,1])

    nodes_index=data.frame(cbind(c(1:length(V(graph_pathway_trim)$name)[1]),V(graph_pathway_trim)$name))
    # Filter and replace elements in tri_scores
    combination_tri_score_filtered = subset(
      tri_scores,
      tri_scores$drug1 %in% nodes_index$X2 &
        tri_scores$drug2 %in% nodes_index$X2
    )
    cluster_names=names(table(cluster_results[,2]))
    ## loop within clusters to get the combination within same cluster.
    cluster_prediction=data.frame()
    # Fetch the edge list
    edge_list <- get.edgelist(graph_pathway_trim)
    # Fetch the weights
    edge_weights <- E(graph_pathway_trim)$weight
    # Create a data frame from the edge list and weights
    edge_info_df <- data.frame(
      from = edge_list[, 1],
      to = edge_list[, 2],
      weight = edge_weights
    )
    edge_info_df=order_2columns(edge_info_df)
    for (i in 1:(length(cluster_names)-1)) {
      cluster_i = unlist(subset(cluster_results[,1], cluster_results[,2] == cluster_names[i]))
      for (j in (i + 1):length(cluster_names)) {
        cluster_j = unlist(subset(cluster_results[,1], cluster_results[,2] == cluster_names[j]))

        combinations = as.data.frame(expand.grid(cluster_i, cluster_j))
        cluster_prediction = rbind(cluster_prediction, combinations)}}

    cluster_prediction <- as.data.frame(lapply(cluster_prediction, function(x) {
      if (is.factor(x)) as.character(x) else x
    }))

    cluster_prediction=order_2columns(cluster_prediction)
    colnames(cluster_prediction)=c("drug1","drug2")

    drug_pair_sheet <- drug_pair_sheet %>%
      rowwise() %>%
      mutate(
        appear = if(any((combination_tri_score_filtered$drug1 == drug1) & (combination_tri_score_filtered$drug2 == drug2))) {
          appear + 1
        } else {
          appear
        }
      )
    drug_pair_sheet <- drug_pair_sheet %>%
      rowwise() %>%
      mutate(
        synergy = if(any((cluster_prediction$drug1 == drug1) & (cluster_prediction$drug2 == drug2))) {
          synergy + 1
        } else {
          synergy
        }
      )
  }


  drug_pair_sheet=drug_pair_sheet[drug_pair_sheet$appear!=0,]
  drug_pair_sheet$probability=drug_pair_sheet$synergy/(100-threshold_graph)
  saveRDS(drug_pair_sheet, paste0(tissue_folder, "RDS/probability_list_all.RDS"))

  return(drug_pair_sheet)
}
