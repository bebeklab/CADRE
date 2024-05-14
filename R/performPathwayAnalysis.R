#' Perform Pathway Analysis
#'
#' Executes differential expression analysis using 'limma' and pathway analysis using 'fgsea'.
#' It uses parallel processing to efficiently handle computations.
#'
#' @param group_info A DataFrame containing normalized group information for constructing model matrices.
#' @param eset_sub_tissue An ExpressionSet to be used for differential expression analysis.
#' @param gene_sets A list of gene sets for performing pathway analysis.
#' @param tissue_folder A string indicating the directory where results should be saved.
#' @return A list containing the results of the pathway analysis for each drug.
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @importFrom fgsea fgsea
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply detectCores
#' @export
performPathwayAnalysis <- function(group_info, eset_sub_tissue,tissue_folder, env=environment(),sensitivity_classification,pathway_gene_set,unique_drugnames) {
  # Set up parallel processing
  cl <- makeCluster(detectCores() - 2)
  varlist <- c("group_info", "eset_sub_tissue","sensitivity_classification","pathway_gene_set","unique_drugnames")
  # Export variables from the local environment of this function
  clusterExport(cl, varlist, envir = environment())
  clusterEvalQ(cl, {
    library(Biobase)
    library(limma)
    library(fgsea)
  })
  # Perform analysis for each drug
  drug_interesting_pathway <- parLapply(cl, unique_drugnames, function(drug) {
    temp_senres=sensitivity_classification[sensitivity_classification$drugname==drug,]
    sub_group_info=cbind(group_info,senres=matrix(ncol=1,nrow=nrow(group_info)))
    sub_group_info$senres <- factor(sub_group_info$senres, levels = c("resistant", "medium", "sensitive"))
    for (cell_line_1 in unique(temp_senres$cellname))
    {sub_group_info[sub_group_info[,2]==cell_line_1,7]=
      temp_senres[temp_senres$cellname==cell_line_1,4]}

    sub_group_info$senres[is.na(sub_group_info$senres)]="medium"
    temp_design=model.matrix( ~0 + sub_group_info[['senres']])
    colnames(temp_design)=levels(as.factor(sub_group_info[['senres']]))
    contrast_matrix_temp <- makeContrasts(sensitive - resistant, levels=temp_design)

    fit_temp1 <- lmFit(eset_sub_tissue,temp_design)
    fit_temp2 <- contrasts.fit(fit_temp1,contrasts=contrast_matrix_temp)
    fit_temp2 <- eBayes(fit_temp2)

    geneRanks <- topTable(fit_temp2,number=Inf)
    geneRanks$rank = -log10(geneRanks$P.Value) * sign(geneRanks$logFC)

    # Order the genes by rank
    geneRanks=geneRanks[order(geneRanks$rank),]
    geneRanks1=geneRanks$rank
    names(geneRanks1) <- rownames(geneRanks)
    set.seed(42)
    fgseaRes_test1 <- fgsea(pathways = pathway_gene_set,
                            stats = geneRanks1,
                            minSize=15,
                            maxSize=500,nPermSimple=100000,nproc=1)

    dysregulated_pathways <- subset(fgseaRes_test1, padj < 0.05)

    return(dysregulated_pathways)
  })
  stopCluster(cl)

  names(drug_interesting_pathway)=unique_drugnames

  # Save results
  saveRDS(drug_interesting_pathway, file = paste0(tissue_folder, "RDS/pathway_analysis_results.RDS"))
  return(drug_interesting_pathway)
}
