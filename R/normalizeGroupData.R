#' Normalize Group Data
#'
#' Extracts and normalizes experimental conditions from the phenotypic data of an ExpressionSet.
#'
#' @param pdata A DataFrame of phenotypic data from an ExpressionSet.
#' @param cancer A character string specifying the cancer type to filter the data by.
#' @return A DataFrame with normalized group information.
#' @importFrom tidyr separate
#' @export
normalizeGroupData <- function(eset_final,pdata, cancer) {


  group_info=cbind(rownames(pdata),pdata$title,pdata$characteristics_ch1.1)
  colnames(group_info)=c("GSM","title","tissue")
  group_info=separate(data.frame(group_info),col=title,into=c("cell_line","drug","dose","time"),sep="_")
  group_info$tissue <- gsub("tissue: ", "", group_info$tissue)
  group_info$time <- gsub("^","t", group_info$time)

  group_info_sub_tissue=subset(group_info,group_info$tissue==cancer&group_info$dose=="0nM")

  title_tissue=list()
  for (cell_line in unique(group_info_sub_tissue[,2]))
  {
    for (time1 in c("2h","6h","24h"))
    {append_temp=paste(cell_line,"_topotecan_0nM_",time1,sep = "")
    title_tissue=append(title_tissue,append_temp)
    }
  }
  eset_sub_tissue=eset_final[,pdata$title %in% title_tissue]

  # Normalize name of certain celllines.
  group_info_sub_tissue$cell_line=gsub("[^a-zA-Z0-9]", "",group_info_sub_tissue$cell_line )
  group_info_sub_tissue$cell_line=gsub("LOX", "LOXIMVI",group_info_sub_tissue$cell_line)
  group_info_sub_tissue$cell_line=gsub("HL60", "HL60TB",group_info_sub_tissue$cell_line)

  return(list(group_info_sub_tissue=group_info_sub_tissue,eset_sub_tissue=eset_sub_tissue))
}
