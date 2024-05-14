#' Collapse Probes and Validate Data
#'
#' Collapses probes in an ExpressionSet using a specified method and ensures gene symbols are correctly assigned.
#'
#' @param eset An ExpressionSet object to be processed.
#' @param method A character string indicating the collapsing method, defaults to 'variance'.
#' @return A modified ExpressionSet with collapsed probes.
#' @import Biobase
#' @export
collapseProbesAndValidate <- function(eset, method = "variance") {
  eset_collapsed <- collapseProbes(eset, colNameOfGeneSymbol = "Gene Symbol", method = method)
  fData <- fData(eset_collapsed)
  eset_collapsed <- eset_collapsed[!is.na(fData$`Gene Symbol`),]
  fData <- fData(eset_collapsed)
  pData <- pData(eset_collapsed)
  mat <- exprs(eset_collapsed)

  if(!(assertthat::are_equal(rownames(fData), rownames(mat)) &
       assertthat::are_equal(rownames(pData), colnames(mat))) ) {
    stop("fData or pData feature names does not match\n")
  } else {
    cat("fData and pData match exprs matrix.")
  }

  rownames(mat) <- fData(eset_collapsed)$`Gene Symbol`
  rownames(fData) <- fData(eset_collapsed)$`Gene Symbol`
  eset_collapsed.symbol  <- ExpressionSet(assayData   = mat,
                                          phenoData   = new("AnnotatedDataFrame",data = pData  ),
                                          featureData = new("AnnotatedDataFrame", data = fData),
                                          annotation  = "illuminaHumanv4")

  eset_collapsed.symbol <- eset_collapsed.symbol[!grepl(pattern = "LOC.*"
                                                        , x = rownames(exprs(eset_collapsed.symbol))), ]
  eset_collapsed.symbol <- eset_collapsed.symbol[!grepl(pattern = "KIA.*"
                                                        , x = rownames(exprs(eset_collapsed.symbol))), ]
  eset_collapsed.symbol <- eset_collapsed.symbol[!grepl(pattern = ".*orf.*"
                                                        , x = rownames(exprs(eset_collapsed.symbol))), ]

  eset_collapsed=eset_collapsed.symbol
  return(eset_collapsed)
}
