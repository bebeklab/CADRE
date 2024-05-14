#' Collapse probes to genes
#'
#' This function summarize several probe measurements corresponding to a single gene.
#'
#' @param x ExpressionSet needs to reduced probes / genes.
#' @return Reduced ExpressionSet.
#' @export
collapseProbes <- function (x, method = "variance", topTable = NULL, colNameOfStat = NULL,
    colNameOfGeneSymbol = "SYMBOL")
{
    if (!(method %in% list("variance", "logFC", "P.Value", "adj.P.Value"))) {
        stop("Method ", paste(method, "not supported"))
    }
    if (is.data.frame(topTable)) {
        if (!(method %in% list("variance", "logFC", "P.Value",
            "adj.P.Value"))) {
            stop("Method ", paste(method, "not supported"))
        }
        if (!(colNameOfGeneSymbol %in% colnames(topTable))) {
            stop(paste(colNameOfGeneSymbol, "must be a column in topTable"))
        }
    }
    if (!is.matrix(x)) {
        mat_exprs <- exprs(x)
    }
    if (method == "variance") {
        cat("Processing eset using:", method, "\n")
        vec_variance <- apply(mat_exprs, 1, var)
        vec_symbol <- fData(x)[rownames(mat_exprs), colNameOfGeneSymbol]
        df_xref <- as.data.frame(cbind(mat_exprs, variance = vec_variance))
        df_xref <- cbind(df_xref, `Gene Symbol` = vec_symbol,
            probeID = rownames(mat_exprs))
        df_xref <- df_xref[order(df_xref[[method]], decreasing = T),
            ]
        uniqueGeneList <- unique(df_xref[["Gene Symbol"]])
        new_matrix <- df_xref[match(uniqueGeneList, table = df_xref[["Gene Symbol"]]),
            ]
        x <- x[as.vector(new_matrix$probeID), ]
        esetToReturn <- x
    }
    else if (method == "logFC") {
        cat("Processing eset using:", method, "\n")
        topTable <- topTable[order(topTable[[method]], decreasing = TRUE),
            ]
        df_xref <- df_xref[rownames(topTable), ]
        if (!assertthat::are_equal(rownames(df_xref), rownames(topTable))) {
            stop("Matrix rows do not match topTable rows in method=logFC")
        }
        df_xref <- cbind(df_xref, topTable)
        df_xref[[colNameOfGeneSymbol]] <- df_xref[[colNameOfGeneSymbol]]
        df_xref <- df_xref[order(df_xref[[method]], decreasing = T),
            ]
        uniqueGeneList <- unique(df_xref[[colNameOfGeneSymbol]])
        new_matrix <- df_xref[match(uniqueGeneList, table = df_xref[[colNameOfGeneSymbol]]),
            ]
        mat_exprs <- new_matrix[, colnames(mat_exprs)]
        rownames(mat_exprs) <- as.vector(new_matrix[[colNameOfGeneSymbol]])
        esetToReturn <- x[as.vector(new_matrix$probeID), ]
        exprs(esetToReturn) <- as.matrix(mat_exprs)
        new_fData <- new_matrix[, c("Gene Symbol", "variance",
            "probeID")]
        rownames(new_fData) <- new_fData[[colNameOfGeneSymbol]]
        fData(esetToReturn) <- new_fData
        assertthat::are_equal(rownames(exprs(esetToReturn)),
            rownames(fData(esetToReturn)))
    }
    return(esetToReturn)
}
