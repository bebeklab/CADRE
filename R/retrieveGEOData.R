#' Retrieve GEO Data
#'
#' Fetches and caches GEO datasets based on the provided GEO Series ID, ensuring data is only downloaded once.
#'
#' @param geoID A character string representing the GEO Series ID.
#' @param tissue_folder A character string indicating the directory path for saving the retrieved data.
#' @return An ExpressionSet containing the GEO data.
#' @import GEOquery
#' @export
retrieveGEOData <- function(geoID, tissue_folder) {
  geo_folder <- file.path(tissue_folder, 'geo', geoID)
  if (!dir.exists(geo_folder)) {
    dir.create(geo_folder, recursive = TRUE)
  }

  eset_path <- file.path(geo_folder, paste0('eset_', geoID, '.RDS'))
  if (!file.exists(eset_path)) {
    selected_gse <- GEOquery::getGEO(geoID, destdir = geo_folder, GSEMatrix = TRUE, getGPL = TRUE)
    eset_selected <- selected_gse[[1]]
    saveRDS(eset_selected, file = eset_path)
  } else {
    eset_selected <- readRDS(eset_path)
  }

  return(eset_selected)
}
