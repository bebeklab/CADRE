#' Reorder Columns If Necessary
#'
#' This function reorders the first two columns of a dataframe if the first column
#' is greater than the second.
#'
#' @param x A dataframe with at least two columns.
#' @return A dataframe with possibly reordered columns.
#' @export
order_2columns <- function(x) {
  idx <- x[,1] > x[,2]
  x[idx, c(1,2)] <- x[idx, c(2,1)]
  return(x)
}
