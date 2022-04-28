#' min-max normalize the matrix, columns are cells
#' @param mat a matrix.
#' @keywords internal
#' @noRd
.normalize <- function(x){
  num <- x - min(x)
  denom <- max(x) - min(x)
  return(num/denom)
}
