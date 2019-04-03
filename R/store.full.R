
#' Transform sparse-stored to full-stored matrix
#'
#' @details The input sparse matrix is a list of 3 vectors, 'ia' for rows, 'ja' for columns, and 'a'
#'   for elements, where \code{length(sparseMatrix$ia)} = the order of the matrix and \code{length(sparseMatrix$ja)}
#'   = the number of nonzeros.
#'
#' @param sparseMatrix is a list created by function \code{readSparse}.
#'
#' @return Full-stored matrix with dimensions = order of input sparse matrix.
#'
#' @export
storeFull <-
function(sparseMatrix) {
  checkSparseList(sparseMatrix)
  N <- sparseMatrix$order
  fullMatrix <- matrix(0, N, N)

  for(i in 1:(N-1)) {
    cols <- sparseMatrix$ja[sparseMatrix$ia[i]:(sparseMatrix$ia[i+1]-1)]
    elems <- sparseMatrix$a[sparseMatrix$ia[i]:(sparseMatrix$ia[i+1]-1)]
    fullMatrix[i,cols] <- elems
    fullMatrix[cols,i] <- elems
  }
  fullMatrix[N,N] <- sparseMatrix$a[sparseMatrix$nze]

  fullMatrix
}

