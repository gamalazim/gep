#' Function to read sparse-stored matrix on disc into an R list
#'
#' @details Read sparse-stored matrix on disc in the supplied file name
#'
#' @param file is a character string for the file name and directory path where the sparse-stored matrix is stored on disc.
#'
#' @return A list including 3 vectors, 'ia' for rows, 'ja' for columns, and 'a'
#'   for elements, where \code{length(sparseMatrix$ia)} = the order of the matrix and \code{length(sparseMatrix$ja)}
#'   = the number of nonzeros; and 2 scalers for the matrix order and number of non-zero elements.
#'
#' @export
readSparse <-
function(file="") {

  #if(order == 0 | nze == 0) stop("provide greater-than-0 order and nze")

  order = nze = 0
  if (is.character(file)) {
    order <- as.numeric(gsub("[^[:digit:].]", "",  read.table("LHS", nrow = 1)$V4))
    nze <- as.numeric(gsub("[^[:digit:].]", "",  read.table("LHS", nrow = 1)$V6))
    out <- .C("readSparse", as.character(file), ia=integer(order),
      ja=integer(nze), a=double(nze),
      as.integer(order), as.integer(nze), PACKAGE="gep")
  }
  else
    stop("'file' must be a character string")

  list(ia = out$ia, ja = out$ja, a = out$a, order = order, nze = nze)
}

