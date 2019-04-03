##
## writeSparse
##

writeSparse <-
function(sparseMatrix, file="") {
  checkSparseList(sparseMatrix)
  order <- length(sparseMatrix$ia)
  nze <- length(sparseMatrix$ja)
  if (file == "")
    print(sparseMatrix)
  else if (is.character(file)) {
    out <- .C("writeSparse", as.character(file), as.integer(sparseMatrix$ia),
      as.integer(sparseMatrix$ja), as.double(sparseMatrix$a),
      as.integer(order), as.integer(nze), PACKAGE = "gep")
  }
  else
    stop("'file' must be a character string")
}

