##
## storeSparse
##
storeSparse <-
function(fullMatrix) {

  N <- dim(fullMatrix)[1]
  if(N != dim(fullMatrix)[2]) stop("Matrix must be square")
  NZE <- sum( fullMatrix[upper.tri(fullMatrix, diag=TRUE)] != 0 )
  sparse01out <-
    .C("makeSparse",
      as.double(t(fullMatrix)[lower.tri(t(fullMatrix), diag=TRUE)]),
      as.integer(N),
      as.integer(NZE),
      ia = integer(N), ja = integer(NZE), a = double(NZE), PACKAGE = "gep")

    list(ia = sparse01out$ia, ja = sparse01out$ja, a = sparse01out$a)
}

