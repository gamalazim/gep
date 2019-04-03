##
## checkSparseList
##
checkSparseList <- function(sl) {
  if(!is.list(sl)) stop("Input sparse matrix must be a list - see 'sparse.store'", call. = FALSE)
  if(length(sl) != 5) stop("Input sparse matrix list must have 5 components", call. = FALSE)
  if(!all(names(sl) == c("ia","ja","a","order","nze"))) stop("Error in names of input sparse matrix", call. = FALSE)
}

