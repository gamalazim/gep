
# To use random():
#    1. covMatrix is supplied as either a sparse list, or
#    2. random(rterm, fvr, storeSparse(covMatrix), data)
#' @export
random <- function(rterm, fvr, covMatrix=NULL, data) {

  mf <- match.call()
  m <- match(c("data", "rterm"), names(mf), 0L)
  mf <- mf[c(1L, m)]

  rterm <- eval(parse(text=as.character(mf$rterm)), envir = data)

## covMatrix is NULL
  if(is.null(covMatrix)) {
    optn <- 0
    covMatrix <- list(ia=0, ja=0, a=0)
  }

## covMatrix list is supplied
  else if(!is.null(covMatrix)) {
    checkSparseList(covMatrix)
    if(length(covMatrix$ia) != length(unique(rterm)))
      stop("Dimentions of supplied covMatrix do not conform with data", call. = FALSE)
    optn <- 1
  }

  list(rterm=factor(rterm)[,drop=TRUE], fvr=fvr, optn=optn, covMatrix=covMatrix, factName=as.character(mf$rterm))
}

##
## Extarct variable name from functions ranim, random, or any other function
## given that variable name comes first and no 'rterm = '
## Did not use, will not work with rterm=factName.
## Added factName item to random and ranim out lists.
##
extractName <- function(str) {
  factName <- strsplit(strsplit(str, split=",")[[1]][1], split='[:(:]')[[1]][2]
  gsub(pattern='[: \n\t:]', replacement="", factName)
}


#' @export
ranim <- function(rterm, fvr, ped, dam=c("dam","mgs"), data) {

  mf <- match.call()
  m <- match(c("ped", "rterm"), names(mf), 0L)
  mf <- mf[c(1L, m)]

  dam <- match.arg(dam)

  ped <- as.data.frame(ped)
  rterm <- eval(parse(text=as.character(mf$rterm)), envir = data)
  rterm <- ped[match(rterm, ped[,1]),1]

  optn <- 2
  if(dam == "mgs") optn <- 1

  list( rterm=rterm, fvr=fvr, optn=optn, ped=ped, n=nrow(ped), factName=as.character(mf$rterm) )
}

## Provide MQTL factors ...
rqtl <- function(loci=NULL, fvr, ped) {

  optn <- 3
  ped <- as.data.frame(ped)

  if(is.null(loci)) loci <- (ncol(ped) - 3)/2
  else  if(!is.numeric(loci)) stop("Argument 'loci' must be numeric, if supplied", call. = FALSE)

  if( length(fvr) != loci ) stop("Number of loci != length(fvr)", call. = FALSE)

  list( loci=loci, fvr=fvr, optn=optn, gped=ped, n=nrow(ped) )
}


######
## R function to wrap mixed.c. mixed.c computes genetic evaluations.
######

##
## 'mixed'
##
#' @export
mixed <-
  function(formula, data, sol.file = "solFile", TM = FALSE, maxit = 100, doNotOrder = FALSE) {

  formula <- as.formula(formula)
  if (length(formula) < 3) stop("formula must be two-sided", call. = FALSE)

## dummy function for formula evaluation
## Conflicts with sequences
##":" <- "%in%" <- function(a, b) { a <- factor(a); b <- factor(b); interaction(a, b, drop=TRUE)}

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- c(as.character(mf[m[1]]), as.character(mf[m[2]]))
  if(TM) mf <- c(mf, "TRUE") else mf <- c(mf, "FALSE")

  respVar <- as.character(formula[2])
  rhsTerms <- attr(terms(formula), "term.labels")
  FI <- attr(terms(formula), "intercept")

  if(dim(data)[2] < 2)
    stop("data must contain at least 2 columns", call. = FALSE)

  N <- dim(data)[1]
  iNC <- 0

  if(TM) {
    if(!FI) warning("an intercept is always fit with threshold model")
    FI <- 1
    iNC <- length(unique(eval(parse(text=respVar), envir=data)))
    if(min(unique(eval(parse(text=respVar), envir=data))) != 1)
      stop("Categories must be a sequence from 1 to k with increment of 1", call. = FALSE)
  }

## Clean up random, ranim, and rqtl
  x <- grep(pattern = "random|ranim|rqtl", rhsTerms)
  if(length(x) > 0) fixedTerms <- rhsTerms[-x]
  else fixedTerms <- rhsTerms

#fixedTerms <- paste("data", fixedTerms, sep="$")

  randomTerms <- rhsTerms[grep(pattern = "random", rhsTerms)]
  ranimTerms <- rhsTerms[grep(pattern = "ranim", rhsTerms)]
  rqtlTerms <- rhsTerms[grep(pattern = "rqtl", rhsTerms)]
  covariateTerms <- fixedTerms[lapply(lapply(fixedTerms, function(ft) {eval(parse(text=ft), envir=data)}), is.factor)==FALSE]

  fixedFactors <- fixedTerms[!is.element(fixedTerms, covariateTerms)]
  nCov <- length(covariateTerms)  ## Number of covariate terms to fit

  optn <-
    c(unlist(lapply(randomTerms, function(rt) {eval(parse(text=paste(rt,"$optn",sep="")))})),0,0)

##
##
##
parametersVector <- c(
  length(fixedFactors),
  length(randomTerms),
  length(ranimTerms),
  length(rqtlTerms),
  nCov)


## TM = TRUE but no random factors is considered error, one reason is, TM iterates on random
## effects until it reaches a desirable stoping criterion
##
  if(TM & sum(parametersVector[2:4]) == 0)
      stop("To fit a threshold model, you need to provide at least one random factor.", call. = FALSE)


  if(parametersVector[1] > 0) parametersVector <- c(parametersVector,
    unlist(lapply(fixedFactors, function(ft) {nlevels(eval(parse(text=ft),envir=data)[,drop=TRUE] )}) ) )

  random.fvr <- 0 ## just in case of no random factors in the model
  if(parametersVector[2] > 0) {
    parametersVector <- c(parametersVector,
      unlist(lapply(randomTerms, function(rt) {nlevels(eval(parse(text=paste(rt,"$rterm",sep="")))[,drop=TRUE])})) )
    random.fvr <- lapply(randomTerms, function(rt) {eval(parse(text=paste(rt,"$fvr",sep="")))})
  }

  no.ind <- 0 ## just in case of no ranim factor

##
## Pedigree and MQTL stuff
##
  if(parametersVector[3] > 1 | parametersVector[4] > 1) stop("repeated 'ranim()/rqtl()' not allowed", call. = FALSE)
  if(parametersVector[3] == 1 | parametersVector[4] == 1) {
    sr <- c(0, unlist( {eval(parse(text=paste(ranimTerms[1],"$ped[,2]",sep="")))}) )
    dm <- c(0, unlist( {eval(parse(text=paste(ranimTerms[1],"$ped[,3]",sep="")))}) )

    if(parametersVector[3] == 1) {
       no.ind <- unlist( {eval(parse(text=paste(ranimTerms[1],"$n",sep="")))})
       random.fvr <- c( eval(parse(text=paste(ranimTerms[1],"$fvr",sep=""))), random.fvr )
       optn[(parametersVector[2]+1)] <- eval(parse(text=paste(ranimTerms[1],"$optn",sep="")))
    }

    if(parametersVector[4] == 1) {
       optn[(parametersVector[2]+2)] <- eval(parse(text=paste(rqtlTerms[1],"$optn",sep="")))
       loci <- eval(parse(text=paste(rqtlTerms[1],"$loci",sep="")))
       no.ind <- unlist( {eval(parse(text=paste(rqtlTerms[1],"$n",sep="")))})
    }
    parametersVector <- c( parametersVector, no.ind )
  }
  if(parametersVector[3] == 0 & parametersVector[4] == 0) sr <- dm <- 0

## covMatrices, associated only with random() not with ranim() nor with rqtl().
  iaj <-
    lapply(randomTerms, function(rt) {eval(parse(text=paste(rt,"$covMatrix",sep="")))})
  iaj <- unlist(iaj)

## Although QTLs could be fitted w/t a polygenic component, I allowed fitting a QTL
## only when a polygenic component is fitted with an animal not a sire model.
##
  if(parametersVector[4] == 1 & optn[parametersVector[2]+1] != 2)
    stop("You need to fit a polygenic component by 'ranim' with complete pedigree", call. = FALSE)

  orderOpt = 0 ## ordering
  if(doNotOrder) orderOpt = 1 ## don't order

##
## if mixedMemOptions list has been changed/set use the new values
## mixedMemOptions can be set by
##   options(mixedMemOptions = list(
##         available = 300000000,
##         max_size  = 30000000,
##         AllocMem  = 895 ))
## Setting any of the above options to a 0 will invoke usage of its default value instead.
##
  memOptions <- c(0,0,0)
  if(!is.null(getOption("mixedMemOptions"))) memOptions <- unlist(getOption("mixedMemOptions"))
  ############################################
  ## if neither ranim() nor rqtl() was used ##
  ############################################
    if(parametersVector[3] == 0 & parametersVector[4] == 0) {
    mixedOut <-
      .C("mixed", as.character(c(fixedFactors,ranimTerms,randomTerms,rqtlTerms,covariateTerms)), as.character(mf),
      as.double(unlist(c( unclass(eval(parse(text=respVar), envir=data)[,drop=TRUE]),
		lapply(fixedFactors, function(ft) {unclass(eval(parse(text=ft),envir=data)[,drop=TRUE])}),
		lapply(randomTerms, function(rt) {unclass(eval(parse(text=paste(rt,"$rterm",sep="")))[,drop=TRUE])}),
		lapply(covariateTerms, function(rt) {unclass(eval(parse(text=rt),envir=data)[,drop=TRUE])})))),

      as.integer(iNC),
      as.integer(N),
      as.integer(FI),
      as.integer(parametersVector),
      as.integer(optn),
      as.integer(sr),
      as.integer(dm),
      as.double(iaj),
      as.integer(0),
      as.double(0),
      as.integer(0),
      as.double(random.fvr),
      as.character(sol.file),
      as.integer(maxit),
      as.integer(orderOpt),
      as.integer(memOptions),
      PACKAGE="gep"
      )
    }

  #############################################
  ## if only ranim() but not rqtl() was used ##
  #############################################
    if(parametersVector[3] == 1 & parametersVector[4] == 0) {
    mixedOut <-
      .C("mixed", as.character(c(fixedFactors,ranimTerms,randomTerms,rqtlTerms,covariateTerms)), as.character(mf),
      as.double(unlist(c( unclass(eval(parse(text=respVar), envir=data)[,drop=TRUE]),
		lapply(fixedFactors, function(ft) {unclass(eval(parse(text=ft),envir=data)[,drop=TRUE])}),
		lapply(ranimTerms, function(rt) {unclass(eval(parse(text=paste(rt,"$rterm",sep="")))[,drop=TRUE])}),
		lapply(randomTerms, function(rt) {unclass(eval(parse(text=paste(rt,"$rterm",sep="")))[,drop=TRUE])}),
		lapply(covariateTerms, function(rt) {unclass(eval(parse(text=rt),envir=data)[,drop=TRUE])})))),

      as.integer(iNC),
      as.integer(N),
      as.integer(FI),
      as.integer(parametersVector),
      as.integer(optn),
      as.integer(sr),
      as.integer(dm),
      as.double(iaj),
      as.integer(0),
      as.double(0),
      as.integer(0),
      as.double(random.fvr),
      as.character(sol.file),
      as.integer(maxit),
      as.integer(orderOpt),
      as.integer(memOptions),
		  PACKAGE="gep"
      )
    }

  ###############################################
  ## if rqtl() were used, regardles of ranim() ##
  ###############################################
    if( parametersVector[4] == 1 ) {
    mixedOut <-
      .C("mixed", as.character(c(fixedFactors,ranimTerms,randomTerms,rqtlTerms,covariateTerms)), as.character(mf),
      as.double(unlist(c( unclass(eval(parse(text=respVar), envir=data)[,drop=TRUE]),
		lapply(fixedFactors, function(ft) {unclass(eval(parse(text=ft),envir=data)[,drop=TRUE])}),
		lapply(ranimTerms, function(rt) {unclass(eval(parse(text=paste(rt,"$rterm",sep="")))[,drop=TRUE])}),
		lapply(randomTerms, function(rt) {unclass(eval(parse(text=paste(rt,"$rterm",sep="")))[,drop=TRUE])}),
		lapply(covariateTerms, function(rt) {unclass(eval(parse(text=rt),envir=data)[,drop=TRUE])})))),

      as.integer(iNC),
      as.integer(N),
      as.integer(FI),
      as.integer(parametersVector),
      as.integer(optn),
      as.integer(sr),
      as.integer(dm),
      as.double(iaj),
      as.integer(loci),
      as.double(eval(parse(text=paste(rqtlTerms[1],"$fvr",sep=""))) ),
      as.integer(unlist(eval(parse(text=paste(rqtlTerms[1],"$gped",sep="")))[,4:(3+2*loci)])),
      as.double(random.fvr),
      as.character(sol.file),
      as.integer(maxit),
      as.integer(orderOpt),
      as.integer(memOptions),
		  PACKAGE="gep"
      )
  }
}














