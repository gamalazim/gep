
#' Read Mixed Model Solutions from disk and Package in an Object of Class mixed
#'
#' Function \code{mixed} saves coefficients on disk. If more work on coefficients is needed, read them back from disk and package them in an objcet of class \sQuote{mixed}.
#'
#' WARNING: Do not alter the original copy of solutions file or \code{readCoef} will fail.
#'
#' @param sol.file path and name of sotored mixed model solutions on disk. Default sol.file = "./solFile".

#' @export
readCoef <- function(sol.file = "solFile")
{
  L4 <- readLines(sol.file, n=4)
    formula <- as.formula(strsplit(L4[1], split = "Formula:")[[1]][2])
    dataName <- gsub(pattern=" ", repl="", strsplit(L4[2], split = "Data:")[[1]][2])
    TM <- gsub(pattern=" ", repl="", strsplit(L4[4], split = "TM:")[[1]][2]) != FALSE

  formula <- as.formula(formula)
  if (length(formula) < 3) stop("formula must be two-sided")

  respVar <- as.character(formula[2])
  rhsTerms <- attr(terms(formula), "term.labels")
  FI <- attr(terms(formula), "intercept")

  data <- eval(parse(text=dataName))
  N <- dim(data)[1]
  iNC <- iNT <- 0

  if(TM) {
    FI <- 1
    iNC <- length(unique(eval(parse(text=respVar), envir=data)))
    iNT <- iNC - 2
  }

## Clean up random, ranim, and rqtl
  x <- grep(pattern = "random|ranim|rqtl", rhsTerms)
  if(length(x) > 0) fixedTerms <- rhsTerms[-x]
  else fixedTerms <- rhsTerms

  randomTerms <- rhsTerms[grep(pattern = "random", rhsTerms)]
  ranimTerms <- rhsTerms[grep(pattern = "ranim", rhsTerms)]
  rqtlTerms <- rhsTerms[grep(pattern = "rqtl", rhsTerms)]
  covariateTerms <- fixedTerms[lapply(lapply(fixedTerms, function(ft) {eval(parse(text=ft), envir=data)}), is.factor)==FALSE]

  fixedFactors <- fixedTerms[!is.element(fixedTerms, covariateTerms)]
  nCov <- length(covariateTerms)  ## Number of covariate terms to fit

  optn <-
    c(unlist(lapply(randomTerms, function(rt) {eval(parse(text=paste(rt,"$optn",sep="")))})),0,0)

parametersVector <- c(
  length(fixedFactors),
  length(randomTerms),
  length(ranimTerms),
  length(rqtlTerms),
  nCov)

  if(parametersVector[1] > 0) parametersVector <- c(parametersVector,
    unlist(lapply(fixedFactors, function(ft) {nlevels(eval(parse(text=ft),envir=data)[,drop=TRUE] )}) ) )

  if(parametersVector[2] > 0) {
    parametersVector <- c(parametersVector,
      unlist(lapply(randomTerms, function(rt) {nlevels(eval(parse(text=paste(rt,"$rterm",sep="")))[,drop=TRUE])})) )
  }

  no.ind <- 0 ## just in case of no ranim factor

##
## Pedigree and MQTL stuff
##
  if(parametersVector[3] == 1 | parametersVector[4] == 1) {
    if(parametersVector[3] == 1) {
       no.ind <- unlist( {eval(parse(text=paste(ranimTerms[1],"$n",sep="")))})
       optn[(parametersVector[2]+1)] <- eval(parse(text=paste(ranimTerms[1],"$optn",sep="")))
    }
    if(parametersVector[4] == 1) {
       optn[(parametersVector[2]+2)] <- eval(parse(text=paste(rqtlTerms[1],"$optn",sep="")))
       loci <- eval(parse(text=paste(rqtlTerms[1],"$loci",sep="")))
       no.ind <- unlist( {eval(parse(text=paste(rqtlTerms[1],"$n",sep="")))})
    }
    parametersVector <- c( parametersVector, no.ind )
  }

  sumLevels <- iNT + FI + sum(parametersVector[5:(5+sum(parametersVector[1:3]))])

  out <- .C("readCoef",
    as.character(sol.file),
    as.integer(parametersVector),
    as.integer(iNC),
    as.integer(FI),

    Est=double(sumLevels),
    Err=double(sumLevels),
    Rel=double(sumLevels),
    PACKAGE="gep"
    )

    if(iNC > 2) { threshEst <- out$Est[1:iNT]; threshSE <- out$Err[1:iNT] }
    else threshEst <- threshSE <- numeric(0)

    if(FI) { interceptEst <- out$Est[1+iNT]; interceptSE <- out$Err[1+iNT] }
    else interceptEst <- interceptSE <- numeric(0)

    fixedLev <- fixedEst <- fixedSE <- list(NULL)
    if(parametersVector[1] > 0) {
      fixedLev <- lapply(fixedFactors, function(ft) {levels(eval(parse(text=ft),envir=data)[,drop=TRUE])})
      names(fixedLev) <- fixedFactors
      fixedEst <- fixedSE <- list()
      csum0 <- iNT+FI + cumsum(parametersVector[6:(5+parametersVector[1])])
      csum0 <- c(iNT+FI, csum0)
      csum1 <- csum0 + 1
      for(i in 1:parametersVector[1]) {
        fixedEst[[i]] <- out$Est[csum1[i]:csum0[i+1]]
        fixedSE[[i]] <- out$Err[csum1[i]:csum0[i+1]]
      }
    }

    randomLev <- randomEst <- randomSE <- randomRel <- list(NULL)
    if(parametersVector[2] > 0) {
      randomLev <- lapply(randomTerms, function(rt) {levels(eval(parse(text=paste(rt,"$rterm",sep="")))[,drop=TRUE])})
##      names(randomLev) <- randomTerms

      names(randomLev) <- lapply(randomTerms, function(rt) {eval(parse(text=paste(rt,"$factName",sep="")))})

      randomEst <- randomSE <- randomRel <- list()
      k <- iNT + FI + no.ind
      if(parametersVector[1] >= 1) k <- k + sum(parametersVector[6:(5+parametersVector[1])])
      csum0 <- k + cumsum(parametersVector[(6+parametersVector[1]): (5+sum(parametersVector[1:2]))])
      csum0 <- c(k, csum0)
      csum1 <- csum0 + 1

      for(i in 1:parametersVector[2]) {
        randomEst[[i]] <- out$Est[csum1[i]:csum0[i+1]]
        randomSE[[i]] <- out$Err[csum1[i]:csum0[i+1]]
        randomRel[[i]] <- out$Rel[csum1[i]:csum0[i+1]]
      }
      names(randomEst) <- names(randomSE) <- names(randomRel) <- names(randomLev)  ## randomTerms
    }

    if(parametersVector[3] > 0) {
      animPlace <- length(randomTerms) + 1
      randomLev[[animPlace]] <- lapply(ranimTerms, function(rt) {eval(parse(text=paste(rt,"$ped[,1]",sep="")))})

      k0 <- iNT + FI
      if(parametersVector[1] >= 1) k0 <- k0 + sum(parametersVector[6:(5+parametersVector[1])])
      k <- k0 + no.ind

      randomEst[[animPlace]] <- out$Est[(k0+1):k]
      randomSE[[animPlace]] <- out$Err[(k0+1):k]
      randomRel[[animPlace]] <- out$Rel[(k0+1):k]
      names(randomEst)[animPlace] <- names(randomSE)[animPlace] <- names(randomRel)[animPlace] <- names(randomLev)[animPlace] <-
        lapply(ranimTerms, function(rt) {eval(parse(text=paste(rt,"$factName",sep="")))})

    }

    covEst <- covSE <- numeric(nCov)
    if(nCov > 0) {
      for(i in 1:nCov) {
        covEst[i] <- out$Est[sumLevels-nCov+i]
	covSE[i] <- out$Err[sumLevels-nCov+i]
      }
    }

    mixedObj <- new("mixed", modelFormula = formula, dataName=dataName,
      TM = TM, FI = FI, numEq = sumLevels, N = NROW(data),
      threshEst = threshEst, threshSE = threshSE,
      interceptEst = interceptEst, interceptSE = interceptSE,
      fixedLev = fixedLev, fixedEst = fixedEst, fixedSE = fixedSE,
      randomLev = randomLev, randomEst = randomEst, randomSE = randomSE, randomRel = randomRel,
      covNames = covariateTerms, covEst = covEst, covSE = covSE)

    mixedObj
}


