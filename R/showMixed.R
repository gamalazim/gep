
#' Get Labels of all Solutions as Stored in Mixed Model Object
#'
#' @param mixedOBJ an object of class \sQuote{mixed}, created by \code{readCoef}.
#' @param optn default option is to get labels for \sQuote{all} variable solutions in the mixed model; other options are \sQuote{fixed} and \sQuote{random}.
#'
# Not @export -ed
getLabel <- function(mixedOBJ, optn = "all") {
  iNT <- length(mixedOBJ@threshEst)

  if(optn == "all") {
    estNames <- c(names(unlist(mixedOBJ@fixedLev)), mixedOBJ@covNames, names(unlist(mixedOBJ@randomLev)))
    estNames <- paste(estNames, c(unlist(mixedOBJ@fixedLev), mixedOBJ@covNames, unlist(mixedOBJ@randomLev)), sep='-')
    if(mixedOBJ@FI == 1 & iNT == 0)  estNames <- c("Intercept", estNames)
    if(iNT > 0) estNames <- c(paste("Thresh",1:iNT,sep=""), "Intercept", estNames)
  }

  if(optn == "fixed") {
    estNames <- c(names(unlist(mixedOBJ@fixedLev)), mixedOBJ@covNames)
    estNames <- paste(estNames, c(unlist(mixedOBJ@fixedLev), mixedOBJ@covNames), sep='-')
    if(mixedOBJ@FI == 1 & iNT == 0)  estNames <- c("Intercept", estNames)
    if(iNT > 0) estNames <- c(paste("Thresh",1:iNT,sep=""), "Intercept", estNames)
  }

  if(optn == "random") {
    estNames <- names(unlist(mixedOBJ@randomLev))
    estNames <- paste(estNames, unlist(mixedOBJ@randomLev), sep='-')
  }

  gsub(pattern=" ", replacement = "", estNames)
}


#' Summary info for a mixed object
#'
#' Interacts with an object of class \sQuote{mixed} to extract model info and random coefficients if specified.
#'
#' @param object is an object of class \sQuote{mixed}, created by \code{readCoef}.
#' @param random is logical for whether to print random solutions. Default is FALSE.
#'
#' @return None. Information is printed to standard output.
#'
#' @export
mixedSummary <- function(object, random=FALSE) {
   Est <- cbind( c(object@threshEst, object@interceptEst, unlist(object@fixedEst), object@covEst),
    c(object@threshSE, object@interceptSE, unlist(object@fixedSE), object@covSE) )

  dimnames(Est)[[2]] <- c("Est", "MSE")
  dimnames(Est)[[1]] <- getLabel(object, optn = "fixed")

  cat("\n")
  cat("Model Formula:", "\n")
  print(object@modelFormula[-1])
  cat("\n")

  info <- c(object@N, object@numEq, object@dataName)
  names(info) <- c("Number of Observations","Number of Solutions", "Data Name")
  print(info, quote=FALSE)

  cat("\n")
  print(Est)

  if(random == TRUE) {
    Est <- cbind(unlist(object@randomEst), unlist(object@randomSE), unlist(object@randomRel) )
    dimnames(Est)[[2]] <- c("Est", "MSE", "Rel")
    dimnames(Est)[[1]] <- getLabel(object, optn = "random")
    cat("\n")
    print(Est)
  }

  invisible(object)
}



#' Extract Model Fixed Coefficients
#'
#' Interacts with an object of class \sQuote{mixed} to extract all fixed factors along with any covariates fit in the model.
#'
#' @param object is an object of class \sQuote{mixed}, created by \code{readCoef}.
#' @param names is logical for whether to print solution lables (TRUE) or not (FALSE). Default is TRUE.
#'
#' @return Numeric vector with solutions labled based on model formula and data.
#'
#' @export
fixedCoef <- function(object, names=TRUE) {

  Est <- c(object@threshEst, object@interceptEst, unlist(object@fixedEst), object@covEst)
  if(names == TRUE) names(Est) <- getLabel(object, optn="fixed")
  if(names == FALSE) names(Est) <- NULL
  Est
}

#' Extract Model Random Coefficients
#'
#' Extract solutions of random factors from an object of class \sQuote{mixed}.
#'
#' @param object is an object of class \sQuote{mixed}, created by \code{readCoef}.
#' @param names is logical for whether to print solution lables (TRUE) or not (FALSE). Default is TRUE.
#' @param list logical for whether to extract the solutions in the form of a list (default) or as a numeric vector by setting \code{list = FALSE}.
#'
#' @return Numeric vector with solutions labled based on model formula and data.
#'
#' @details Used in situations where it is important to package together solutions of all random effects in one numeric vector with matching lables; this is done by setting \code{list = FALSE}. In cases where a particular random effect is needed, extract random solutions as a list using the option \code{list = TRUE}.

#' @export
randomCoef <- function(object, names=TRUE, list = TRUE) {

  if(is.null(unlist(object@randomLev))) print(NULL)

  else {
    if(list) object@randomEst
    else {
      Est <- unlist(object@randomEst)
      if(names == TRUE) names(Est) <- getLabel(object, optn = "random")
      if(names == FALSE) names(Est) <- NULL
      Est
    }
  }
}

#' Extract Reliabilities of Random Effects
#'
#' Extract reliabilities of random efects from an object of class \sQuote{mixed}.
#'
#' @param object is an object of class \sQuote{mixed}, created by \code{readCoef}.
#' @param names is logical for whether to print solution lables (TRUE) or not (FALSE). Default is TRUE.
#' @param list logical for whether to extract reliabilities in the form of a list (default) or as a numeric vector by setting \code{list = FALSE}.
#'
#' @return Numeric vector with reliabilities labled based on model formula and data.
#'
#' @details Used in situations where it is important to package together reliabilities of all random effects in one numeric vector with matching lables; this is done by setting \code{list = FALSE}. In cases where reliabilities of a particular random effect are needed, extract reliabilities as a list using the option \code{list = TRUE}.

#' @export
getRel <- function(object, names=TRUE, list = TRUE) {

  if(is.null(unlist(object@randomLev))) print(NULL)

  else {
    if(list) object@randomRel
    else {
      Rel <- unlist(object@randomRel)
      if(names == TRUE) names(Rel) <- getLabel(object, optn = "random")
      if(names == FALSE) names(Rel) <- NULL
      Rel
    }
  }
}


#' Predict Function using an Object of Class \sQuote{mixed}.
#'
#' @param mixedOBJ is an object of class \sQuote{mixed}, created by \code{readCoef}.
#' @param newdata data frame with new data to predict, default is to use training data frame associated with the mixed object
#'
#' @return Numeric vector with predictions of length equal to number of rows in newdata
#'
#' @details In the case of \code{TM = TRUE}, predictions are calculated on the uunderlying continuous scale with no transformation ot the probability scale yet.
#' @note New to gep 0.1.1 (added January 2019)
#' @export
#===============================================
# Predict Function with Object of Class 'mixed'
# January 2019
#===============================================
mixedPredict <- function(mixedOBJ, newdata = eval(parse(text=as.character(mixedOBJ@dataName)))) {

  iNT <- length(mixedOBJ@threshEst)
  df <- newdata

  df.sub <- df[,c(names(mixedOBJ@fixedLev), names(mixedOBJ@randomLev))]
  for(i in 1:ncol(df.sub)) {
    df.sub[,i] <- df.sub[,i][,drop=T]
    df.sub[,i] <- match(df.sub[,i], c(mixedOBJ@fixedLev, mixedOBJ@randomLev)[[i]])
  }

  pred <- numeric(nrow(df))
  for(i in 1:ncol(df.sub)) pred <- pred + unlist(c(mixedOBJ@fixedEst, mixedOBJ@randomEst)[i])[df.sub[,i]]
  if(mixedOBJ@FI == 1) pred <- pred + mixedOBJ@interceptEst
  k <- 0; for(i in mixedOBJ@covNames) { k <- k+1; pred <- pred + df[,i] * mixedOBJ@covEst[k] }
  pred
}






