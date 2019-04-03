##
## This file sets required classes
##


setClass("mixed",
         representation(modelFormula = "formula", dataName = "character",
                        TM = "logical", FI = "numeric", numEq = "numeric", N = "numeric",
                        threshEst = "numeric", threshSE = "numeric",
                        interceptEst = "numeric", interceptSE = "numeric",
                        fixedLev = "list", fixedEst = "list", fixedSE = "list",
                        randomLev = "list", randomEst = "list", randomSE = "list", randomRel = "list",
                        covNames = "character", covEst = "numeric", covSE = "numeric" ))


