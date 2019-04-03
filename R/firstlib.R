#' @import methods
#' @import utils
#' @useDynLib gep


# The following initializes mixedMemOptions to the same values every time the package is loaded:
#options(mixedMemOptions = list(available=300000000,max_size = 30000000,AllocMem=895))
#}
#

.onLoad <- function(libname, pkgname){
     #
     # do whatever needs to be done when the package is loaded
     #

#     library.dynam("gep", pkgname, libname)
     packageStartupMessage( "Genetic Evaluation Package, gep 0.1.1, loaded" )
}
