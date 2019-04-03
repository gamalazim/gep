% File /man/mixed.Rd
% gep R package
% Distributed under GPL 2 or later

\name{mixed}
\alias{mixed}
\title{Fitting Generalized and Linear Mixed Models}
\description{
\code{mixed} is used to fit generalized and linear mixed models.
It is written for large sparse systems of equations. It handles random
effects with correlated levels. The inverse of covarince matrices can be
either built or supplied.
}
\usage{
mixed(formula, data, sol.file = "solFile", TM = FALSE, maxit = 500)
}
\arguments{
\item{formula}{an object of class \code{"\link{formula}"} (or one that
can be coerced to that class): a symbolic description of the
model to be fitted.  The details of model specification are given
under \sQuote{Details}.}

\item{data}{a data frame, list or environment (or object
coercible by \code{\link{as.data.frame}} to a data frame) containing
the variables in the model. An error message is generated if argument
\code{data} is missing or if it contains less than two variables.}

\item{sol.file}{an optional argument to specify a path and file name for writing
solutions on disk. If \code{sol.file} is missing, solutions are written to
\sQuote{./solFile}, by default.}

\item{TM}{a logical value indicating whether a threshold model analysis is to be
carried out. }

\item{maxit}{maximum number of iterations on threshold model solutions.
Default is 500.}
}

\details{
Models for \code{mixed} are specified symbolically.  A typical model has
the form \code{response ~ terms} where \code{response} is either categorical
or continuous response vector and \code{terms} is a series of terms which specifies a
linear predictor for \code{response}.  A terms specification of the form
\code{first + second} indicates all the terms in \code{first} together
with all the terms in \code{second} with duplicates removed. A
specification of the form \code{first:second} indicates the set of
terms obtained by taking the interactions of all terms in \code{first}
with all terms in \code{second}.  The specification \code{first*second}
indicates the \emph{cross} of \code{first} and \code{second}.  This is
the same as \code{first + second + first:second}.

A formula has an implied intercept term.  To remove this use either
\code{y ~ x - 1} or \code{y ~ 0 + x}. If a threshold model is fit, the code to remove
the intercept is ignored and an intercept is always fit with \code{TM = TRUE}.

A specification of the form \code{random(rterm, fvr, covMatrix=NULL, data)} indicates that \code{rterm} is
to be fit as a random factor with \sQuote{fvr = fixed variance ratio}, where
\sQuote{fvr = var(rterm)/var(e)} and \sQuote{I*fvr} is the identity matrix of \code{dim(I) = nlevels(rterm)}
multiplied by \sQuote{fvr}. Notice that a \code{data} argument is required inside \code{random} to indicate
where \code{rterm} is to be taken from. \code{covMatrix} argument may be used to specify any variance
structure besides \sQuote{I*fvr}. See function \code{\link{random}} for details.

A specification of the form \code{ranim(rterm, fvr, ped, dam=c("dam","mgs"), data)}
indicates that \code{rterm} is to be fit as a random factor with \sQuote{fvr = fixed variance ratio}, where
\sQuote{fvr = var(rterm)/var(e)} and \sQuote{G*fvr} is the inverse of the numerator relationship matrix of
\code{dim(G) = nlevels(rterm)} multiplied by \sQuote{fvr}. The \code{ped} argument inside \code{ranim} is
the numerically coded and chronologically ordered pedigree frame and \code{data}  is the data frame where
\code{rterm} is to be taken from. If the argument \code{dam} is set to \sQuote{dam}, default, the inverse
of the complete numerator relationship matrix is built, and if it is set to \sQuote{mgs}, the a sire maternal
grandsire inverse is built instead (i.e., a sire model). See function \code{\link{ranim}} for details.

A specification of the form \code{rqtl(loci, fvr, data)} fits two random allele effects per locus per individual.
A conditional gametic relationship matrix inverse is built for each quantitaive trait locus given marker data and
multiplied by one of the elements in the vector, \sQuote{fvr}, which contains the additive genetic variance of
each locus fitted. The \code{data} argument inside \code{rqtl} is a data frame with the numerically coded and
chronologically ordered pedigree in the first 3 columns and marker loci genotypes listed as two columns per locus
starting at column 4. See function \code{\link{rqtl}} for details.

}


\value{
\code{mixed} does not return an object, it only writes solutions on disk in a file name determined by the
argument \code{sol.file} as explained above. This behaviour is intended to minimize memory consumption with large
number of effects. The following functions operate on the solutions file to extract the desired results.

The function \code{summary} is used to obtain and print a summary of the analysis.
The generic accessor functions \code{coefficients}, \code{fitted.values} and
\code{residuals} extract various useful features of the solutions file generated on disk by \code{mixed}.
}

\seealso{
\code{\link{random}} for details on how to include random effects of various variance structures.
\code{\link{ranim}} for details on how to fit a genetic model.
\code{\link{rqtl}} for details on how to fit random effects of quantitative trait loci.
}
