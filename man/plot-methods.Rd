\name{plot.tucker3}
\alias{plot.tucker3}
\alias{plot.parafac}
\title{ Plot a tucker3 object
}
\description{
    Different plots for the results of Tucker3 analysis, stored in a 
    \code{tucker3} object, see Details.
}

\usage{
    \method{plot}{tucker3}(x, which = c("dd", "comp", "allcomp", "jbplot", 
    "tjplot", "all"), ask = (which == "all" && dev.interactive(TRUE)), id.n, ...)
}
\arguments{
  \item{x}{A \code{tucker3} object}
  \item{which}{Which plot to select (see Details). Default is \code{dd}, Distance-distance plot.}
  \item{ask}{Generates all plots in interactive mode}
  \item{id.n}{Number of items to highlight}
  \item{\dots}{Other parameters to be passed to the lower level functions.}
}
\details{
Different plots for a \code{tucker3} object will be produced. Use the parameter 
    \code{which} to select which plot to produce:
\describe{
    \item{\code{dd}}{Distance-distance plot}
    \item{\code{comp}}{Paired component plot}
    \item{\code{jbplot}}{Joint biplot}
    \item{\code{tjplot}}{Trajectory plot}
}
}
\references{
Kiers, H.A. (2000).Some procedures for displaying results from three-way methods. 
    \emph{Journal of Chemometrics}. 14(3): 151-170.
	
Kroonenberg, P.M. (1983).Three-mode principal component analysis: 
    Theory and applications (Vol. 2), DSWO press.
}
\author{
Valentin Todorov \email{valentin.todorov@chello.at} and 
Maria Anna Di Palma \email{madipalma@unior.it} and
Michele Gallo \email{mgallo@unior.it}
}
\examples{
##
## The Bus data from package ThreeWay
##
data(Bus)       # from ThreeWay: 7x5x37, unfolded to 7x185
dim(Bus)

## Preprocessing
##
max.scale <- c(47, 10, 10, 15, 15)
maxBus <- rep(max.scale, 37)
BusN <- t(t(Bus)/maxBus)

## Use rarray to turn back to a 3-way array
##
X <- rarray(BusN, 7, 5, 37)
dimnames(X) <- list(dimnames(Bus)[[1]], c("L-Occ", "P-Occ", "Q-Occ", 
    "S-Occ", "R-Occ"), paste("V",1:37,sep=""))

## Perform classical Tucker3 (nor robust option and no 
##  ilr transformation for compositional data)
res <- Tucker3(X, robust=FALSE, ilr=FALSE)

## Distance-distance plot
plot(res, which="dd", main="Distance-distance plot")

## Paired component plot, mode A
plot(res, which="comp", main="Paired component plot (mode A)")

## Paired component plot, mode B
plot(res, which="comp", mode="B", main="Paired component plot (mode B)") 

## Joint biplot
plot(res, which="jbplot", main="Joint biplot") 

## Trajectory plot
plot(res, which="tjplot", main="Trajectory biplot") 


}
\keyword{ Multivariate }
\keyword{ Robust }

