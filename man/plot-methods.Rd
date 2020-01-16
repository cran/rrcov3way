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
    \method{plot}{parafac}(x, which = c("dd", "comp", "percomp", "allcomp", 
    "all"), ask = (which == "all" && dev.interactive(TRUE)), id.n, ...)
}
\arguments{
  \item{x}{A \code{tucker3} object}
  \item{which}{Which plot to select (see Details). Default is \code{dd}, Distance-distance plot.}
  \item{ask}{Generates all plots in interactive mode}
  \item{id.n}{Number of items to highlight}
  \item{\dots}{Other parameters to be passed to the lower level functions.}
}
\details{
Different plots for a \code{tucker3} or \code{parafac} object will be produced. Use the parameter 
    \code{which} to select which plot to produce:
\describe{
    \item{\code{dd}}{Distance-distance plot}
    \item{\code{comp}}{Paired components plot}
    \item{\code{percomp}}{Per-component plot}
    \item{\code{allcomp}}{All components plot}
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
#############
##
## Example with the UNIDO Manufacturing value added data

data(va3way)
dim(va3way)

## Treat quickly and dirty the zeros in the data set (if any)
va3way[va3way==0] <- 0.001

##
res <- Tucker3(va3way)
res
print(res$fit)
print(res$A)

## Print the core matrix
print(res$GA)

## Distance-distance plot
plot(res, which="dd", main="Distance-distance plot")

## Paired component plot, mode A
plot(res, which="comp", main="Paired component plot (mode A)")

## Paired component plot, mode B
plot(res, which="comp", mode="B", main="Paired component plot (mode B)")

## Joint biplot
plot(res, which="jbplot", main="Joint biplot")

## Trajectory
plot(res, which="tjplot", choices=c(1:4), arrows=FALSE, main="Trajectory biplot")

}
\keyword{ Multivariate }
\keyword{ Robust }

