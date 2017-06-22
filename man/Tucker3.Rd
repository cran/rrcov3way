\name{Tucker3}
\alias{Tucker3}
\alias{print.tucker3}
\title{
    Robust Tucker3 estimator for compositional data
}
\description{
Compute a robust Tucker3 model for compositional data
}
\usage{
Tucker3(X, P = 2, Q = 2, R = 2, conv = 1e-06, 
    center = FALSE, scale = FALSE, scale.mode = c("B", "A", "C"), 
    robust = FALSE, ilr = FALSE, ncomp.rpca = 2, alpha = 0.75, 
    maxiter=100, crit=0.975, trace = FALSE)
}
\arguments{
  \item{X}{3-way array of data}
  \item{P}{Number of A-mode components}
  \item{Q}{Number of B-mode components}
  \item{R}{Number of C-mode components}
  \item{conv}{Convergence criterion, defaults to \code{1e-6}}
  \item{center}{Whether to center the data}
  \item{scale}{Whether to scale the data}
  \item{scale.mode}{If centering and/or scaling the data, on which mode to do this}
  \item{robust}{Whether to apply a robust estimation}
  \item{ilr}{If the data are a composition, use an ilr transformation}
  \item{ncomp.rpca}{Number of components for robust PCA}
  \item{alpha}{Measures the fraction of outliers the algorithm should
    resist. Allowed values are between 0.5 and 1 and the default is 0.75}
  \item{maxiter}{Maximal number of iterations}
  \item{crit}{Cut-off for identifying outliers, default \code{crit=0.975}}
  \item{trace}{Logical, provide trace output}
}
\details{
The function can compute four
    versions of the Tucker3 model:

\enumerate{
    \item Classical Tucker3,
    \item Tucker3 for compositional data,
    \item Robust Tucker3 and
    \item Robust Tucker3 for compositional data.
}

This is controlled through the parameters \code{robust=TRUE} and \code{ilr=TRUE}.
}
\value{
An object of class "tucker3" which is basically a list with components:
\item{fit}{Fit value}
\item{fp}{Fit percentage}
\item{A}{Orthogonal loading matrix for the A-mode}
\item{B}{Orthogonal loading matrix for the B-mode}
\item{Bclr}{Orthogonal loading matrix for the B-mode, clr transformed.
    Available only if \code{ilr=TRUE} (default), otherwise NULL}
\item{C}{Orthogonal loading matrix for the C-mode}
\item{GA}{Core matrix, which describes the relation between \code{A}, \code{B} and \code{C},
    unfolded in A-form. The largest squared elements of the core matrix indicate the
    most important factors in the model of \code{X}.}
\item{iter}{Number of iterations}
\item{RD}{Residual distances}
\item{flag}{The observations whose residual distance \code{RD} is larger than cutoff.RD can be
    considered as outliers and receive a flag equal to zero.
    The regular observations receive a flag 1}
}

\references{
    Tucker, L.R. (1966). Some mathematical notes on three-mode 
        factor analysis. \emph{Psychometrika}, 31: 279--311.

    Egozcue J.J., Pawlowsky-Glahn, V.,  Mateu-Figueras G. and 
        Barcel'o-Vidal, C. (2003). Isometric logratio transformations 
        for compositional data analysis. \emph{Mathematical Geology}, 35(3): 279--300.
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

## Use toArray to turn to a 3-way array
##
X <- toArray(BusN, 7, 5, 37)
dimnames(X) <- list(dimnames(Bus)[[1]], 
                    substr(dimnames(Bus)[[2]][1:5], 1, 1), 
                    substr(dimnames(Bus)[[2]][seq(1, ncol(Bus), 5)], 3, 8))
dimnames(X)

## Perform classical Tucker3 (non-robust option and no
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
plot(res, which="tjplot", main="Trajectory biplot")

}
\keyword{ Multivariate }
\keyword{ Robust }
