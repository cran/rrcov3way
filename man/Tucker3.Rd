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
Tucker3(X, P = 2, Q = 2, R = 2, 
    center = FALSE, center.mode = c("A", "B", "C", "AB", "AC", "BC", "ABC"), 
    scale = FALSE, scale.mode = c("B", "A", "C"), 
    conv = 1e-06, start="svd",
    robust = FALSE, coda.transform=c("none", "ilr", "clr"), 
    ncomp.rpca = 0, alpha = 0.75, robiter=100, crit=0.975, trace = FALSE)
}
\arguments{
  \item{X}{3-way array of data}
  \item{P}{Number of A-mode components}
  \item{Q}{Number of B-mode components}
  \item{R}{Number of C-mode components}
  \item{center}{Whether to center the data}
  \item{center.mode}{If scaling the data, on which mode to do this}
  \item{scale}{Whether to scale the data}
  \item{scale.mode}{If centering the data, on which mode to do this}
  \item{conv}{Convergence criterion, defaults to \code{1e-6}}
  \item{start}{Initial values for the A, B and C components. Can be \code{"svd"} 
    for starting point of the algorithm from SVD's, \code{"random"} for random
    starting point (orthonormalized component matrices), or a list containing user specified components.}
  \item{robust}{Whether to apply a robust estimation}
  \item{coda.transform}{If the data are a composition, use an \emph{ilr} or \emph{clr} transformation. 
    Default is non-compositional data, i.e. \code{coda.transform="none"}}
  \item{ncomp.rpca}{Number of components for robust PCA}
  \item{alpha}{Measures the fraction of outliers the algorithm should
    resist. Allowed values are between 0.5 and 1 and the default is 0.75}
  \item{robiter}{Maximal number of iterations for robust estimation}
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

This is controlled through the parameters \code{robust=TRUE} and \code{coda.transform="ilr"}.
}
\value{
An object of class "tucker3" which is basically a list with components:
\item{fit}{Fit value}
\item{fp}{Fit percentage}
\item{A}{Orthogonal loading matrix for the A-mode}
\item{B}{Orthogonal loading matrix for the B-mode}
\item{Bclr}{Orthogonal loading matrix for the B-mode, clr transformed.
    Available only if \code{coda.transform="ilr"}, otherwise NULL}
\item{C}{Orthogonal loading matrix for the C-mode}
\item{GA}{Core matrix, which describes the relation between \code{A}, \code{B} and \code{C},
    unfolded in A-form. The largest squared elements of the core matrix indicate the
    most important factors in the model of \code{X}.}    
\item{iter}{Number of iterations}
\item{rd}{Residual distances}
\item{sd}{Score distances}
\item{flag}{The observations whose residual distance \code{RD} is larger than cutoff.RD can be
    considered as outliers and receive a flag equal to zero.
    The regular observations receive a flag 1}
\item{robust}{The paramater \code{robust}, whether robust method is used or not}    
\item{coda.transform}{The input paramater \code{coda.transform}, what trnasofrmation for compositional data was used}    
\item{La}{Diagonal matrix containing the \emph{intrinsic eigenvalues} for A-mode}
\item{Lb}{Diagonal matrix containing the \emph{intrinsic eigenvalues} for B-mode}
\item{Lc}{Diagonal matrix containing the \emph{intrinsic eigenvalues} for C-mode}
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
