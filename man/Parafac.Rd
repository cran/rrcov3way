\name{Parafac}
\alias{Parafac}
\alias{print.parafac}
\title{
    Robust Parafac estimator for compositional data
}
\description{
Compute a robust Parafac model for compositional data
}
\usage{
Parafac(X, ncomp = 2, conv = 1e-06, center = FALSE, 
    scale=FALSE, scale.mode=c("B", "A", "C"), 
    orth=c(), robust = FALSE, ilr = FALSE, ncomp.rpca = 2, 
    alpha = 0.75, maxiter = 100, crit=0.975, trace = FALSE)
}
\arguments{
  \item{X}{3-way array of data}
  \item{ncomp}{Number of components}
  \item{conv}{Convergence criterion, defaults to \code{1e-6}}
  \item{center}{Whether to center the data}
  \item{scale}{Whether to scale the data}
  \item{scale.mode}{If centering and/or scaling the data, on which mode to do this}
  \item{orth}{Orthogonality constraints}
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
    versions of the Parafac model: 
    
\enumerate{
    \item Classical Parafac,
    \item Parafac for compositional data,
    \item Robust Parafac and
    \item Robust Parafac for compositional data.
}
  
This is controlled though the paramters \code{robust=TRUE} and \code{ilr=TRUE}.
}
\value{
An object of class "parafac" which is basically a list with components:
\item{fit}{Fit value}
\item{fp}{Fit percentage}
\item{A }{Orthogonal loading matrix for the A-mode}
\item{B}{Orthogonal loading matrix for the A-mode}
\item{Bclr}{Orthogonal loading matrix for the B-mode, clr transformed. 
    Available only if \code{ilr=TRUE} (default), otherwise NULL}
\item{C}{Orthogonal loading matrix for the C-mode}
\item{Xhat}{Robust reconstructed array}
\item{iter}{Number of iterations}
\item{RD}{Residual distances}
\item{flag}{The observations whose residual distance \code{RD} is larger than cutoff.RD can be 
    considered as outliers and receive a flag equal to zero. 
    The regular observations receive a flag 1}
\item{robust}{The paramater \code{robust}, whether robust method is used or not}    
\item{ilr}{The paramater \code{ilr}, whether ilr transformation is used or not}    
}
\references{
	Harshman, R.A. (1970). Foundations of Parafac procedure: 
       models and conditions for an "explanatory" multi-mode factor 
       analysis. \emph{UCLA Working Papers in Phonetics}, 16: 1--84.

	Engelen, S., Frosch, S. and Jorgensen, B.M. (2009). A fully 
       robust PARAFAC method analyzing fluorescence data. 
       \emph{Journal of Chemometrics}, 23(3): 124--131.

	Kroonenberg, P.M. (1983).Three-mode principal component analysis: 
       Theory and applications (Vol. 2), DSWO press.

	Rousseeuw, P.J. and Driessen, K.V. (1999). A fast algorithm for 
        the minimum covariance determinant estimator. 
        \emph{Technometrics}, 41(3): 212--223.

	Egozcue J.J., Pawlowsky-Glahn V., Mateu-Figueras G. and 
        Barcel'o-Vidal, C. (2003). Isometric logratio transformations 
        for compositional data analysis. \emph{Mathematical Geology}, 35(3): 279-300
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
res <- Parafac(va3way)
res
print(res$fit)
print(res$A)

## Distance-distance plot
plot(res, which="dd", main="Distance-distance plot")

data(ulabor)
res <- Parafac(ulabor, robust=TRUE, ilr=TRUE)
res

## Plot Orthonormalized A-mode component plot
plot(res, which="comp", mode="A", main="Component plot, A-mode")

## Plot Orthonormalized B-mode component plot
plot(res, which="comp", mode="B", main="Component plot, B-mode")

## Plot Orthonormalized B-mode component plot
plot(res, which="comp", mode="C", main="Component plot, C-mode")


}
\keyword{ Multivariate }
\keyword{ Robust }

