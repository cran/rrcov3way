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
Parafac(X, ncomp = 2, center = FALSE, 
    center.mode = c("A", "B", "C", "AB", "AC", "BC", "ABC"),
    scale=FALSE, scale.mode=c("B", "A", "C"), 
    const="none", conv = 1e-06, start="svd", maxit=10000, 
    robust = FALSE, coda.transform=c("none", "ilr", "clr"), 
    ncomp.rpca = 0, alpha = 0.75, robiter = 100, crit=0.975, trace = FALSE)
}
\arguments{
  \item{X}{3-way array of data}
  \item{ncomp}{Number of components}
  \item{center}{Whether to center the data}
  \item{center.mode}{If centering the data, on which mode to do this}
  \item{scale}{Whether to scale the data}
  \item{scale.mode}{If scaling the data, on which mode to do this}
    \item{const}{Optional constraints for each mode. Can be a three element character
    vector or a single character, one of \code{"none"} for no constraints (default),
    \code{"orth"} for orthogonality constraints, \code{"nonneg"} for nonnegativity constraints or
    \code{"zerocor"} for zero correlation between the extracted factors. For example,
    \code{const="orth"} means orthogonality constraints for all modes,
    while \code{const=c("orth", "none", "none")} sets the orthogonality constraint
    only for mode A.}
  \item{conv}{Convergence criterion, defaults to \code{1e-6}}
  \item{start}{Initial values for the A, B and C components. Can be \code{"svd"} 
    for starting point of the algorithm from SVD's, \code{"random"} for random
    starting point (orthonormalized component matrices or nonnegative matrices in
    case of nonnegativity constraint), or a list containing user specified components.}
  \item{maxit}{Maximum number of iterations, default is \code{maxit=10000}.}
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
    versions of the Parafac model: 
    
\enumerate{
    \item Classical Parafac,
    \item Parafac for compositional data,
    \item Robust Parafac and
    \item Robust Parafac for compositional data.
}
  
This is controlled though the paramters \code{robust=TRUE} and coda.transform=c("none", "ilr").
}
\value{
An object of class "parafac" which is basically a list with components:
\item{fit}{Fit value}
\item{fp}{Fit percentage}
\item{A }{Orthogonal loading matrix for the A-mode}
\item{B}{Orthogonal loading matrix for the A-mode}
\item{Bclr}{Orthogonal loading matrix for the B-mode, clr transformed. 
    Available only if coda.transform="ilr", otherwise NULL}
\item{C}{Orthogonal loading matrix for the C-mode}
\item{Xhat}{(Robustly) reconstructed array}
\item{const}{Optional constraints (same as the input parameter)}
\item{iter}{Number of iterations}
\item{rd}{Residual distances}
\item{sd}{Score distances}
\item{flag}{The observations whose residual distance \code{rd} is larger than \code{cutoff.rd} or 
    score distance \code{sd} is larger than \code{cutoff.sd}, can be 
    considered outliers and receive a flag equal to zero. 
    The regular observations receive a flag 1}
\item{robust}{The paramater \code{robust}, whether robust method is used or not}    
\item{coda.transform}{Which coda transformation is used, can be \code{coda.transform=c("none", "ilr", "clr")}.}    
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
res <- Parafac(ulabor, robust=TRUE, coda.transform="ilr")
res

## Plot Orthonormalized A-mode component plot
plot(res, which="comp", mode="A", main="Component plot, A-mode")

## Plot Orthonormalized B-mode component plot
plot(res, which="comp", mode="B", main="Component plot, B-mode")

## Plot Orthonormalized C-mode component plot
plot(res, which="comp", mode="C", main="Component plot, C-mode")


}
\keyword{ Multivariate }
\keyword{ Robust }

