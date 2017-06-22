\name{permute} 
\alias{permute}
\title{Permutation of a matricized array}
\description{
 Permutes the matricized (\code{n} \code{x} \code{m} \code{x} \code{p}) array \code{X} to the matricized array \code{Y} of order (\code{m} \code{x} \code{p} \code{x} \code{n}).
} 
\usage{
 permute(X,n,m,p)
}
\arguments{
  \item{X}{Matrix (or data.frame coerced to a matrix) containing the matricized array}
  \item{n}{Number of \code{A}-mode entities of the array \code{X}}
  \item{m}{Number of \code{B}-mode entities of the array \code{X}}
  \item{p}{Number of \code{C}-mode entities of the array \code{X}}
}
\value{
  \item{Y}{Matrix containing the permuted matricized array}
} 
\references{H.A.L. Kiers (2000). Towards a standardized notation and terminology in multiway analysis. \emph{Journal of Chemometrics 14:105--122}.}
\examples{
X <- array(1:24, c(4,3,2))
dim(X)

## matricize the array

##matricized X with the A-mode entities in its rows
Xa <- unfold(X)
dim(Xa)
Xa

## matricized X with the B-mode entities in its rows
Xb <- permute(Xa, 4, 3, 2)
dim(Xb)
Xb

## matricized X with the C-mode entities in its rows
Xc <- permute(Xb, 3, 2, 4)
dim(Xc)
Xc
}
\keyword{array}
\keyword{multivariate}
\keyword{algebra}
