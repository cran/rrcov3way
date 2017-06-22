\name{unfold} 
\alias{unfold}
\title{Matrix unfolding}
\description{
 Conducts matricizations of a three-way array into matrices according to the selected mode.
}
\usage{
    unfold(x, mode=c("A", "B", "C"))
}
\arguments{
    \item{x}{Array to be unfolded}
    \item{mode}{the selected mode for unfolding}
}
\value{
A matrix represnting the input array, according to the selected mode:
\itemize{
\item Mode=A: \code{B}-mode entities are nested within \code{C}-mode entities (all the frontal slices of the array next to each other)
item Mode=B: \code{C}-mode entities nested within \code{A}-mode entities (all the horizontal slices of the array next to each other)
item Mode C: \code{A}-mode entities nested within \code{B}-mode entities (all the lateral slices of the array next to each other)
}
}

\references{H.A.L. Kiers (2000). Towards a standardized notation and terminology in multiway analysis. \emph{Journal of Chemometrics 14:105--122}.}

\author{
Valentin Todorov \email{valentin.todorov@chello.at}
}
\examples{

(X <- array(1:24, c(4,3,2)))
dim(X)

## matricize the array

##  matricized X with the A-mode entities in its rows
##  all the frontal slices of the array next to each other
##
(Xa <- unfold(X))
dim(Xa)

## matricized X with the B-mode entities in its rows
##  all the horizontal slices of the array next to each other
##
(Xb <- unfold(X, mode="B"))
dim(Xb)

## matricized X with the C-mode entities in its rows
##  all the lateral slices of the array next to each other
##
(Xc <- unfold(X, mode="C"))
dim(Xc)
}
\keyword{array}
\keyword{multivariate}
\keyword{algebra}
