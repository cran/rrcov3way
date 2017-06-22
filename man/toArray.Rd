\name{toArray}
\alias{toArray}
\title{
Matrix to array conversion
}
\description{
Restore an array from its matricization with all the frontal slices of the array next to each other (mode="A")
}
\usage{
toArray(x, n, m, r, mode = c("A", "B", "C"))
}
\arguments{
  \item{x}{ Matrix (or data.frame coerced to a matrix) containing the elements of the frontal slices of an array}
  \item{n}{number of A-mode elements}
  \item{m}{number of B-mode elements}
  \item{r}{number of C-mode elements}
  \item{mode}{in which mode is the matricized array}
}
\value{
Three way array
}
\references{
H.A.L. Kiers (2000). Towards a standardized notation and terminology in multiway analysis.  \emph{Journal of Chemometrics}, 14: 105--122.
}
\author{
Valentin Todorov \email{valentin.todorov@chello.at} and
Maria Anna Di Palma \email{madipalma@unior.it} and
Michele Gallo \email{mgallo@unior.it}
}

\examples{
    data(elind)
    di <- dim(elind)
    toArray(unfold(elind), di[1], di[2], di[3])
}
\keyword{array}
\keyword{multivariate}
\keyword{algebra}
