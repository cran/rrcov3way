\name{do3Postprocess}
\alias{do3Postprocess}
\alias{do3Postprocess.default}
\alias{do3Postprocess.parafac}
\alias{do3Postprocess.tucker3}
\alias{tallArray}
\alias{wideArray}
\alias{tall2wide}
\alias{coordinates}
\alias{coordinates.parafac}
\alias{coordinates.tucker3}
\alias{weights}
\alias{weights.parafac}
\alias{weights.tucker3}
\alias{reflect}
\alias{reflect.parafac}
\alias{reflect.tucker3}
\alias{reorder}
\alias{reorder.parafac}
\alias{reorder.tucker3}
\alias{is.orthogonal}
\alias{is.orthonormal}
\title{
    Postprocessing: renormalization, reflection and reordering; access to some of the components of the model.
}
\description{
    The estimated model will be renormalized, reflected (change of sign) or the components will be reordered.
    Functions that provide access to some components of the model: coordinates, weights.
}
\usage{
    \method{do3Postprocess}{tucker3}(x, reflectA, reflectB, reflectC, reorderA, reorderB, reorderC, ...)
    \method{do3Postprocess}{parafac}(x, reflectA, reflectB, reflectC, reorder, ...)
    \method{coordinates}{parafac}(x, mode = c("A", "B", "C"), type = c("normalized", "unit", "principal"), ...)
    \method{coordinates}{tucker3}(x, mode = c("A", "B", "C"), type = c("normalized", "unit", "principal"), ...)
    \method{weights}{parafac}(object, ...)
    \method{weights}{tucker3}(object, mode = c("A", "B", "C"), ...)
    \method{reflect}{parafac}(x, mode = c("A", "B", "C"), rsign = 1, ...)
    \method{reflect}{tucker3}(x, mode = c("A", "B", "C"), rsign = 1, ...)
}
\arguments{
  \item{x}{Tucker3 or Parafac solution}
  \item{object}{Tucker3 or Parafac solution (alternative of \code{x} for the generic function \code{weights()})}
  \item{reflectA}{How to handle the signs of the components of mode A - can be a single number or a vector with length of the number of components of A}
  \item{reflectB}{How to handle the signs of the components of mode B - can be a single number or a vector with length of the number of components of B}
  \item{reflectC}{How to handle the signs of the components of mode C - can be a single number or a vector with length of the number of components of C}
  \item{reorder}{How to reorder the components of a Parafac solution - a vector with length of the number of components}
  \item{reorderA}{How to reorder the components of mode A - a vector with length of the number of components of A giving the new order}
  \item{reorderB}{How to reorder the components of mode B - a vector with length of the number of components of B giving the new order}
  \item{reorderC}{How to reorder the components of mode C - a vector with length of the number of components of C giving the new order}
  \item{mode}{For which mode to provide the coordinates or weights. Default is mode A}
  \item{type}{Which type of coordinates to provide. Default is "normalized"}
  \item{rsign}{How to change the sign of the components of the given mode. Can be a single number or a vector with length of the number of components of the corresponding mode.}
  \item{\dots}{Potential further arguments passed to lower level functions.}
}
\value{
The output value of do3Postproces() is the postprocessed solution, Parafac or Tucker3. The output of \code{weights()} and \code{coordinates()} are the respoective values.
}
\references{
REFERENCES
}
\author{
Valentin Todorov \email{valentin.todorov@chello.at} and
Maria Anna Di Palma \email{madipalma@unior.it} and
Michele Gallo \email{mgallo@unior.it}
}

\examples{
    data(elind)
    x1 <- do3Scale(elind, center=TRUE, scale=TRUE)
    cp <- Parafac(x1, ncomp=3, const=c("orth", "none", "none"))
    cp$B
    cp1 <- do3Postprocess(cp, reflectB=-1)      # change the sign of all components of B
    cp$B
    weights(cp1)
    coordinates(cp1)
    coordinates(cp1, type="principal")
}
\keyword{array}
\keyword{multivariate}
\keyword{algebra}

