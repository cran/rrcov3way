\name{do3Scale}
\alias{do3Scale}
\alias{do3Scale.default}
\alias{do3Scale.parafac}
\alias{do3Scale.tucker3}
\title{
Centering and scaling
}
\description{
Centering and/or normalization of a three way array or a matricized 
array across one mode (modes indicated by "A", "B" or "C").
}
\usage{
    \method{do3Scale}{tucker3}(x, renorm.mode = c("A", "B", "C"), \dots)
    \method{do3Scale}{parafac}(x, renorm.mode = c("A", "B", "C"), \dots)
    \method{do3Scale}{default}(x, center = FALSE, scale = FALSE, 
        center.mode = c("A", "B", "C", "AB", "AC", "BC", "ABC"), 
        scale.mode = c("B", "A", "C"), 
        only.data=TRUE, \dots)
}
\arguments{
  \item{x}{Three dimensional array of order (I x J x K) or a matrix 
    (or data.frame coerced to a matrix) of order (I x JK) containing 
    the matricized array (frontal slices)}
  \item{center}{Whether and how to center the data. Can be \code{NULL}, 
  logical \code{TRUE} or \code{FALSE}, function or a numeric vector with 
  length corresponding to the number of columns in the corresponding mode.
  If \code{center=TRUE}, \code{mean()} is used; default is \code{center=FALSE}.}
  \item{scale}{Whether and how to scale the data. Can be \code{NULL}, 
  logical \code{TRUE} or \code{FALSE}, function or a numeric vector with 
  length corresponding to the number of columns in the corresponding mode. 
  If \code{scale=TRUE}, \code{sd()} is used; default is \code{scale=FALSE}.}
  \item{center.mode}{Across which mode to center. Default is \code{center.mode="A"}}
  \item{scale.mode}{Within which mode to scale. Default is \code{scale.mode="B"}}
  \item{renorm.mode}{Within which mode to renormalize a Parafac or Tucker3 solution. 
  See in Details how this is performed for the different models. 
  Default is \code{renorm.mode="A"}}
  \item{only.data}{Whether to return only the centered/scaled data or also 
    the center and the scale themselves. Default is \code{only.data=TRUE}}
  \item{\dots}{potential further arguments passed to lower level functions.}
}
\value{
A named list, consisting of the centered and/or scaled data, a center vector, 
a scale vector and the mode in which the data were centered/scaled. 
}
\references{
    Kiers, H.A.L. (2000).Towards a standardizrd notationand terminology 
        in multiway analysis. \emph{Journal of Chemometrics}, 14:105-122.
        
    Kroonenberg, P.M. (1983).Three-mode principal component analysis: 
        Theory and applications (Vol. 2), DSWO press.
}
\author{
Valentin Todorov \email{valentin.todorov@chello.at} and
Maria Anna Di Palma \email{madipalma@unior.it} and
Michele Gallo \email{mgallo@unior.it}
}

\examples{
    data(elind)
    (x1 <- do3Scale(elind, center=TRUE, scale=TRUE))
    (x2 <- do3Scale(elind, center=TRUE, scale=TRUE, center.mode="B"))
    (x3 <- do3Scale(elind, center=TRUE, scale=TRUE, center.mode="C", scale.mode="C"))
}
\keyword{array}
\keyword{multivariate}
\keyword{algebra}

