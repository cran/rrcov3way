\name{congruence}
\alias{congruence}
\alias{congruence}
\title{Coefficient of factor congruence (phi)}
\description{
Computes the Tucker's congruence (phi) coefficients among two sets of factors.
}
\usage{
congruence(x, y=NULL)
}
\arguments{
  \item{x}{ A vector or matrix of factor loadings}
  \item{y}{ A vector or matrix of factor loadings (may be empty)}
}
\details{Find the Tucker's coefficient of congruence between two sets of factor loadings. 

Factor congruences are the cosines of pairs of vectors defined by the loadings matrix and based at the origin.  Thus, for loadings that differ only by a scaler (e.g. the size of the eigen value), the factor congruences will be 1.

For factor loading vectors of X and Y the measure of factor congruence, phi, is 
\deqn{
\phi = \frac{\sum X Y}{\sqrt{\sum(X^2)\sum(Y^2)}} 
.}{phi = sum(X Y)/sqrt(sum(X^2) sum(X^2)) }

}
\value{A matrix of factor congruences.
}

\references{L.R Tucker (1951). A method for synthesis of factor analysis studies. \emph{Personnel Research Section Report No. 984}. Department of the Army, Washington, DC.
}
\author{
 \email{valentin.todorov@chello.t}
}
 
\examples{
X <- getLoadings(PcaClassic(delivery))
Y <- getLoadings(PcaHubert(delivery, k=3))
round(congruence(X,Y),3)
}
\keyword{ multivariate }
\keyword{ models }
