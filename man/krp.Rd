\name{krp}
\alias{krp}
\title{
Function to compute Khatri-Rao product
}
\description{
The function \code{krp(A,B)} returns the Khatri-Rao product of two matrices \code{A} and \code{B}, of
dimensions I x K and J x K respectively. The result is an IJ x K matrix formed by the matching 
column-wise Kronecker products, i.e. the k-th column of the Khatri-Rao product is defined as \code{kronecker(A[, k], B[, k])}.
}
\usage{
krp(A, B)
}
\arguments{
  \item{A}{
	Matrix of order I x K.
}
  \item{B}{
	Matrix of order J x K.
}
}
\value{
	The IJ x K matrix of columnwise Kronecker products.
}
\author{
Valentin Todorov \email{valentin.todorov@chello.at}
}
\note{
	\code{A} and \code{B} must have the same number of columns.
}
\references{
	Khatri, C.G. and Radhakrishna Rao, C. (1968). Solutions to Some Functional Equations and Their Applications to Characterization of Probability Distributions. \emph{Sankhya: Indian J. Statistics (Series A)}, 30:167-180.
	Smilde, A., Bro R. and Gelardi, P. (2004). Multi-way Analysis: Applications in Chemical Sciences, Chichester:Wiley
}
\examples{
A <- matrix(1, nrow=5, ncol=2)
B <- matrix(1:4, ncol=2)
krp(A,B)

A <- matrix(1:10, ncol=2, byrow=TRUE)
B <- diag(1,2)
krp(A,B)
}

