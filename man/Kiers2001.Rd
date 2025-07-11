% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datadoc.R
\docType{data}
\name{Kiers2001}
\alias{Kiers2001}
\title{Fictious data set of Kiers and Mechelen (2001).}
\format{
A three-way array with dimension \code{6 x 5 x 4}.
 The first dimension refers to the 6 persons. The second dimension
 refers to the five responce variables.
 The third dimension refers to four different situations.
}
\source{
Kiers HAL, Mechelen IV (2001) Three-way component analysis: Principles and
illustrative application. Psychological Methods 6(1):84--110
}
\usage{
data(Kiers2001)
}
\description{
A fictitious data presented by Kiers and Mechelen (2001) to illustrate Tucker3 
 model: A set of six persons with scores on five response variables for four 
 different situations. The response variables indicate to what extent each 
 individual displays an emotional, sensitive, caring, thorough, or accurate 
 behavior. The data set is represented as a 6 x 5 x 4 array. The data are 
 chosen such that they correspond to a Tucker3 model with \code{P=2, Q=2, R=2}.
}
\examples{

 data(Kiers2001)

 t3 <- Tucker3(Kiers2001, P=2, Q=2, R=2)
 t3
 t3$A
 t3$B
 t3$C
 t3$G
 plot(t3)
 plot(t3, which="jbplot", xlim=c(0, 2))

}
\references{
Todorov, V., Simonacci, V., Gallo, M., and Di Palma, M. (2025). Robust tools
 for three-way component analysis of compositional data: The R package
 rrcov3way. Behaviormetrika. In press.
}
\keyword{datasets}
