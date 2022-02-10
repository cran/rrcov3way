\name{ulabor}
\alias{ulabor}
\docType{data}
\title{
Undeclared labor by region in Italy
}
\description{
The dataset contains the undeclared labor in thousands work units. The data
originate from Italy and are recorded at a regional level over a certain 
time horizon for five macroeconomic activities defined according to 
NACE Rev. 1.1 classification.
}
\usage{data("ulabor")}
\format{
A three-way array with dimension 22x5x5. The first dimension 
refers to 22 regions in Italy. The second dimension refers to 
the 5 economic activities. The third dimension refers to 
the years in the period 2001-2009.
}
%\details{}
\source{
    ISTAT (2011). Note metodologiche, la misura dell'occupazione non 
    regolare nelle stime di contabilita nazionale [online]. Roma. 
    %%%Available at: \url{http://www.istat.it/it/archivio/39522}.
}
\references{
    ISTAT (2011). Note metodologiche, la misura dell'occupazione non 
    regolare nelle stime di contabilita nazionale [online]. Roma. 
    %%%Available at: \url{http://www.istat.it/it/archivio/39522}.
    
    Di Palma M.A., Filzmoser P., Gallo M. and Hron, K. (2016). 
        A robust CP model for compositional data, submitted.
}
\examples{

data(ulabor)
dim(ulabor)
str(ulabor)

## Plot robust and non-robust DD-plots of the ilr-transformed data
usr <- par(mfrow=c(1,2))
res1 <- Parafac(ulabor, robust=TRUE, coda.transform="ilr")
res2 <- Parafac(ulabor, coda.transform="ilr")
plot(res1)
plot(res2)
par(usr)

\dontrun{

## Plot Orthonormalized A-mode component plot
res <- Parafac(ulabor, robust=TRUE, coda.transform="ilr")
plot(res, which="comp", mode="A", main="Component plot, A-mode")

## Plot Orthonormalized B-mode component plot
plot(res, which="comp", mode="B", main="Component plot, B-mode")

## Plot Orthonormalized B-mode component plot
plot(res, which="comp", mode="C", main="Component plot, C-mode")


## Per component plot
## adapted for the example and only for robust, ilr transformed model
##
##
res <- Parafac(ulabor, robust=TRUE, coda.transform="ilr")

plot(res, which="percomp")              # component 1
plot(res, which="percomp", comp=2)      # component 2
}
}
\keyword{datasets}
