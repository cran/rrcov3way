\name{waterquality}
\alias{waterquality}
\docType{data}
\title{
    Water quality data in Wyoming, USA
}
\description{
    Water quality data for three years of seasonal compositional groundwater 
    chemistry data for 14 wells at a study site in Wyoming, USA. 
    Routine water quality monitoring typically involves measurement of J parameters and 
    constituents measured at I number of static locations at K sets of seasonal occurrences.
}
\usage{data("waterquality")}
\format{
A three-way array with dimension 14x12x10. 
The first dimension refers to 14 wells at a study site in Wyoming, USA. 
The second dimension refers to the ten most reactive and indicative dissolved
constituents at the site: B, Ba, Ca, Cl, K, Mg, Na, Si, Sr, and SO4. In addition,
the concentration of water in each sample was calculated. 
The third dimension refers to the time of collection - ten occasions.
}
%\details{}
%\source{}
\references{
Engle, M.A., Gallo, M., Schroeder, K.T., Geboy, N.J., Zupancic, J.W., (2014). Three-way compositional analysis of water quality monitoring data. \emph{Environmental and Ecological Statistics}, 21(3):565-581.
}
\examples{
data(waterquality)
dim(waterquality)               # [1] 14 12 10
dim(waterquality[,,1])          # [1] 14 12
rownames(waterquality[,,1])     # the 14 wells
colnames(waterquality[,,1])     # the 12 chemical compositions
dim(waterquality[,1,])          # [1] 14 10
colnames(waterquality[,1,])     # the ten occasions

(res <- Tucker3(waterquality, robust=FALSE, ilr=TRUE))

## Distance-distance plot
plot(res, which="dd", main="Distance-distance plot")
## Paired component plot, mode A
plot(res, which="comp", main="Paired component plot (mode A)")

## Paired component plot, mode B
plot(res, which="comp", mode="B", main="Paired component plot (mode B)") 

## Joint biplot
plot(res, which="jbplot", main="Joint biplot") 

## Trajectory
plot(res, which="tjplot", main="Trajectory biplot") 
}
\keyword{datasets}
