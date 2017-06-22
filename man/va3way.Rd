\name{va3way}
\alias{va3way}
\docType{data}
\title{
Manufacturing value added by technology intensity for several years
}
\description{
A three-way array containing manufacturing value added by 
technology intensity for 55 countries in the period 2000--2010.
UNIDO maintains a unique database containing key industrial 
statistics indicators for more than 160 countries in the 
world in the period 1963-2011: INDSTAT 2, available at 
    \url{stat.unido.org}. The data are organized according 
to the International Standard Industrial Classification of all 
economic activities (ISIC) Revision 3 at 2-digit level. The present data 
set was created by aggregating the 23 2-digit divisions into five groups 
    according to technology intensity, 
using the UNIDO derived classification (Upadhyaya, 2011). Then 55 countries were
selected which have relatively complete data in the period 2000--2010.

}
\usage{data(va3way)}
\format{
A three-way array with dimension 55x5x11. 
The first dimension refers to 55 countries. 
The second dimension refers to the five categories of technology intensity 
    described above. 
The third dimension refers to the years in the period 2000--2010.

}
\details{
Note that the values in the second mode (sectors) sum up to 
a constant - the total manufacturing value added of a country 
in a given year and thus the data set has a compositional character.
}
\source{
    \url{stat.unido.org}
}
\references{
    Upahdyaya S (2011). Derived classifications for industrial performance indicators.
    In \emph{Int. Statistical Inst.: Proc. 58th World Statistical Congress, 
    2011, Dublin (Session STS022)}.
    
    Upadhyaya S, Todorov V (2008). UNIDO Data Quality. UNIDO Staff Working Paper, Vienna.
}

\examples{
data(va3way)
ct <- 2
x <- va3way[ct,,]/1000000
plot(colnames(x), x[1,], ylim=c(min(x), max(x)), type="n", ylab="Manufacturing Value 
    Added in million USD", xlab="Years")
for(i in 1:nrow(x))
    lines(colnames(x), x[i,], col=i)
legend("topleft", legend=rownames(x), col=1:nrow(x), lwd=1)
title(paste("Coutnry: ", rownames(va3way[,,1])[ct]))


## Treat quickly and dirty the zeros in the data set (if any) 
## in order to be able to perform ilr transformation:

va3way[va3way==0] <- 0.001

res <- Tucker3(va3way)

##
## Not yet a print function
##
print(res$fit)
print(res$A)

## Print the core matrix
print(res$GA)

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
