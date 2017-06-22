\name{Arno}
\alias{Arno}
\docType{data}
\title{
    Chemical composition of water in the main stream of Arno river
}
\description{
    Chemical composition of water in the main stream of Arno river.
}
\usage{data("Arno")}
\format{
A three-way array with dimension 23x11x4. 
The first dimension refers to 23 distances from the spring. 
The second dimension refers to the 11 chemical compositions. 
The third dimension refers to the time of collection - four occasions.
}
\details{
    The Arno data example was used in Gallo and Buccinati (2013) to illustrate 
    a particular version of the Tucker model, known as the weighted principal component analysis.
    The Tucker3 results are usually given in the form of tables or plots and in this work 
    for the representation of the Tucker3 results of logratio data, is proposed to use 
    one-mode plots, clr-joint biplots (Gallo, 2015), and trajectory plots.
}
\source{
    Nisi B., Vaselli O., Buccianti A., Minissale A., Delgado-Huertas A., Tassi F., Montegrossi G. (2008). Geochemical and isotopic investigation of the dissolved load in the running waters from 
    the Arno valley: evaluation of the natural and anthropogenic input. 
    In Memorie Descrittive della Carta Geologica d'Italia, Nisi (eds.), 79: 1-91.
    
    Nisi B., Buccianti A., Vaselli O., Perini G., Tassi F., Minissale A., Montegrossi G. (2008) 
    Hydrogeochemistry and strontium isotopes in the Arno river basin (Tuscany, Italy): 
    Constraints on natural controls by statistical modeling. \emph{Journal of Hydrology} 360: 166-183.
}
\references{

    Gallo M. and Buccianti A. (2013). Weighted principal component analysis for
    compositional data: application example for the water chemistry of 
    the Arno river (Tuscany, central Italy), \emph{Environmetrics}, 24(4):269-277.


Gallo M. (2015). Tucker3 model for compositional data. 
    \emph{Communications in Statistics-Theory and Methods}, 44(21):4441-4453.
}
\examples{
data(Arno)
dim(Arno)               # [1] 23 11  4
dim(Arno[,,1])          # [1] 23 11
rownames(Arno[,,1])     # the 23 distances from the spring
colnames(Arno[,,1])     # the 11 chemical compositions
dim(Arno[,1,])          # [1] 23  4
colnames(Arno[,1,])     # the four occasions

res <- Tucker3(Arno, robust=FALSE, ilr=TRUE)
res

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
