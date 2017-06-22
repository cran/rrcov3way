\name{elind}
\alias{elind}
\docType{data}
\title{
    OECD Electronics Industries Data
}
\description{

OECD publishes comparative statistics of the export size of various 
sectors of the electronics industry: 
\enumerate{
  \item information science, 
  \item telecommunication products, 
  \item radio and television equipment, 
  \item components and parts, 
  \item electromedical equipment, and 
  \item scientific equipment. 
}

The data consist of specialisation indices of electronics industries 
of 23 European countries for the years 1973--1979. The specialization 
index is defined as the proportion of the monetary value of an electronic 
industry compared to the total export value of manufactured goods of 
a country compared to the similar proportion for the world as 
a whole (see D'Ambra, 1985, p. 249 and Kroonenberg, 2008, p.282).   

}
\usage{data(elind)}
\format{
A three-way array with dimension 23x6x7. 
The first dimension refers to 23 countries. 
The second dimension refers to the six indices of electronics industries. 
The third dimension refers to the years in the period 1978--1985.
}
\source{
	The data set is available from Pieter Kroonenberg's web site at: 
    \url{http://three-mode.leidenuniv.nl/data/electronicindustriesinfo.htm}
}
\references{
    D'Ambra, L. (1985). Alcune estensione dell'analisi in componenti principali 
        per lo studio dei sistemi evolutivi. Uno studio sul commercio internazionale 
        dell'elettronica. 
        In: Ricerche Economiche. 2. del Dipartimento di Scienze Economiche Ca'Foscari, Venezia.
    
    Kroonenberg PM (2008). Applied multiway data analysis. Wiley series in probability and
        statistics. John Wiley and Sons, Hoboken, NJ, p.282.
}
\examples{
data(elind)

res <- Parafac(elind, robust=FALSE, ilr=FALSE)

## Distance-distance plot
plot(res, which="dd", main="Distance-distance plot")

## Paired component plot, mode A
plot(res, which="comp", main="Paired component plot (mode A)")

## Paired component plot, mode B
plot(res, which="comp", mode="B", main="Paired component plot (mode B)") 

## Per-component plot
plot(res, which="percomp", comp=1, main="Per component plot") 

## all components plot
plot(res, which="allcomp", main="All components plot", legend.position="topright") 

}
\keyword{datasets}
