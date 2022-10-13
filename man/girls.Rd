\name{girls}
\alias{girls}
\docType{data}
\title{Sempe girls' growth curves data
}
\description{
Thirty girls selected from a French auxiological study (1953-1975) to get 
insight into the physical growth patterns of children from ages four to fifteen, 
Sempe (1987). They were measured yearly between the ages 4 and 15 on the following eight variables:
\enumerate{
\item 	weight = Weight
\item 	length = Length
\item	crump = Crown-rump length
\item	head = Head circumference
\item	chest = Chest circumference
\item	arm = Arm
\item	calf = Calf
\item	pelvis = Pelvis
}

The data set is three way data array of size 30 (girls) x 8 (variables) x 12 (years). 
}
\usage{data("girls")}
\format{
  The format is a three way array with the following dimensions:
    The first dimension refers to 30 girls. 
    The second dimension refers to the eight variables measured on the girls. 
    The third dimension refers to the years -- 4 to 15.  
}
\details{
The data are generally preprocessed as standard multiway profile data. 
For details see Kroonenberg (2008), Chapters 6 and 15.
}
\source{
	The data sets are available from Pieter Kroonenberg's web site at:
    \url{https://web.universiteitleiden.nl/fsw/three-mode/data/girlsgrowthcurvesinfo.htm}
}
\references{
    Sempe, M. (1987). Multivariate and longitudinal data on growing children: 
    Presentation of the French auxiological survey. In J.Janssen et al. Data analysis. 
    The Ins and Outs of solving real problems (pp. 3-6). New York: Plenum Press.
    
    Kroonenberg (2008). Applied multiway data analysis. Wiley series 
        in probability and statistics. Hoboken NJ, Wiley.
}
\examples{
data(girls)
str(girls)
## Center the data in mode A and find the "average girl"
center.girls <- do3Scale(girls, center=TRUE, only.data=FALSE)
X <- center.girls$x
center <- center.girls$center
average.girl <- as.data.frame(matrix(center, ncol=8, byrow=TRUE))
dimnames(average.girl) <- list(dimnames(X)[[3]], dimnames(X)[[2]])

## Divide these variables by 10 to reduce their range
average.girl$weight <- average.girl$weight/10
average.girl$length <- average.girl$length/10
average.girl$crrump <- average.girl$crrump/10

average.girl 
p <- ncol(average.girl)
plot(rownames(average.girl), average.girl[,1], ylim=c(min(average.girl), 
    max(average.girl)), type="n", xlab="Age", ylab="")
for(i in 1: p)
{
    lines(rownames(average.girl), average.girl[,i], lty=i, col=i) 
    points(rownames(average.girl), average.girl[,i], pch=i, col=i) 
}    
legend <- colnames(average.girl)
legend[1] <- paste0(legend[1], "*")
legend[2] <- paste0(legend[3], "*")
legend[3] <- paste0(legend[4], "*")
legend("topleft", legend=legend, col=1:p, lty=1:p, pch=1:p)
}
\keyword{datasets}
