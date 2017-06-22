\name{Kojima}
\alias{Kojima.boys}
\alias{Kojima.girls}
\docType{data}
\title{
Parental behaviour in Japan
}
\description{
The data are drawn from a study (Kojima, 1975) of the perception of 
parental behaviour by parents and their children. Two data sets, boys and girls are available
as Kojima.boys and Kojima.girls. 

\itemize{
\item Boys data were analysed in Kroonenberg (2008)
\item Girls data were analysed in Kroonenberg, Harshman, & Murakami (2009).
}
}
\usage{
data(Kojima)

}
\format{
Both data sets are three dimensional arrays:
\itemize{
    \item boys: 150 x 18 x 4
    \item girls: 153 x 18 x 4
}
The rows (1st mode) are 150 Japanese sons/153 Japanese daughters.
The columns (2nd mode) are 18 scales (Acceptance, Child centerness, Possesiveness, etc.).
The slices (3rd mode) are the 4 judgements (See Details for explanation).

}
\details{
The \code{boys} data are ratings expressing the judgments of parents with respect to their 
own behaviour toward their sons, and the judgments of their sons with respect to their 
parents. Thus, there are four conditions: 
\itemize{
\item Father-Own behaviour (F-F), 
\item Mother-Own behaviour (M-M), 
\item Son-Father (B-F), 
\item Son-Mother (B-M). 
}
The judgments involved 150 middle-class Japanese eighth-grade boys on the 18 subscales of the inventory. 
Thus, the data set consists of a 150 (Sons) x 18 (scales) x 4 (judgment combinations) 
data array.

Similarly, the girls data are ratings expressing the judgments of parents with respect to 
their own behaviour toward their daughters, and the judgments of their daughters with 
respect to their parents. Thus, there are four conditions: 
\itemize{
\item Father-Own behaviour (F-F), 
\item Mother-Own behaviour (M-M),
\item Daughter-Father (G-F), 
\item Daughter-Mother (G-M).
} 
The judgments involved 153 middle-class Japanese eighth-grade girls on the 18 subscales of 
the inventory. Thus, the data set consists of a 153 (Daughters) x 18 (scales) x 4 (judgment combinations) 
data array.

\emph{Preprocessing}
Given that the data are three-way profile data they are treated in the standard manner:
centering per occasion-variable combination and by normalising the data after centring 
per lateral slice i.e. per scale over all sons/daughters x judges combinations. 
For details see Kroonenberg (2008), Chapter 13.
}
\source{
The data sets are available from the Pieter Kroonenberg's web site at \url{http://three-mode.leidenuniv.nl/}.
}
\references{

    Kojima, H. (1975). Inter-battery factor analysis of parents' and children's reports of 
    parental behavior. Japanese Psychological Bulletin, 17, 33-48 (in Japanese).
    
    Kroonenberg, P. M. (2008). Applied multiway data analysis. Wiley series in probability and statistics. Wiley, Hoboken NJ.
    
    Kroonenberg, P. M., Harshman, R. A, & Murakami, T. (2009). Analysing three-way profile data using 
    the Parafac and Tucker3 models illustrated with views on parenting. Applied Multivariate 
    Research, 13:5-41. PDF available at: http://www.phaenex.uwindsor.ca/ojs/leddy/index.php/AMR/article/viewFile/2833/2271
}
\examples{
    data(Kojima)
    dim(Kojima.boys)
    dim(Kojima.girls)
}
\keyword{datasets}
