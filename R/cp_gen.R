##
##
##  VT::11.01.2024
##
##  I, J, K - number of rows, columns and (obsevrations, variables and occasions)
##  nsim    - Number of generated data sets
##  nf      - number of PARAFAC components
##  noise
##  noise1
##  Acol, Bcol, Ccol    -
##  congA, congB, congC -
##  eps     - percent contamination
##  type    - type of outliers
##  c1, c2  - parameters for outlier generation
##  silent  - do not issue warnings
##
##  roxygen2::roxygenise("c:/Users/valen/OneDrive/MyRepo/R/rrcov3way", load_code=roxygen2:::load_installed)


#'  Generate PARAFAC data sets, optionally with outliers
#'
#' @description Generates \code{nsim} data sets according to the given parameters.
#'  If \code{eps > 0}, the specified fraction of random outliers of the identified by the 
#'  parameter \code{type} type are added to the data sets.
#' @param I number of observations
#' @param J number of variables
#' @param K number of occasions
#' @param nsim number of data sets to generate
#' @param nf number of PARAFAC components
#' @param noise level of homoscedastic (HO) noise
#' @param noise1 level of heteroscedastic (HE) noise
#' @param Acol whether to apply collinearity with factor congA to mode A
#' @param Bcol whether to apply collinearity with factor congB to mode B
#' @param Ccol whether to apply collinearity with factor congC to mode C
#' @param congA collinearity factor for mode A
#' @param congB collinearity factor for mode B
#' @param congC collinearity factor for mode C
#' @param eps fraction of outliers (percent contamination)
#' @param type type of outliers: one of \code{"none"} for no outliers (possible only of \code{eps==0}), 
#'  \code{"bl"} for bad leverage points, \code{"gl"} for good leverage points and
#'  \code{"og"} for orthogonal outliers
#' @param c1 parameter for outlier generation (\code{c1=10} for \code{type="bl"} 
#'  or  \code{type="gl"} and \code{c1=1} for \code{type="og"})
#' @param c2 parameter for outlier generation (\code{c2=0.1} for \code{type="bl"}
#'  or \code{type="og"} and \code{c2=0} for \code{type="gl"})
#' @param silent whether to issue warnings
#' @return A list consisting of the following lists:
#'  \itemize{
#'  \item As list of \code{nsim} matrices for the mode A
#'  \item Bs list of \code{nsim} matrices for the mode B
#'  \item Cs list of \code{nsim} matrices for the mode C
#'  \item Xs list of \code{nsim} PARAFAC data sets, each with dimension IxJxK
#'  \item Os list of \code{nsim} vectors containing the added outliers (if any)
#'  \item param list of parameters used for generation of the data sets
#'  }
#'
#' @references
#' Todorov, V. and Simonacci, V. and Gallo, M. and Trendafilov, N. (2023). A novel 
#'  estimation procedure for robust CANDECOMP/PARAFAC model fitting. 
#'  \emph{Econometrics and Statistics}. In press.
#'
#' Tomasi, G. and Bro, R., (2006). A comparison of algorithms for fitting the PARAFAC model. 
#'  \emph{Computational Statistics & Data Analysis} \bold{50} (7), 1700--1734.
#'
#' Faber, N.M. and Bro, R. and Hopke, P.K. (2003). Recent developments in CANDECOMP/PARAFAC algorithms: 
#'  A critical review. \emph{Chemometrics and Intelligent Laboratory Systems} \bold{65}, 119--137.
#'
#' @examples
#'  ##  Generate one PARAFAC data set (nsim=1) with R=2 components (nf=2) and dimensions
#'  ##  50x10x10. Apply 0.15 homoscedastic noise and 0.10 heteroscedastic noise, apply 
#'  ##  collinearity with congruence factor 0.5 to all modes. Add 20% bad leverage points.
#'
#'  library(rrcov3way)
#'  xdat <- cp_gen(I=50, J=100, K=10, nsim=1, nf=2,
#'      noise=0.15, noise1=0.10, Acol=TRUE, Bcol=TRUE, Ccol=TRUE,
#'      congA=0.5, congB=0.5, congC=0.5,
#'      eps=0.2, type="bl")
#'  names(xdat)
#'    
cp_gen <- function(I=20, J=20, K=20, nsim=200, nf=3,
    noise=0.05, noise1=0, 
    Acol=TRUE, Bcol=TRUE, Ccol=TRUE,
    congA=0.5, congB=0.5, congC=0.5,
    eps=0, type=c("none", "bl", "gl", "og"), c1=10, c2=0.1,
    silent=FALSE) {

    ## Add factor collinearity c to the matrix A
    addcong <- function(A, c) {
        nf <- ncol(A)
        congr <- matrix(c, nf, nf)          # congruence matrix we want, e.g. t(A) %*% (A)
        diag(congr) <- 1
        R <- chol(congr)                    # upper triangular matrix

        A.qr <- qr(A)
        Q <- qr.Q(A.qr)
        R1 <- qr.R(A.qr)
        sgn <- sign(diag(R1))               # assicura elementi positivi
        R.new <- diag(sgn) %*% R            # assicura elementi positivi
        A <- Q %*% R.new
        A
    }

    type <- match.arg(type)         # type of outliers
    if(eps == 0 && type != "none")
        stop("Incorrect type of outliers specified for eps==0. Must be type='none'!")
    if(eps != 0 && type == "none")
        stop("Type of outliers not specified for eps>0. Must be one of 'bl', 'gl' or 'og'!")
        
    if(missing(c1) | missing(c2)) {
        if(type == "bl") {
            c1 <- 10; c2 <- 0.1
        } else if(type == "gl") {
            c1 <- 10; c2 <- 0.0
        } else if(type == "og") {
            c1 <- 1; c2 <- 0.1
        } else
            stop("Undefined outlier type")
    }

    param <- list(I=I, J=J, K=K, nsim=nsim, nf=nf, noise=noise, noise1=noise1,
        Acol=Acol, Bcol=Bcol, Ccol=Ccol, congA=congA, congB=congB, congC=congC,
        eps=eps, type=type, c1=c1, c2=c2)

    Xmat = array(NA, c(I, J, K))
    Xlist <- Alist <- Blist <- Clist <- Olist <- vector("list", nsim)

    for(i in 1:nsim) {
        Amat <- matrix(runif(I*nf), I, nf)
        Bmat <- matrix(runif(J*nf), J, nf)
        Cmat <- matrix(runif(K*nf), K, nf)

        ## Add factor collinearity as requested
        if(Acol)
            Amat <- addcong(Amat, congA)
        if(Bcol)
            Bmat <- addcong(Bmat, congB)
        if(Ccol)
            Cmat <- addcong(Cmat, congC)

        ##  Generate the "pure" (i.e. no noise) X array
        Xmat <- toArray(Amat %*% t(krp(Cmat, Bmat)), I, J, K)

        ## Generate homoscedastic noise
        # the coefficient sqrt(n/(1-n)) adds appropriate level of noise in terms of total variability
        # e.g. 5% -> 0.2294, 10% -> 0.3333, 20% -> 0.5
        E <- array(rnorm(I*J*K), c(I, J, K))
        for(k in 1:K) {
            E[,, k] <- sqrt(noise/(1-noise)) * E[,, k] %*% diag(sqrt(diag(t(Xmat[,,k]) %*% Xmat[,,k])) / sqrt(diag(t(E[,,k]) %*% E[,,k])))
        }

        ## Generate heteroschedastic noise
        E1 <- array(rnorm(I*J*K), c(I, J, K))
        for(k in 1:K) {
            E1[,,k] <- E1[,,k] * Xmat[,, k]
            E1[,,k] <- sqrt(noise1/(1-noise1)) * E1[,,k] %*% diag(sqrt(diag(t(Xmat[,,k]) %*% Xmat[,,k]))/sqrt(diag(t(E1[,,k]) %*% E1[,,k])))
        }

        ## Add outliers
        iout <- c()
        if(eps > 0) {
            iout <- sample(1:I, size=eps*I)

            if(type == "og") {
                if(c1 != 1 && !silent)
                    warning("Bad leverage points are generated instead of residual outliers!")
                Xmat <- Xmat + E
                Xmat <- Xmat + E1
                Xmat[iout,,] <- Xmat[iout,,] + c2
            } else if(type == "gl") {
                if(c2 != 0 && !silent)
                    warning("Bad leverage points are generated instead of good leverage points")
                Xmat[iout,,] <- Xmat[iout,,] * c1
                Xmat <- Xmat + E
                Xmat <- Xmat + E1
            } else if(type == "bl") {
                Xmat[iout,,] <- Xmat[iout,,] * c1 + c2
                Xmat <- Xmat + E
                Xmat <- Xmat + E1
            }
        } else {
            Xmat <- Xmat + E
            Xmat <- Xmat + E1
        }

        ## Prepare the output: A, B, C, X and O[utliers]
        Alist[[i]] <- Amat
        Blist[[i]] <- Bmat
        Clist[[i]] <- Cmat

        Xlist[[i]] <- Xmat

        if(length(iout) > 0)
            Olist[[i]] <- iout
    }

    list(As=Alist, Bs=Blist, Cs=Clist, Xs=Xlist, Os=Olist, param=param)
}
