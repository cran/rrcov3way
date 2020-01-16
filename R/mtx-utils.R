######
##  VT::16.09.2019
##
##  - mtrace:       The trace of a square numeric matrix
##  - krp:          The Khatri-Rao product of two matrices
##  - congruence:   Tuker's congruence coefficient
##  - orth:
##  - orthmax2:
##  - .pinv:        Moore-Penrose pseudoinverse
##
##  roxygen2::roxygenise("C:/projects/statproj/R/rrcov3way")
##

#'  The trace of a square numeric matrix
#'
#' @description Computes the trace of a square numeric matrix. If \code{A} is not numeric and square matrix,
#'  the function terminates with an error message.
#'
#' @param A A square numeric matrix.
#'
#' @return the sum of the values on the diagonal of the matrix \code{A}, i.e. \code{sum(diag(A))}.
#'
#' @examples
#' (a <- matrix(c(5,2,3, 4,-3,7, 4,1,2), ncol=3))
#' (b <- matrix(c(1,0,1, 0,1,2, 1,0,3), ncol=3))
#'
#' mtrace(a)
#' mtrace(b)
#'
#' ## tr(A+B)=tr(A)+tr(B)
#' all.equal(mtrace(a) + mtrace(b), mtrace(a+b))
#'
#' ## tr(A)=tr(A')
#' all.equal(mtrace(a), mtrace(t(a)))
#'
#' ## tr(alphA)=alphatr(A)
#' alpha <- 0.5
#' all.equal(mtrace(alpha*a), alpha*mtrace(a))
#'
#' ##  tr(AB)=tr(BA)
#' all.equal(mtrace(a %*% b), mtrace(b %*% a))
#'
#'
#' ##  tr(A)=tr(BAB-1)
#' all.equal(mtrace(a), mtrace(b %*% a %*% solve(b)))

#' @export
#' @author Valentin Todorov, \email{valentin.todorov@@chello.at}
mtrace <- function(A)
{
    if(!is.matrix(A))
        A <- as.matrix(A)
    if(nrow(A) != ncol(A))
        stop("'A' is not a square matrix.")

    if(!is.numeric(A))
        stop("'A' is not a numeric matrix.")

    sum(diag(A))
}

#'  The Khatri-Rao product of two matrices
#'
#' @description The function \code{krp(A,B)} returns the Khatri-Rao product of two matrices \code{A} and \code{B}, of
#'  dimensions I x K and J x K respectively. The result is an IJ x K matrix formed by the matching
#'  column-wise Kronecker products, i.e. the k-th column of the Khatri-Rao product is
#'  defined as \code{kronecker(A[, k], B[, k])}.
#'
#' @param A Matrix of order I x K.
#' @param B Matrix of order J x K.
#'
#' @return The IJ x K matrix of columnwise Kronecker products.
#'
#' @references
#'      Khatri, C. G., and Rao, C. Radhakrishna (1968).
#'      Solutions to Some Functional Equations and Their Applications to Characterization of Probability Distributions.
#'      Sankhya: Indian J. Statistics, Series A 30, 167-180.
#'
#'  	Smilde, A., Bro R. and Gelardi, P. (2004). Multi-way Analysis: Applications in Chemical Sciences, Chichester:Wiley
#'
#' @examples
#' a <- matrix(1:12, 3, 4)
#' b <- diag(1:4)
#' krp(a, b)
#' krp(b, a)

#' @export
#' @author Valentin Todorov, \email{valentin.todorov@@chello.at}

krp <- function (A, B)
{
    d1 <- dim(A)
    d2 <- dim(B)
    if(d1[2] != d2[2])
        stop("A and B must have same number of columns.")

    AoB <- matrix(0, d2[1] * d1[1], d1[2])
    for(i in 1:d1[2])
        AoB[, i] <- kronecker(A[, i], B[, i])

    AoB
}

#'  Coefficient of factor congruence (phi)
#'
#' @description The function \code{congruence(x, y)} computes the Tucker's
#'  congruence (phi) coefficients among two sets of factors.
#'
#' @details Find the Tucker's coefficient of congruence between two sets of factor loadings.
#'  Factor congruences are the cosines of pairs of vectors defined by the loadings matrix
#'  and based at the origin. Thus, for loadings that differ only by a scaler
#'  (e.g. the size of the eigen value), the factor congruences will be 1.
#'
#'  For factor loading vectors of X and Y the measure of factor congruence, phi, is
#'  \deqn{
#'  \phi = \frac{\sum X Y}{\sqrt{\sum(X^2)\sum(Y^2)}}
#'  .}{phi = sum(X Y)/sqrt(sum(X^2) sum(X^2)) }
#'
#'  If \code{y=NULL} and \code{x} is a numeric matrix, the congruence
#'  coefficients between the columns of the matrix \code{x} are returned.
#'  The result is a symmetric matrix with ones on the diagonal. If two matrices
#'  are provided, they must have the same size and the result is a square matrix containing the
#'  congruence coefficients between all pairs of columns of the two matrices.
#'
#' @param x A vector or matrix of factor loadings.
#' @param y A vector or matrix of factor loadings (may be NULL).
#'
#' @return A matrix of factor congruences.
#'
#' @references
#'      L.R Tucker (1951). A method for synthesis of factor analysis studies. Personnel Research Section Report No. 984. Department of the Army, Washington, DC.
#'
#' @examples
#'  X <- getLoadings(PcaClassic(delivery))
#'  Y <- getLoadings(PcaHubert(delivery, k=3))
#'  round(congruence(X,Y),3)
#'
#' @export
#' @author Valentin Todorov, \email{valentin.todorov@@chello.at}
congruence <- function(x, y = NULL)
{
    x <- as.matrix(x)
    y <- if(is.null(y)) x else as.matrix(y)

    if(nrow(x) != nrow(y))
        stop("Both 'x' and 'y' must have the same length (number of rows)")

    dx <- if(ncol(x) == 1) matrix(1/sqrt(sum(x^2))) else diag(1/sqrt(colSums(x^2)))
    dy <- if(ncol(y) == 1) matrix(1/sqrt(sum(y^2))) else diag(1/sqrt(colSums(y^2)))

    ret <- dx %*% crossprod(x,y) %*% dy
    colnames(ret) <- colnames(y)
    rownames(ret) <- colnames(x)
    ret
}

## Moore-Penrose pseudoinverse
##
.pinv <- function(X, tol=.Machine$double.eps)
{
    X <- as.matrix(X)
    xsvd <- svd(X)
    nze <- sum( xsvd$d > (tol*xsvd$d[1]) )

    return (if(nze > 1L) xsvd$v[,1:nze] %*% diag(1/xsvd$d[1:nze]) %*% t(xsvd$u[,1:nze])
            else        outer(xsvd$v[,1],xsvd$u[,1]) / xsvd$d[1]
            )
}

#'  Orthonormal basis for the column space of matrix
#'
#' @description Computes orthonormal basis for the column space of matrix (range space, image of a matrix)
#'
#' @param A A numeric matrix.
#' @return B orthonormal basis for the column space of \code{A}.
#'
#' @details \code{orth(A)} returns an orthonormal basis for the column space of \code{A}.
#'  The columns of the result matrix \code{B} span the same space as the columns of \code{A},
#'  and the columns of \code{B} are orthogonal to each other, i.e. \code{t(B) %*% B == I}. The number of columns
#'  in \code{B} is equal to the rank of \code{A}. The orthonormal basis is obtained
#'  from \code{U} in the singular value decomposition. If \code{r = rank(A)},
#'  the first \code{r} columns of \code{U} form an orthonormal basis
#'  for the column space of \code{A}.
#'
#' @examples
#'  hilbert <- function(n) { i <- seq_len(n); 1/outer(i - 1L, i, "+") }
#'  H12 <- hilbert(12)
#'  rankMM(H12)             # -> 11 - numerically more realistic
#'  rankMM(H12, tol=0)      # -> 12
#'  B <- orth(H12)
#'
#'  t(B) %*% B
#'  ## pracma::subspace(H12, B)
#'
orth <- function (A)
{
    if(length(A) == 0)
        return(c())
    if(!is.numeric(A))
        stop("Argument 'A' must be a numeric matrix.")
    if (is.vector(A))
        A <- matrix(c(A), nrow=length(A))

    svd <- svd(A)
    U <- svd$u
    s <- svd$d
    tol <- max(dim(A)) * max(s) * .Machine$double.eps
    r <- sum(s > tol)

    U[, 1:r, drop = FALSE]
}

#'  Orthomax Rotation
#'
#' @description Performs simultaneous orthomax rotation of two matrices
#'  (using one rotation matrix).
#'
#' @details The function to be maximized is
#'  \code{sum((A1^2) - 1/nrow(A1) * gamma1 * sum((sum(A1^2))^2))^2 + sum((A2^2) - 1/nrow(A2) * gamma2 * sum((sum(A2^2))^2))^2}.
#' @param A1 A numeric matrix.
#' @param A2 A numeric matrix, with the same number of columns as \code{A1}
#' @param gamma1 orthmax parameter for A1
#' @param gamma2 orthmax parameter for A1
#' @param conv Convergence criterion (default is \code{conv=1e-6})
#' @return  A list with the following elements will be returned:
#'
#'    \itemize{
#'    \item \code{B1} rotated version of \code{A1}
#'    \item \code{B2} rotated version of \code{A2}
#'    \item \code{T} rotation matrix
#'    \item \code{f} orthomax function value
#'    }
#'

orthmax2 <- function(A1, A2, gamma1, gamma2, conv=1e-6)
{				

    m1 <- nrow(A1)
    m2 <- nrow(A2)
    r <- ncol(A1)
    if(r != ncol(A2))
    	stop("Column orders of A1 and A2 must be equal!")

    T <- diag(r)
    B1 <- A1
    B2 <- A2
    f <- sum(B1^4) - gamma1/m1 * sum((colSums(B1^2))^2) + sum(B2^4) - gamma2/m2 * sum((colSums(B2^2))^2)

    if(r>1)
    {
    	fold <- f - 2*conv*abs(f)
    	if(f == 0)
    		fold <- -conv
    	iter <- 0
	
        while(f-fold > abs(f)*conv)
        {
    		fold=f
    		iter=iter+1
    		for (i in 1:(r-1)){
    			for (j in (i+1):r){
    		
    				# Jennrich & Clarkson
    				xx=T[,i]
    				yy=T[,j]
    				a=0
    				b=0
    			
    				# for A1
    				x=B1[,i]
    				y=B1[,j]
    				x2=x^2
    				y2=y^2
    				a=a+(-gamma1/m1)*(.25*(sum(x2-y2))^2 - (sum(x*y))^2) + .25*sum(x2^2 + y2^2 - 6*x2*y2)
    				b=b+(-gamma1/m1)*sum(x*y)*sum(x2-y2) + sum((x^3)*y - x*(y^3))

    				# for A2
    				x=B2[,i]
    				y=B2[,j]
    				x2=x^2
    				y2=y^2
    				a=a+(-gamma2/m2)*(.25*(sum(x2-y2))^2 - (sum(x*y))^2) + .25*sum(x2^2 + y2^2 - 6*x2*y2)
    				b=b+(-gamma2/m2)*sum(x*y)*sum(x2-y2) + sum((x^3)*y - x*(y^3))
    				theta=0.25*atan2(b,a)
    				cs=cos(theta)
    				sn=sin(theta)
    				x=B1[,i]
    				y=B1[,j]
    				B1[,i]=cs*x+sn*y
    				B1[,j]=cs*y-sn*x
    				x=B2[,i]
    				y=B2[,j]
    				B2[,i]=cs*x+sn*y
    				B2[,j]=cs*y-sn*x
    				T[,i]=cs*xx+sn*yy
    				T[,j]=cs*yy-sn*xx
    			}
    		}	

    		f <- sum(B1^4) - gamma1/m1 * sum((colSums(B1^2))^2) + sum(B2^4) - gamma2/m2 * sum((colSums(B2^2))^2)
    	}
    }

    list(B1=B1, B2=B2, T=T, f=f)
}
