######
##  VT::16.09.2019
##
##  - mtrace:       The trace of a square numeric matrix
##  - krp:          The Khatri-Rao product of two matrices
##  - congruence:   Tuker's congruence coefficient
##
##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/R/rrcov3way", load_code=roxygen2:::load_installed)
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
#'  require(rrcov)
#'
#'  data(delivery, package="robustbase")
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

##  Find the general inverse of a (rectangular) matrix a, i.e. find
##  a matrix G such that AGA == A.
##  First try to find   t(A %*% solve(crossprod(A))) and if it does not work,
##  calculate the general inverse using the function pinv() from package pracma.
##
do_inv <- function(a, tolerance=1e-12) {
    tryCatch(expr={t(a %*% solve(crossprod(a)))},
             error=function(msg) pracma::pinv(a, tol=tolerance))
}
