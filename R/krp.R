##  Khatri-Rao product.
##
##  krp(A,B) returns the Khatri-Rao product of two matrices A and B, of
##     dimensions I-by-K and J-by-K respectively. The result is an I*J-by-K
##     matrix formed by the matching column-wise Kronecker products, i.e.
##     the k-th column of the Khatri-Rao product is defined as
##     kronecker(A[, k], B[, k]).
##

krp <- function (A, B)
{
    xdim <- dim(A)
    ydim <- dim(B)
    if(xdim[2] != ydim[2])
        stop("A and B must have same number of columns.")

    AoB <- matrix(0, ydim[1] * xdim[1], xdim[2])
    for(u in 1:xdim[2])
        AoB[, u] <- kronecker(A[, u], B[, u])

    AoB
}
