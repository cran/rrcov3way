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
