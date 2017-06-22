## Tuker's congruence coefficient between:
##  - two factos (x and y are vectors or)
##  - the columns of x (y = NULL)
##  - the columns of x and y
##
##  If x and y are matrices, they must have the same number of rows
##
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
