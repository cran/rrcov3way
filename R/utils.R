######
##  VT::16.09.2019
##
##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/R/rrcov3way", load_code=roxygen2:::load_installed)
##

do3Scale <- function (x, ...) UseMethod("do3Scale")

#'  Varimax Rotation for Tucker3 models
#'
#' @description Computes \emph{varimax} rotation of the core and component matrix of a Tucker3 model to simple structure.
#'
#' @param x A Tucker 3 object
#' @param \dots Potential further arguments passed to called functions.
#' @return A list including the following components:
#'
#' @examples
#'  ## Rotation of a Tucker3 solution
#'  data(elind)
#'  (t3 <- Tucker3(elind, 3, 2, 2))
#'  xout <- do3Rotate(t3, c(3, 3, 3), rotate=c("A", "B", "C"))
#'  xout$vvalue
#'
#' @rdname do3Rotate
#' @export
#' @author Valentin Todorov, \email{valentin.todorov@@chello.at}
do3Rotate <- function (x, ...) UseMethod("do3Rotate")

do3Postprocess <- function (x, ...) UseMethod("do3Postprocess")
reflect <- function (x, ...) UseMethod("reflect")
reorder <- function (x, ...) UseMethod("reorder")
coordinates <- function (x, ...) UseMethod("coordinates")

permute <- function(X, n, m, p) {
    matrix(as.vector(t(as.matrix(X))), m, n*p)
}

toArray <- function(x, n, m, r, mode=c("A", "B", "C")) {
    xa2arr<-function(x, n, m, p)
    {
        x <- as.matrix(x)
        X <- array(0, c(n, m, p))
        for(k in 1:p)
            for(j in 1:m)
                X[, j, k] <- x[, (k - 1) * m + j]
        X
    }

    mode <- match.arg(mode)
    ret <- if(mode == "A") xa2arr(x, n=n, m=m, p=r)
           else if(mode == "B") xa2arr(permute(permute(x, m, r, n), r, n, m), n, m, r)
           else xa2arr(permute(x, r, n, m), n, m, r)
    ret
}

unfold <- function(x, mode=c("A", "B", "C")) {
    mode <- match.arg(mode)
    if(length(dm <- dim(x)) < 3)
        stop("\n'X' is not a multyway array!")
    n <- dm[1]
    m <- dm[2]
    p <- dm[3]

    switch(mode,
        A = {
            mat <- matrix(NA, nrow = n, ncol = m * p)
            for(k in 1:p)
                mat[, ((k - 1) * m + 1):(k * m)] <- x[, , k]
            mat
        },
        B = {
            mat <- matrix(NA, nrow = m, ncol = p * n)
            for(i in 1:n)
                mat[, ((i - 1) * p + 1):(i * p)] <- x[i, , ]
            mat
        },
        C = {
            mat <- matrix(NA, nrow = p, ncol = n * m)
            for(j in 1:m)
                mat[, ((j - 1) * n + 1):(j * n)] <- t(x[, j, ])
            mat
        }
    )
}

##  Represent a 3-dimensional array as a tall matrix
##
tallArray <- function(X)
    apply(X, 2, rbind)

##  Represent a 3-dimensional array as a wide matrix (unfold in mode A)
##
wideArray <- function(X)
    unfold(X)

##  Turn a tall IKxJ matrix into an wide IxJK matrix
##
tall2wide <- function(Xtall, I, J, K)
{
    Xwide <- NULL
    for (i in 1:K)
        Xwide <- cbind(Xwide, Xtall[((i-1)*I + 1):(i*I) ,])
    Xwide
}

## ILR transformation (see package 'chemometrics')
.ilrV <- function(x)
{
    dn <- dimnames(x)
    x.ilr <- matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
    if(!is.null(dn) && is.list(dn))
        rownames(x.ilr) <- dn[[1]]
    colnames(x.ilr) <- paste0("Z", 1:(ncol(x) - 1))
    for (i in seq_len(ncol(x.ilr)))
    {
        x.ilr[, i] <- sqrt((i)/(i + 1)) * log(((apply(as.matrix(x[,1:i]), 1, prod))^(1/i))/(x[, i + 1]))
    }

    if (is.data.frame(x))
        x.ilr <- data.frame(x.ilr)
    return(x.ilr)
}

.clrV <- function(x)
{
    res <- if(dim(x)[2] == 1) x
           else{
            	gm <- apply(x, 1, .gm)
            	res <- log(x/gm, exp(1))
            }
    res
}

## Geometric mean
.gm <- function (x)
{
    if (!is.numeric(x))
        stop("x has to be a vector of class numeric")
    if(any(na.omit(x == 0))) 0 else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
}

ilrArray <- function(x) {
    di <- dim(x)
    I <- di[1]
    J <- di[2]
    K <- di[3]

    dn <- dimnames(x)
    Xtall <- tallArray(x)
    Xilr <- .ilrV(Xtall)
    Xwideilr <- tall2wide(Xilr, I, J, K)

    ret <- toArray(Xwideilr, I, J-1, K)
    dimnames(ret)[[1]] <- dn[[1]]
    dimnames(ret)[[2]] <- paste0("Coord-", 1:(J-1))
    dimnames(ret)[[3]] <- dn[[3]]
    ret
}

clrArray <- function(x) {
    di <- dim(x)
    I <- di[1]
    J <- di[2]
    K <- di[3]

    dn <- dimnames(x)
    Xtall <- tallArray(x)
    Xclr <- .clrV(Xtall)
    Xwideclr <- tall2wide(Xclr, I, J, K)

    ret <- toArray(Xwideclr, I, J, K)
    dimnames(ret) <- dn
    ret
}

## Check if the matrix a is with orthonormal columns:
## - a matrix has orthonormal columns if its Gram matrix is the identity
#   i.e. A'A == I
##
is.orthogonal <- function(a)
{
    ata <- t(a) %*% a
    diag(ata) <- 1
    ret <- all.equal(ata, diag(ncol(a)), check.attributes=FALSE)
    is.logical(ret) && ret
}

is.orthonormal <- function(a)
{
    ret <- all.equal(t(a) %*% a, diag(ncol(a)), check.attributes=FALSE)
    is.logical(ret) && ret
}

do3Scale.tucker3 <- function(x, renorm.mode=c("A", "B", "C"), ...)
{
    renorm.mode <- match.arg(renorm.mode)
    P <- ncol(x$A)
    Q <- ncol(x$B)
    R <- ncol(x$C)

    if(renorm.mode == "A")
    {
        ss <- diag(sqrt(rowSums(x$GA^2)), P)
        dimnames(ss) <- list(colnames(x$A), colnames(x$A))
    	x$A <- x$A %*% ss
    	x$GA <- solve(ss) %*% x$GA
    } else if(renorm.mode == "B")
    {
    	x$GA <- permute(x$GA, P, Q, R)
        ss <- diag(sqrt(rowSums(x$GA^2)), Q)
        dimnames(ss) <- list(colnames(x$B), colnames(x$B))
    	x$B <- x$B %*% ss
    	x$GA <- solve(ss) %*% x$GA
    	x$GA <- permute(x$GA, Q, R, P)
    	x$GA <- permute(x$GA, R, P, Q)
    } else
    {
    	x$GA <- permute(x$GA, P, Q, R)
    	x$GA <- permute(x$GA, Q, R, P)
        ss <- diag(sqrt(rowSums(x$GA^2)), R)
        dimnames(ss) <- list(colnames(x$C), colnames(x$C))
    	x$C <- x$C %*% ss
    	x$GA <- solve(ss) %*% x$GA
    	x$GA <- permute(x$GA, R, P, Q)
    }
    x
}

do3Scale.parafac <- function(x, renorm.mode=c("A", "B", "C"), ...)
{
    renorm.mode <- match.arg(renorm.mode)
    R <- ncol(x$A)

    ssa <- diag(1/sqrt(colSums(x$A^2)), R)
    ssb <- diag(1/sqrt(colSums(x$B^2)), R)
    ssc <- diag(1/sqrt(colSums(x$C^2)), R)
    dimnames(ssa) <- list(colnames(x$A), colnames(x$A))
    dimnames(ssb) <- list(colnames(x$B), colnames(x$B))
    dimnames(ssc) <- list(colnames(x$C), colnames(x$C))
    if(renorm.mode == "A")
    {
    	x$B <- x$B %*% ssb
    	x$C <- x$C %*% ssc
    	x$A <- x$A %*% solve(ssb) %*% solve(ssc)

    } else if(renorm.mode == "B")
    {
    	x$A <- x$A %*% ssa
    	x$C <- x$C %*% ssc
    	x$B <- x$B %*% solve(ssa) %*% solve(ssc)

    } else
    {
    	x$A <- x$A %*% ssa
    	x$B <- x$B %*% ssb
    	x$C <- x$C %*% solve(ssa) %*% solve(ssb)
    }

    x
}

do3Scale.default <- function(x, center=FALSE, scale=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale.mode=c("B", "A", "C"), only.data=TRUE, ...)
{
    ss <- function(x) sqrt(sum(x^2))

    center.mode <- match.arg(center.mode)
    cmode <- vector("character", length=3)
    for(i in 1:nchar(center.mode))
        cmode[i] <- substr(center.mode, i, i)
    scale.mode <- match.arg(scale.mode)
    di <- dim(x)
    dn <- dimnames(x)

    stopifnot(length(di) == 3L)
    xcenter <- NULL
    xscale <- NULL

    if(!is.logical(center) || center == TRUE)
    {
        for(i in seq_len(length(cmode)))
        {
            cmi <- cmode[i]
            if(nchar(cmi) > 0)
            {
                ## print(paste("Centering mode ", cmi))
                x <- unfold(x, cmi)
                if(is.logical(center) && center == TRUE)
                    center <- mean

                x1 <- robustbase::doScale(x, center=center, scale=NULL)
                xcenter <- x1$center
                x <- toArray(x1$x, di[1], di[2], di[3], mode=cmi)
                dimnames(x) <- dn
            }
        }
    }

    if(!is.null(scale) && (!is.logical(scale) || scale == TRUE))
    {
        x <- t(unfold(x, scale.mode))
        if(is.logical(scale) && scale == TRUE)
            scale <- ss

        x1 <- robustbase::doScale(x, center=NULL, scale=scale)
        xscale <- x1$scale
        x <- toArray(t(x1$x), di[1], di[2], di[3], mode=scale.mode)
        dimnames(x) <- dn
    }

    ret <- if(only.data) x else list(x=x, center=xcenter, center.mode=center.mode, scale=xscale, scale.mode=scale.mode)
    ret
}

#' @param weights A numeric vector with length 3: relative weights (greater or
#'  equal 0) for the simplicity of the component matrices \code{A}, \code{B}
#'  and \code{C} respectively.
#' @param rotate Within which mode to rotate the Tucker3 solution:
#'  \code{rotate="A"} means to rotate the component matrix \code{A} of mode A;
#'  \code{rotate=c("A", "B")} means to rotate the component matrices \code{A}
#'  and \code{B} of modes A and B respectively. Default is to rotate all modes,
#'  i.e. \code{rotate=c("A", "B", "C")}.
#'
#' @rdname do3Rotate
#' @export
do3Rotate.tucker3 <- function(x, weights=c(0, 0, 0), rotate=c("A", "B", "C"), ...)
{
    rot1 <- rot2 <- rot3 <- 0
    rot1 <- ifelse("A" %in% rotate, 1, rot1)
    rot2 <- ifelse("B" %in% rotate, 1, rot2)
    rot3 <- ifelse("C" %in% rotate, 1, rot3)

    rot <- ThreeWay::varimcoco(x$A, x$B, x$C, x$G, weights[1], weights[2], weights[3], rot1, rot2, rot3, nanal=1)
    dn1 <- dimnames(x$A)
    dn2 <- dimnames(x$B)
    dn3 <- dimnames(x$C)
    dn <- dimnames(x$GA)
    x$A <- rot$AS
    dimnames(x$A) <- dn1
    x$B <- rot$BT
    dimnames(x$B) <- dn2
    x$C <- rot$CU
    dimnames(x$C) <- dn3
    x$GA <- rot$K
    dimnames(x$GA) <- dn

    vvalue <- rep(0, 4)
    names(vvalue) <- c("GA", "A", "B", "C")
    vvalue[1] <- rot$f1
    vvalue[2] <- rot$f2a
    vvalue[3] <- rot$f2b
    vvalue[4] <- rot$f2c

    ans <- list(x=x, S=rot$S, T=rot$T, U=rot$U, vvalue=vvalue, f=rot$f)
    class(ans) <- "rotation"

    ans
}

coordinates.parafac <- function(x, mode=c("A", "B", "C"), type=c("normalized", "unit", "principal"), ...)
{
    mode <- match.arg(mode)
    type <- match.arg(type)
    a <- switch(mode,
            "A"=x$A,
            "B"=x$B,
            "C"=x$C)
    ncomp <- ncol(x$A)
    nelem <- nrow(a)

    eq <- all.equal(sum(colSums(a^2)), ncomp)
    if(!is.logical(eq) || !eq)
        stop("Solution is not normalized for mode ", mode)

    switch(type,
        "normalized"=       # normalized coordinates - what is returned by renormalize()
                a,
        "unit"=             # Unit mean-square coordinates
                a*sqrt(nelem),
        "principal"=        # Principal coordinates (loadings)
        {
            pc <- a %*% diag(sqrt(weights(x) * nelem))
            colnames(pc) <- colnames(a)                 # diag() loses the names
            pc
        })
}

coordinates.tucker3 <- function(x, mode=c("A", "B", "C"), type=c("normalized", "unit", "principal"), ...)
{
    mode <- match.arg(mode)
    type <- match.arg(type)
    a <- switch(mode,
            "A"=x$A,
            "B"=x$B,
            "C"=x$C)
    ncomp <- ncol(a)
    nelem <- nrow(a)

    switch(type,
        "normalized"=       # normalized coordinates - what is returned by renormalize()
                a,
        "unit"=             # Unit mean-square coordinates
                a*sqrt(nelem),
        "principal"=        # Principal coordinates (loadings)
        {
            pc <- a %*% diag(sqrt(weights(x) * nelem))
            colnames(pc) <- colnames(a)                 # diag() loses the names
            pc
        })
}

reflect.tucker3 <- function(x, mode=c("A", "B", "C"), rsign=1, ...)
{
    mode <- match.arg(mode)
    P <- ncol(x$A)
    Q <- ncol(x$B)
    R <- ncol(x$C)
    ncomp <- if(mode=="A") P else if(mode=="B") Q else R
    xnames <- if(mode=="A") dimnames(x$A) else if(mode=="B") dimnames(x$B) else dimnames(x$C)
    GAnames <- dimnames(x$GA)

    stopifnot(length(rsign) == 1 | length(rsign) == ncomp)

    if(length(rsign) == 1)
        rsign <- rep(rsign, ncomp)
    rsign <- sign(rsign)
    rdiag <- diag(rsign, nrow=ncomp)

    if(mode == "A")
    {
        x$A <- x$A %*% rdiag
        x$GA <- rdiag %*% x$GA
        dimnames(x$A) <- xnames
        dimnames(x$GA) <- GAnames
    } else if(mode=="B")
    {
        x$B <- x$B %*% rdiag
        x$GA <- permute(x$GA, P, Q, R)
        x$GA <- rdiag %*% x$GA
        x$GA <- permute(x$GA, Q, R, P)
        x$GA <- permute(x$GA, R, P, Q)
        dimnames(x$B) <- xnames
        dimnames(x$GA) <- GAnames
    }else
    {
        x$C <- x$C %*% rdiag
        x$GA <- permute(x$GA, P, Q, R)
        x$GA <- permute(x$GA, Q, R, P)
        x$GA <- rdiag %*% x$GA
        x$GA <- permute(x$GA, R, P, Q)
        dimnames(x$C) <- xnames
        dimnames(x$GA) <- GAnames
    }

    x
}

reflect.parafac <- function(x, mode=c("A", "B", "C"), rsign=1, ...)
{
    mode <- match.arg(mode)
    ncomp <- x$ncomp
    anames <- dimnames(x$A)
    bnames <- dimnames(x$B)
    cnames <- dimnames(x$C)

    stopifnot(length(rsign) == 1 | length(rsign) == ncomp)

    if(length(rsign) == 1)
        rsign <- rep(rsign, ncomp)
    rsign <- sign(rsign)
    rdiag <- diag(rsign, nrow=ncomp)

    if(mode == "A")
    {
        x$A <- x$A %*% rdiag
        x$B <- x$B %*% rdiag
        dimnames(x$A) <- anames
        dimnames(x$B) <- bnames
    } else if(mode=="B")
    {
        x$B <- x$B %*% rdiag
        x$C <- x$C %*% rdiag
        dimnames(x$B) <- bnames
        dimnames(x$C) <- cnames
    }else
    {
        x$C <- x$C %*% rdiag
        x$A <- x$A %*% rdiag
        dimnames(x$C) <- cnames
        dimnames(x$A) <- anames
    }

    x
}

##  Calculate the explained variance (standardized weights)
##  by component for a PARAFAC model.
##
##  Note that the stdandard weights can be extracted only if at least one of the
##  modes has orthogonal components (i.e. the model was
##  estimated with orthogonality constraint).
weights.parafac <- function(object, ...)
{
    if(!is.orthogonal(object$A) && !is.orthogonal(object$B) && !is.orthogonal(object$C))
        warning(paste0("It is not possible to obtain a partitioning of the ",
                       "total variability by components since none of the ",
                       "modes has orthogonal components!"))

    object$GA^2/object$ss
}

##  Calculate the explained variance (standardized weights)
##  by component for a TUCKER 3 model.
##
weights.tucker3 <- function(object, mode=c("A", "B", "C"), ...)
{
    mode <- match.arg(mode)
    P <- ncol(object$A)
    Q <- ncol(object$B)
    R <- ncol(object$C)

    if(mode=="A")
        ret <- diag(object$GA %*% t(object$GA))/object$ss
    else if(mode=="B")
    {
        Y <- permute(object$GA, P, Q, R)
        ret <- diag(Y %*% t(Y))/object$ss
    } else      # mode C
    {
        Y <- permute(object$GA, P, Q, R)
        Y <- permute(Y, Q, R, P)
        ret <- diag(Y %*% t(Y))/object$ss
    }
    ret
}

## Reorder the components of a PARAFAC model
##
reorder.parafac <- function(x, order=TRUE, ...)
{
    ncomp <- x$ncomp
    if(is.logical(order))
    {
        ## Order the components in decreasing order of the
        ##  explained variance (standardized weights). Note that the
        ##  standardized weights can be extracted only if at least one of the
        ##  modes has orthonormal components (i.e. the model was
        ##  estimated with orthogonality constraint).
        if(order)
            order <- order(weights(x), decreasing=TRUE)
        else
            return(x)
    } else
        order <- as.integer(order)

    if(length(unique(order)) != ncomp | max(order) > ncomp | min(order) < 1)
        stop("Incorrect new order specified. Must have length equal to ", ncomp, " and contain unique integers in the range 1 to ", ncomp, ".")

    x$A <- x$A[, order]
    x$B <- x$B[, order]
    x$C <- x$C[, order]
    x$GA <- x$GA[order]
    x
}

## Reorder the components of a Tucker 3 model
##
reorder.tucker3 <- function(x, mode=c("A", "B", "C"), order=TRUE, ...)
{
    mode <- match.arg(mode)
    P <- ncol(x$A)
    Q <- ncol(x$B)
    R <- ncol(x$C)
    ncomp <- if(mode=="A") P else if(mode=="B") Q else R

    if(is.logical(order))
    {
        ## Order the components in decreasing order of the
        ##  explained variance (standardized weights).
        if(order)
            order <- order(weights(x, mode=mode), decreasing=TRUE)
        else
            return(x)
    } else
        order <- as.integer(order)

    order <- as.integer(order)
    if(length(unique(order)) != ncomp | max(order) > ncomp | min(order) < 1)
        stop("Incorrect new order specified. Must have length equal to ", ncomp, " and contain unique integers in the range 1 to ", ncomp, ".")

    dn <- dimnames(x$GA)
    if(mode == "A")
    {
        x$A <- x$A[, order]
        ga <- toArray(x$GA, P, Q, R)
        ga <- ga[order,,]
        x$GA <- unfold(ga)
        rownames(x$GA) <- colnames(x$A)
        colnames(x$GA) <- dn[[2]]
    } else if(mode=="B")
    {
        x$B <- x$B[, order]
        ga <- toArray(x$GA, P, Q, R)
        ga <- ga[,order,]
        x$GA <- unfold(ga)
        rownames(x$GA) <- dn[[1]]
    }else
    {
        x$C <- x$C[, order]
        ga <- toArray(x$GA, P, Q, R)
        ga <- ga[,,order]
        x$GA <- unfold(ga)
        rownames(x$GA) <- dn[[1]]
    }

    x
}

do3Postprocess.tucker3 <- function(x, reflectA, reflectB, reflectC, reorderA, reorderB, reorderC, ...)
{
    if(!missing(reflectA))
        x <- reflect(x, mode="A", rsign=reflectA)
    if(!missing(reflectB))
        x <- reflect(x, mode="B", rsign=reflectB)
    if(!missing(reflectC))
        x <- reflect(x, mode="C", rsign=reflectC)

    if(!missing(reorderA))
        x <- reorder(x, mode="A", order=reorderA)
    if(!missing(reorderB))
        x <- reorder(x, mode="B", order=reorderB)
    if(!missing(reorderC))
        x <- reorder(x, mode="C", order=reorderC)
    x
}

do3Postprocess.parafac <- function(x, reflectA, reflectB, reflectC, reorder, ...)
{
    if(!missing(reflectA))
        x <- reflect(x, mode="A", rsign=reflectA)
    if(!missing(reflectB))
        x <- reflect(x, mode="B", rsign=reflectB)
    if(!missing(reflectC))
        x <- reflect(x, mode="C", rsign=reflectC)

    if(!missing(reorder))
        x <- reorder(x, reorder)
    x
}
