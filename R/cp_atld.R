##  VT::08.03.2023
##
##  roxygen2::roxygenise("c:/Users/valen/OneDrive/MyRepo/R/rrcov3way", load_code=roxygen2:::load_installed)
##
#'  Alternating Trilinear Decomposition (ATLD) for Candecomp/Parafac (CP)
#'
#' @description Alternating Trilinear Decomposition algorithm for estimating
#'  the Candecomp/Parafac (CP) model. It is based on an alternating least squares
#'  principle and replaces the iterative procedure used in the traditional
#'  PARAFAC algorithm by an improved procedure. This improved procedure
#'  contains Moore-Penrose pseudoinverse computations based on singular
#'  value decomposition (SVD) which should be theoretically more robust
#'  to similarities in spectra and time profiles
#'
#' @param X A three-way array or a matrix. If \code{X} is a matrix
#'  (matricised threeway array), \code{n}, \code{m} and \code{p} must be
#'  given and are the number of A-, B- and C-mode entities respectively
#'
#' @param n Number of A-mode entities
#' @param m Number of B-mode entities
#' @param p Number of C-mode entities
#' @param ncomp Number of components to extract
#' @param start Initial values for the A, B and C components. Can be
#'  \code{"svd"} for starting point of the algorithm from SVD's,
#'  \code{"random"} for random starting point (orthonormalized
#'  component matrices or nonnegative matrices in case of nonnegativity
#'  constraint), or a list containing user specified components.
#' @param conv Convergence criterion, default is \code{conv=1e-6}.
#' @param maxit Maximum number of iterations, default is \code{maxit=10000}.
#' @param trace Logical, provide trace output.
#' @return The result of the decomposition as a list with the following
#'  elements:
#'    \itemize{
#'    \item \code{fit} Value of the loss function
#'    \item \code{fp} Fit value expressed as a percentage
#'    \item \code{ss} Sum of squares
#'    \item \code{A} Component matrix for the A-mode
#'    \item \code{B} Component matrix for the B-mode
#'    \item \code{C} Component matrix for the C-mode
#'    \item \code{iter} Number of iterations
#'    \item \code{tripcos} Minimal triple cosine between two components
#'        across the three component matrices, used to inspect degeneracy
#'    \item \code{mintripcos} Minimal triple cosine during the iterative
#'        algorithm observed at every 10 iterations, used to inspect degeneracy
#'    \item \code{ftiter} Matrix containing in each row the function value
#'        and the minimal triple cosine at every 10 iterations
#'    }
#'
#' @references
#'  H.-L. Wu, M. Shibukawa, K. Oguma, An alternating trilinear decomposition
#'  algorithm with application to calibration of HPLC-DAD for
#'  simultaneous determination of overlapped chlorinated aromatic hydrocarbons,
#'  \emph{Journal of Chemometrics} \bold{12} (1998) 1--26.
#'
#' @examples
#'
#' \dontrun{
#' ## Example with the OECD data
#'  data(elind)
#'  dim(elind)
#'
#'  res <- cp_atld(elind, ncomp=3)
#'  res$fp
#'  res$fp
#'  res$iter
#' }
#' @export
#' @author Valentin Todorov, \email{valentin.todorov@@chello.at}; Violetta Simonacci, \email{violetta.simonacci@unina.it}
#'
cp_atld <- function(X, n, m, p, ncomp, conv=1e-06, start="random", maxit=5000, trace=FALSE)
{
    if(missing(ncomp))
        stop("Number of factors to extract 'ncomp' must be provided!")
    r <- ncomp

    if(length(dim(X)) == 3)
    {
        dn <- dimnames(X)

        Xa <- unfold(X)
        n <- dim(X)[1]
        m <- dim(X)[2]
        p <- dim(X)[3]

        Xbb <- matrix(aperm(X, perm=c(2,3,1)), m, n*p)          # unfold(X, mode="B")
        Xc <- matrix(aperm(X, perm=c(3,1,2)), p, n*m)           # unfold(X, mode="C")
    }else if(length(dim(X)) == 2)
    {
        Xa <- as.matrix(X)
        if(missing(n) | missing(m) | missing(p))
        stop("The three dimensions of the matricisized array must be provided!")
        if(n != dim(Xa)[1])
        stop("'n' must be equal to the first dimension of the matrix 'X'!")
        if(m*p != dim(Xa)[2])
        stop("'m*p' must be equal to the second dimension of the matrix 'X'!")

        data <- toArray(Xa, n, m, p)
        Xbb <- matrix(aperm(data, perm=c(2,3,1)), m, n*p)          # unfold(X, mode="B")
        Xc <- matrix(aperm(data, perm=c(3,1,2)), p, n*m)        # unfold(X, mode="C")
        rm(data)
    }else
    stop("'X' must be three dimensional array or a matrix!")

    ## Check initial values
    if(!is.list(start) && length(start) != 1)
        stop("'start' must be either a list with elements A, B and C or a single character - one of 'random' or 'svd'!")

    if(!is.list(start) && !(start %in% c("random", "svd")))
        stop("'start' must be either a list with elements A, B and C or one of 'random' or 'svd'!")

    single_iter <- matrix(0, maxit)
    ftiter <- matrix(0, maxit/10, 2)
    mintripcos <- 0
    epsilon <- 10 * .Machine$double.eps * max(n, m, p)

        ssx <- sum(X^2)

        ## If start == "random" OR start == "svd" and n or m or p < r
        ##  - not supported:    OR const[i]=="nonneg"
        ## Random start (with orthonormalized component matrices),
        ##  -- not supported: nonnegative if const == "nonneg"
        A <- pracma::orth(matrix(rnorm(max(n,r) * r), max(n,r)))[1:n, , drop=FALSE]
        B <- pracma::orth(matrix(rnorm(max(m,r) * r), max(m,r)))[1:m, , drop=FALSE]
        C <- pracma::orth(matrix(rnorm(max(p,r) * r), max(p,r)))[1:p, , drop=FALSE]

        if(is.list(start)) {
            A <- start$A
            B <- start$B
            C <- start$C
        }else if(start == "svd") {
            if(n >= r)
                A <- eigen(Xa %*% t(Xa))$vectors[, 1:r, drop=FALSE]
            if(m >= r)
                B <- eigen(Xbb %*% t(Xbb))$vectors[, 1:r, drop=FALSE]
            if(p >= r)
                C <- eigen(Xc %*% t(Xc))$vectors[, 1:r, drop=FALSE]
        }

        PC <- do_inv(C, epsilon)

        f <- sum((Xa - A %*% t(rrcov3way::krp(C,B)))^2)
        if(trace)
            cat(paste("\nCandecomp/Parafac function value at Start is ", f, sep = " "), fill = TRUE)

        fold <- f + 2 * conv * f
        iter <- 0

        while(abs((f - fold)/f) > conv && iter < maxit)
        {
            iter <- iter + 1
            fold <- f

            A <- A %*% diag(1/sqrt(diag(t(A) %*% A))) %*% diag(sqrt(diag(t(C) %*% C)))
            PA <- do_inv(A, epsilon)

            B <- Xbb %*% rrcov3way::krp(t(PA), t(PC))
            B <- B %*% diag(1/sqrt(diag(t(B) %*% B))) %*% diag(sqrt(diag(t(A) %*% A)))
            PB <- do_inv(B, epsilon)

            C <- Xc %*% rrcov3way::krp(t(PB),t(PA))
            C <- C %*% diag(1/sqrt(diag(t(C) %*% C))) %*% diag(sqrt(diag(t(B) %*% B)))
            PC <- do_inv(C, epsilon)

            A <- Xa %*% rrcov3way::krp(t(PC), t(PB))

            f <- sum((Xa - A %*% t(rrcov3way::krp(C,B)))^2)

            ## Record Relative Fit
            single_iter[iter] <- abs((f - fold)/f)
        }

    Rsq <- 1 - f/ssx
    fp <- 100 * Rsq

    ## Degeneracy problem if |tripcos| > 0.5
    tripcos <- min(rrcov3way::congruence(A, A) * rrcov3way::congruence(B, B) * rrcov3way::congruence(C, C))
    names(tripcos) <- c("Minimal triple cosine")
    if(iter < 10) {
        mintripcos <- tripcos
    }

    if(trace) {
        cat(paste("\nCandecomp/Parafac function value is", f, "after",
            iter, "iterations", sep = " "), fill = TRUE)
        cat(paste("Fit percentage is", round(fp, 2), "%", sep = " "), fill = TRUE)
    }

    out <- list(fit=f, fp=fp, ss=ssx, A=A, B=B, C=C, iter=iter, 
        tripcos=tripcos, mintripcos=mintripcos, single_iter=single_iter[1:iter])
    out
}
