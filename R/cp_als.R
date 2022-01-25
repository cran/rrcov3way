##
##
##  VT::20.10.2019
##
##  X   - 3-way array or a matrix. If X is a matrix (matricised 3-way array),
##      n, m and p are number of A-, B- and C-mode entities respectively
##
##  ncomp - Number of extracted components
##
##  const - Type of constraints on A, B and C respectively
##      ("none" for no constraints, "orth" for orthogonality constraints,
##      "nonneg" for nonnegativity or "zerocor" for zero correlations
##          constraints)
##
##  start: svd, random or given matrices A, B and C
##
##  roxygen2::roxygenise("c:/Users/valen/OneDrive/MyRepo/R/rrcov3way",
##      load_code=roxygen2:::load_installed, clean=TRUE)


#'  Alternating Least Squares (ALS) for Candecomp/Parafac (CP)
#'
#' @description Alternating Least Squares (ALS) algorithm with optional
#'  constraints for the minimization of the Candecomp/Parafac (CP) loss
#'  function.
#'
#' @param X A three-way array or a matrix. If \code{X} is a matrix
#'  (matricised threeway array), \code{n}, \code{m} and \code{p} must be
#'  given and are the number of A-, B- and C-mode entities respectively
#'
#' @param n Number of A-mode entities
#' @param m Number of B-mode entities
#' @param p Number of C-mode entities
#' @param ncomp Number of components to extract
#' @param const Optional constraints for each mode. Can be a three element
#'  character vector or a single character, one of \code{"none"} for no
#'  constraints (default), \code{"orth"} for orthogonality constraints,
#'  \code{"nonneg"} for nonnegativity constraints or
#'  \code{"zerocor"} for zero correlation between the extracted factors.
#'  For example, \code{const="orth"} means orthogonality constraints for
#'  all modes, while \code{const=c("orth", "none", "none")} sets the
#'  orthogonality constraint only for mode A.
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
#'    \item \code{A} Component matrix for the A-mode
#'    \item \code{B} Component matrix for the B-mode
#'    \item \code{C} Component matrix for the C-mode
#'    \item \code{f} Value of the loss function
#'    \item \code{fp} Fit value expressed as a percentage
#'    \item \code{iter} Number of iterations
#'    \item \code{tripcos} Minimal triple cosine between two components
#'        across the three component matrices, used to inspect degeneracy
#'    \item \code{mintripcos} Minimal triple cosine during the iterative
#'        algorithm observed at every 10 iterations, used to inspect degeneracy
#'    \item \code{ftiter} Matrix containing in each row the function value
#'        and the minimal triple cosine at every 10 iterations
#'    \item \code{const} Optional constraints (same as the input parameter
#'        \code{const})
#'    }
#'
#' @references
#' 	Harshman, R.A. (1970). Foundations of Parafac procedure:
#'  models and conditions for an "explanatory" multi-mode factor
#'  analysis. \emph{UCLA Working Papers in Phonetics}, 16: 1--84.
#'
#' Harshman, R. A., & Lundy, M. E. (1994). PARAFAC: Parallel factor analysis.
#'  Computational Statistics and Data Analysis, 18, 39--72.
#'
#' Lawson CL, Hanson RJ (1974). Solving Least Squares Problems.
#'  Prentice Hall, Englewood Cliffs, NJ.
#'
#' @note The argument \code{const} should be a three element character vector.
#'  Set \code{const[j]="none"} for unconstrained update in j-th mode weight
#'  matrix (the default),
#'  \code{const[j]="orth"} for orthogonal update in j-th mode weight matrix,
#'  \code{const[j]="nonneg"} for non-negative constraint on j-th mode or
#'  \code{const[j]="zerocor"} for zero correlation between the extracted
#'      factors.
#'  The default is unconstrained update for all modes.
#'
#'  The loss function to be minimized is \eqn{sum(k)|| X(k) - A D(k) B' ||^2},
#'  where \eqn{D(k)} is a diagonal matrix holding the \code{k}-th row of
#'  \code{C}.
#'
#' @examples
#'
#' \dontrun{
#' ## Example with the OECD data
#'  data(elind)
#'  dim(elind)
#'
#'  res <- cp_als(elind, ncomp=3)
#'  res$fp
#'  res$fp
#'  res$iter
#'
#'  res <- cp_als(elind, ncomp=3, const="nonneg")
#'  res$A
#' }
#' @export
#' @author Valentin Todorov, \email{valentin.todorov@@chello.at}
#'
cp_als <- function (X, n, m, p, ncomp, const="none", start="random",
    conv=1e-6, maxit=10000, trace=FALSE)
{
    if(missing(ncomp))
        stop("Number of factors to extract 'ncomp' must be provided!")
    r <- ncomp

    if(length(dim(X)) == 3)
    {
        Xa <- unfold(X)
        n <- dim(X)[1]
        m <- dim(X)[2]
        p <- dim(X)[3]

        Xb <- matrix(aperm(X, perm=c(2,1,3)), m, n*p)
        Xc <- matrix(aperm(X, perm=c(3,1,2)), p, n*m)
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
        Xb <- matrix(aperm(data, perm=c(2,1,3)), m, n*p)
        Xc <- matrix(aperm(data, perm=c(3,1,2)), p, n*m)
        rm(data)
    }else
    stop("'X' must be three dimensional array or a matrix!")

    ## Check constraints
    if(length(const) == 1)
        const <- rep(const, 3)
    if(!all(const %in% c("none", "orth", "nonneg", "zerocor")))
    stop("All elements of 'const' must be one of 'none', 'orth', 'nonneg' or 'zerocor'")

    ## If length of const is less than 3, pad it with "none"
    if(length(const) < 3)
        const <- c(const, rep("none", 3))[1:3]

    ## Check initial values
    if(!is.list(start) && length(start) != 1)
        stop("'start' must be either a list with elements A, B and C or a single character - one of 'random' or 'svd'!")

    if(!is.list(start) && !(start %in% c("random", "svd")))
        stop("'start' must be either a list with elements A, B and C or one of 'random' or 'svd'!")

    ftiter <- matrix(0, maxit/10, 2)
    mintripcos <- 0
    ssx <- sum(Xa^2)

    ## If start == "random" OR start == "svd" and n or m or p < r
    ##      OR const[i]=="nonneg"
    ## Random start (with orthonormalized component matrices),
    ##  nonnegative if const == "nonneg"
    A <- if(const[1] == "nonneg") matrix(runif(max(n,r) * r), max(n,r))[1:n, , drop=FALSE]
         else                     pracma::orth(matrix(rnorm(max(n,r) * r), max(n,r)))[1:n, , drop=FALSE]
    B <- if(const[2] == "nonneg") matrix(runif(max(m,r) * r), max(m,r))[1:m, , drop=FALSE]
         else                     pracma::orth(matrix(rnorm(max(m,r) * r), max(m,r)))[1:m, , drop=FALSE]
    C <- if(const[3] == "nonneg") matrix(runif(max(p,r) * r), max(p,r))[1:p, , drop=FALSE]
         else                     pracma::orth(matrix(rnorm(max(p,r) * r), max(p,r)))[1:p, , drop=FALSE]

    if(is.list(start)) {
        A <- start$A
        B <- start$B
        C <- start$C
    }else if(start == "svd") {
        if(n >= r && const[1] != "nonneg")
            A <- eigen(Xa %*% t(Xa))$vectors[, 1:r, drop=FALSE]
        if(m >= r && const[2] != "nonneg")
            B <- eigen(Xb %*% t(Xb))$vectors[, 1:r, drop=FALSE]
        if(p >= r && const[3] != "nonneg")
            C <- eigen(Xc %*% t(Xc))$vectors[, 1:r, drop=FALSE]
    }

    f <- sum((Xa - tcrossprod(A, krp(C, B)))^2)
    if(trace)
        cat(paste("Candecomp/Parafac function value at Start is ", f, sep = " "), fill = TRUE)

    fold <- f + 2 * conv * f
    iter <- 0
    BB <- t(B) %*% B
    CC <- t(C) %*% C

    while((fold - f > conv * f | iter < 2) & f > conv^2 & iter < maxit) {

        fold <- f

        ## Step 1: Update mode A
        if(const[1] == "none") {                # A - no constraint
            A <- Xa %*% krp(C, B) %*% solve(crossprod(C) * crossprod(B))
        } else if (const[1] == "orth")          # A - orthogonality
        {
            svd <- svd(Xa %*% krp(C, B))
            A <- svd$u %*% t(svd$v)
        } else if (const[1] == "nonneg")        # A - non-negativity
        {
            Atemp <- A
            cob <- krp(C, B)
            xtx <- crossprod(cob)
            for(i in 1:n){
                xty <- crossprod(cob, Xa[i,])
                A[i, ] <- coef(nnls(cob, Xa[i,]))
            }
            if(any(colSums(A) == 0)){
              A <- Atemp
              stop(paste("Error in nonnegative LS for mode A at iter=", iter))
            }
        } else if(const[1] == "zerocor")        # A - zero correlations constraints
        {
            XF <- Xa %*% krp(C, B)
            XFmean <- apply(XF, 2, mean)
            svd <- svd(sweep(XF, 2, XFmean))
            A <- sweep(svd$u %*% t(svd$v), 2, XFmean %*% solve(BB * CC), "+")
        }
        AA <- t(A) %*% A

        ## Step 2: Update mode B
        if (const[2] == "none") {           # B - no constraint
            B <- Xb %*% krp(C, A) %*% solve(crossprod(C) * crossprod(A))
        } else if(const[2] == "orth")       # B - orthogonality
        {
            svd <- svd(Xb %*% krp(C, A))
            B <- svd$u %*% t(svd$v)
        } else if(const[2] == "nonneg")     # B - non-negativity
        {
            Btemp <- B
            coa <- krp(C, A)
            xtx <- crossprod(coa)
            for(i in 1:m) {
                xty <- crossprod(coa, Xb[i,])
                B[i,] <- coef(nnls(coa, Xb[i,]))
            }
            if(any(colSums(B) == 0)){
              B <- Btemp
              stop(paste("Error in nonnegative LS for mode B at iter=", iter))
            }
        } else if(const[2] == "zerocor")   # B - zero correlations constraints
        {
            XF <- Xb %*% krp(C, A)
            XFmean <- apply(XF, 2, mean)
            svd <- svd(sweep(XF, 2, XFmean))
            B <- sweep(svd$u %*% t(svd$v), 2, XFmean %*% solve(AA * CC), "+")
        }
        BB <- t(B) %*% B

        ## Step 3: Update mode C
        if(const[3] == "none") {                # C - no constraint
            C <- Xc %*% krp(B, A) %*% solve(crossprod(B) * crossprod(A))
        }else if(const[3] == "orth")  {         # C - orthogonality
            svd <- svd(Xc %*% krp(B, A))
            C <- svd$u %*% t(svd$v)
        }else if(const[3] == "nonneg") {        # C - non-negativity
            Ctemp <- C
            boa <- krp(B, A)
            xtx <- crossprod(boa)
            for(i in 1:p){
                xty <- crossprod(boa, Xc[i,])
                C[i, ] <- coef(nnls(boa, Xc[i,]))
            }
            if(any(colSums(C) == 0)){
              C <- Ctemp
              stop(paste("Error in nonnegative LS for mode C at iter=", iter))
            }
        }else if(const[3] == "zerocor") {       # C - zero correlation constraint
            XF <- Xc %*% krp(B, A)
            XFmean <- apply(XF, 2, mean)
            svd <- svd(sweep(XF, 2, XFmean))
            C <- sweep(svd$u %*% t(svd$v), 2, XFmean %*% solve(AA * BB), "+")
        }
        CC <- t(C) %*% C

        ## Step 4: Convergence check
        f <- sum((Xa - tcrossprod(A, krp(C, B)))^2)
        iter <- iter + 1

        if(iter %% 10 == 0)
        {
            tripcos <- min(congruence(A, A) * congruence(B, B) * congruence(C, C))
            if (iter == 10)
                mintripcos <- tripcos
            if (tripcos < mintripcos)
                mintripcos <-  tripcos
            if ((iter %% 1000) == 0 & trace)
                cat(paste("Minimal Triple cosine =", tripcos), fill = TRUE)

            ftiter[iter/10, ] <- c(f, tripcos)
        }

        if(iter %% 50 == 0 && trace) {
            cat(paste("f=", f, "after", iter, "iters; diff.=", (fold - f),
                sep = " "), fill = TRUE)
        }
    }

    ftiter <- ftiter[1:iter/10, ]           # take only the first iter/10 rows
    Rsq <- 1 - f/ssx
    fp <- 100 * Rsq

    tripcos <- min(congruence(A, A) * congruence(B, B) * congruence(C, C))
    names(tripcos) <- c("Minimal triple cosine")

    if(iter < 10)
        mintripcos <- tripcos

    if(trace) {
        cat(paste("\nCandecomp/Parafac function value is", f, "after",
            iter, "iterations", sep = " "), fill = TRUE)
        cat(paste("Fit percentage is", round(fp, 2), "%", sep = " "), fill = TRUE)
    }

    out <- list(A=A, B=B, C=C, f=f, fp=fp, ss=ssx, iter=iter, tripcos=tripcos,
            mintripcos=mintripcos, ftiter=ftiter, const=const)
    out
}
