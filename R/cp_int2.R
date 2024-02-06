##  VT::08.03.2023
##
##  roxygen2::roxygenise("c:/Users/valen/OneDrive/MyRepo/R/rrcov3way", load_code=roxygen2:::load_installed)

#'  ATLD-ALS algorithm for Candecomp/Parafac (CP)
#'
#' @description Integrated algorithm combining ATLD and ALS for the
#'  minimization of the Candecomp/Parafac (CP) loss function.
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
#'  constraints (default), \code{"orth"} for orthogonality constraints or
#'  \code{"zerocor"} for zero correlation between the extracted factors.
#'  For example, \code{const="orth"} means orthogonality constraints for
#'  all modes, while \code{const=c("orth", "none", "none")} sets the
#'  orthogonality constraint only for mode A.
#' @param start Initial values for the A, B and C components. Can be
#'  \code{"svd"} for starting point of the algorithm from SVD's,
#'  \code{"random"} for random starting point (orthonormalized
#'  component matrices or nonnegative matrices in case of nonnegativity
#'  constraint), or a list containing user specified components.
#' @param initconv Convergence criterion for the initialization phase (ATLD),
#'  default is \code{conv=1e-2}.
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
#'    \item \code{const} Optional constraints (same as the input parameter
#'        \code{const})
#'    }
#'
#' @references
#'  H.-L. Wu, M. Shibukawa, K. Oguma, An alternating trilinear decomposition
#'  algorithm with application to calibration of HPLC-DAD for
#'  simultaneous determination of overlapped chlorinated aromatic hydrocarbons,
#'  \emph{Journal of Chemometrics} \bold{12} (1998) 1--26.
#'
#' Simonacci, V. and Gallo, M. (2020). An ATLD--ALS method for the trilinear decomposition
#'  of large third-order tensors, \emph{Soft Computing} 24 13535--13546.
#'
#' Todorov, V. and Simonacci, V. and Gallo, M. and Trendafilov, N. (2023). A novel 
#'  estimation procedure for robust CANDECOMP/PARAFAC model fitting. 
#'  \emph{Econometrics and Statistics}. In press.
#'
#' @note The argument \code{const} should be a three element character vector.
#'  Set \code{const[j]="none"} for unconstrained update in j-th mode weight
#'  matrix (the default),
#'  \code{const[j]="orth"} for orthogonal update in j-th mode weight matrix or
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
#'  res <- cp_int2(elind, ncomp=3)
#'  res$fp
#'  res$fp
#'  res$iter
#'
#'  res <- cp_int2(elind, ncomp=3, const="orth")
#'  res$A
#' }
#' @export
#' @author Valentin Todorov, \email{valentin.todorov@@chello.at}; Violetta Simonacci, \email{violetta.simonacci@unina.it}
#'
cp_int2 <- function (X, n, m, p, ncomp, initconv=1e-02, conv=1e-06,
    const="none", start="random", maxit=5000, trace=FALSE)
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

        Xb <- matrix(aperm(X, perm=c(2,1,3)), m, n*p)
        Xbb <- matrix(aperm(X, perm=c(2,3,1)), m, n*p)      # unfold(X, mode="B")
        Xc <- matrix(aperm(X, perm=c(3,1,2)), p, n*m)       # unfold(X, mode="C")
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
        Xbb <- matrix(aperm(data, perm=c(2,3,1)), m, n*p)   # unfold(X, mode="B")
        Xc <- matrix(aperm(data, perm=c(3,1,2)), p, n*m)    # unfold(X, mode="C")
        rm(data)
    }else
    stop("'X' must be three dimensional array or a matrix!")

    ## Check constraints
    if(length(const) == 1)
        const <- rep(const, 3)
    if(!all(const %in% c("none", "orth", "zerocor")))
    stop("All elements of 'const' must be one of 'none', 'orth' or 'zerocor'")

    ## If length of const is less than 3, pad it with "none"
    if(length(const) < 3)
        const <- c(const, rep("none", 3))[1:3]

    ## Check initial values
    if(!is.list(start) && length(start) != 1)
        stop("'start' must be either a list with elements A, B and C or a single character - one of 'random' or 'svd'!")

    if(!is.list(start) && !(start %in% c("random", "svd")))
        stop("'start' must be either a list with elements A, B and C or one of 'random' or 'svd'!")

    single_iter <- matrix(0, maxit)
    ftiter <- matrix(0, maxit/10, 2)
    mintripcos <- 0
    eps <- .Machine$double.eps
    epsilon <- 10 * eps * norm(Xa, "1") * max(n, m, p)

        ssx <- sum(X^2)

        ## If start == "random" OR start == "svd" and n or m or p < r
        ##  -- not supported:    OR const[i]=="nonneg"
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

        if(trace)
            cat("\n Initialization of ATLD starting ...\n")

        PC <- do_inv(C, epsilon)
        f = sum((Xa - A %*% t(krp(C, B)))^2)
        fold = f + 2 * initconv * f
        iter = 0

        if(trace)
            cat("\n Initialization of ATLD ready. Starting ATLD iteration ...\n")

        ## ATLD Iterative steps
        while(abs((f - fold)/f) > initconv && iter < maxit) {
            iter <-  iter + 1
            fold <- f

            if(trace)
                cat("\n iter=", iter, "estimating A...\n")

            A <- A %*% diag(1/sqrt(diag(t(A) %*% A))) %*% diag(sqrt(diag(t(C) %*% C)))
            PA <- do_inv(A, epsilon)

            B <- Xbb %*% krp(t(PA), t(PC))

            if(trace)
                cat("\n iter=", iter, "estimating B...\n")

            B <- B %*% diag(1/sqrt(diag(t(B) %*% B))) %*% diag(sqrt(diag(t(A) %*% A)))
            PB <- do_inv(B, epsilon)

            if(trace)
                cat("\n iter=", iter, "estimating C...\n")

            C <- Xc %*% krp(t(PB), t(PA))
            C <- C %*% diag(1/sqrt(diag(t(C) %*% C))) %*% diag(sqrt(diag(t(B) %*% B)))
            PC <- do_inv(C, epsilon)

            A <- Xa %*% krp(t(PC), t(PB))

            f <- sum((Xa - A %*% t(krp(C, B)))^2)

            ## Record Relative Fit
            single_iter[iter] <- abs((f - fold)/f)

            if(trace)
                cat("\n iter=", iter, "f=", f, abs((f - fold)/f), "\n")
        }
        iter_opt <- iter

        BB <- t(B) %*% B
        CC <- t(C) %*% C

        ## ALS iterations
        while((fold - f > conv * f | iter <= iter_opt+1) & f > conv^2 & iter < maxit) {
            iter = iter + 1
            fold = f

            ## Step 1: Update mode A
            if(const[1] == "none") {               # A - no constraint
                A <- Xa %*% krp(C, B) %*% solve(crossprod(C) * crossprod(B))
            }else if(const[1] == "orth") {         # A - orthogonality
                svd <- svd(Xa %*% krp(C, B))
                A <- svd$u %*% t(svd$v)
            }else if(const[1] == "zerocor") {      # A - zero correlation
                XF <- Xa %*% krp(C, B)
                XFmean <- apply(XF, 2, mean)
                svd <- svd(sweep(XF, 2, XFmean))
                A <- sweep(svd$u %*% t(svd$v), 2, XFmean %*% solve(BB * CC), "+")
            }
            AA = t(A) %*% A

           ## Step 2: Update mode B
            if(const[2] == "none") {              # B: no constraint
                B <- Xb %*% krp(C, A) %*% solve(crossprod(C) * crossprod(A))
            }else if(const[2] == "orth") {        # B: orthogonality
                svd <- svd(Xb %*% krp(C, A))
                B <- svd$u %*% t(svd$v)
            }else if(const[2] == "zerocor") {     # B: zero correlation
                XF <- Xb %*% krp(C, A)
                XFmean <- apply(XF, 2, mean)
                svd <- svd(sweep(XF, 2, XFmean))
                B <- sweep(svd$u %*% t(svd$v), 2, XFmean %*% solve(AA * CC), "+")
            }
            BB = t(B) %*% B

            ## Step 3: Update mode C
            if(const[3] == "none") {              # C: no constraint
               C <- Xc %*% krp(B, A) %*% solve(crossprod(B) * crossprod(A))
            }else if(const[3] == "orth") {        # C: orthonormality
                svd <- svd(Xc %*% krp(B, A))
                C <- svd$u %*% t(svd$v)
            }else if(const[3] == "zerocor") {     # C: zero correlation
                XF <- Xc %*% krp(B, A)
                XFmean <- apply(XF, 2, mean)
                svd <- svd(sweep(XF, 2, XFmean))
                C <- sweep(svd$u %*% t(svd$v), 2, XFmean %*% solve(AA * BB), "+")
            }
            CC = t(C) %*% C

            ## Step 4: Convergence check
            if (const[3] == "none") {
                FF <- AA * BB
                f <- ssx - rrcov3way::mtrace(CC %*% FF)
            }else {
                f <- sum((Xa - tcrossprod(A, krp(C, B)))^2)
            }

            ##  Record Relative Fit
            single_iter[iter] <- abs((fold - f)/f)

            ##  TRIPLE COSINE
            if(iter%%10 == 0) {
                tripcos = min(rrcov3way::congruence(A, A) * rrcov3way::congruence(B, B) * rrcov3way::congruence(C, C))
                if (iter == 10)
                    mintripcos = tripcos
                if (tripcos < mintripcos)
                    mintripcos = tripcos
                if (trace && iter %% 1000 == 0)
                    cat(paste("Minimal Triple cosine =", tripcos, sep = " "), fill = TRUE)

                ftiter[iter/10, ] = c(f, tripcos)
            }
            if(trace && iter %% 50 == 0)
                cat(paste("f=", f, "after", iter, "iters; diff.=", fold - f, sep = " "), fill = TRUE)
        }

    ftiter <- ftiter[1:(iter/10), , drop=FALSE]           # take only the first iter/10 rows

    ## Fit percentage
    Rsq <- 1 - f/ssx
    fp <- 100 * Rsq

    ## Degeneracy problem if |tripcos|>.5
    tripcos = min(rrcov3way::congruence(A, A) * rrcov3way::congruence(B, B) * rrcov3way::congruence(C, C))
    names(tripcos) = c("Minimal triple cosine")
    if (iter < 10)
        mintripcos = tripcos

    out <- list(fit=f, fp = fp, ss=ssx, A = A, B = B, C = C, iter = iter,
        iter_opt=iter_opt, tripcos=tripcos, mintripcos=mintripcos, ftiter=ftiter,
        single_iter=single_iter[1:iter])

    out
}
