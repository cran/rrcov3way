Parafac <- function(X, ncomp=2,
    center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
    scale=FALSE, scale.mode=c("B", "A", "C"),
    const="none", conv=1e-6, start="svd", maxit=10000,
    robust=FALSE, coda.transform=c("none", "ilr", "clr"),
    ncomp.rpca=0, alpha=0.75, robiter=100, crit=0.975,      # arguments for the robust parafac
    trace=FALSE)
{
    call <- match.call()

    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)

    coda.transform <- match.arg(coda.transform)
    ilr <- coda.transform != "none"
    if(coda.transform == "clr" & robust)
        stop("The robust option is not possible with 'clr' transform compositional data. Please use 'ilr'.")

    if(length(crit) != 1 || crit <= 0 || crit >= 1)
        stop("'crit' has to be a single positive number less than 1!")

    if(crit < 0.5)
        crit <- 1 - crit

    stopifnot(alpha <= 1 & alpha >= 0.5)

    if(robust & ilr)
    {
        ret <- .Parafac.rob.ilr(X=X, ncomp=ncomp, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode, const=const, conv=conv, start=start, maxit=maxit, coda.transform=coda.transform, ncomp.rpca=ncomp.rpca, alpha=alpha, robiter=robiter, crit=crit, trace=trace)
    }
    else if(!robust & !ilr)
    {
        ret <- .Parafac(X=X, ncomp=ncomp, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode, const=const, conv=conv, start=start, maxit=maxit, crit=crit, trace=trace)
    }
    else if(!robust & ilr)                  # classical for compositional data
    {
        ret <- .Parafac.ilr(X=X, ncomp=ncomp, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode, const=const, conv=conv, start=start, maxit=maxit, coda.transform=coda.transform, crit=crit, trace=trace)
    }
    else if(robust & !ilr)                  # robust, for non-compositional data
    {
        ret <- .Parafac.rob(X=X, ncomp=ncomp, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode, const=const, conv=conv, start=start, maxit=maxit, ncomp.rpca=ncomp.rpca, alpha=alpha, robiter=robiter, crit=crit, trace=trace)
    }
    else
        stop("Not yet implemented!")

    ## Total sum of squares, PARAFAC fit and fit percentage:
    ## ret$ss <- sum(X^2)
    ## ret$fit will be ||X_A - A G_A kron(C',B')||^2 where X_A and G_A denote the matricized (frontal slices) data array and core array
    ## ret$fp is equal to: 100*(ss-ret$fit)/ss

    ## Calculate the superdiagonal core g_sss and store it in GA
    ret$GA <- sqrt(colSums(ret$A^2)) * sqrt(colSums(ret$B^2)) * sqrt(colSums(ret$C^2))
    names(ret$GA) <- colnames(ret$A)

    ret$ncomp <- ncomp
    ret$call <- call
    ret
}

.cutoff.rd <- function(rd, h, crit=0.975, robust=TRUE)
{
    if(robust)
    {
        ## a) Using median and MAD
        ##  ret <-  (median(rd^(2/3)) + mad(rd^(2/3)) * qnorm(crit))^(3/2)
        ##
        ## b) Using MASS cov.mcd (will need to import - library(MASS))
        ## unimcd <- cov.mcd(rd^(2/3),quantile.used=h)
        ## ret <- sqrt(qnorm(0.975, unimcd$center, sqrt(unimcd$cov))^3)
        ##
        ## c) Using UNIMCD
         unimcd <- rrcov:::unimcd(rd^(2/3), quan=h)
         ret <- sqrt(qnorm(crit, unimcd$tmcd, unimcd$smcd)^3)

        ## d) Using UNIMCD by CovMcd
        ## unimcd <- CovMcd(rd, alpha=h/length(rd))
        ## ret <- sqrt(qnorm(crit, getCenter(unimcd), sqrt(getCov(unimcd)))^3)
        ##
    } else
        ret <- (mean(rd^(2/3)) + sd(rd^(2/3)) * qnorm(crit))^(3/2)

    ret
}

.cutoff.sd <- function(A, alpha, crit, robust=TRUE)
{
    if(robust)
    {
## VT::24.04.2020 Workaround - PcaHubert with mcd=FALSE will crash if A is one-dimensional matrix.
##                  now this is fixed in rrcov
##
##        pc <- PcaHubert(A, alpha=alpha, k=ncol(A), mcd=FALSE, crit.pca.distances=crit)
        pc <- PcaHubert(A, alpha=alpha, k=ncol(A), mcd=if(ncol(A) == 1) TRUE else FALSE, crit.pca.distances=crit)
        SD <- pc@sd
        cutoff.sd <- pc@cutoff.sd

##        A.cov <- CovMcd(A)
##        SD <- sqrt(mahalanobis(A, center=getCenter(A.cov), cov=getCov(A.cov)))
##        cutoff.sd <- sqrt(qchisq(crit, ncol(A)))
    }else
    {
        A.cov <- CovClassic(A)
        SD <- sqrt(mahalanobis(A, center=getCenter(A.cov), cov=getCov(A.cov)))
        cutoff.sd <- sqrt(qchisq(crit, ncol(A)))
    }

    list(sd=SD, cutoff.sd=cutoff.sd)
}

## Classical PARAFAC
##
##  - const: constraints (defaultes to "none"), orth=orthogonality constraints,
##      nonneg=nonnegativity, zerocor=zero correlation constraints
.Parafac <- function (X, ncomp,
    center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
    scale=FALSE, scale.mode=c("B", "A", "C"),
    const="none", conv=1e-6, start="svd", maxit=10000,
    crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)

    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]
    dn <- dimnames(X)

    X <- do3Scale(X, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwide <- unfold(X)

    ret <- cp_als(X, ncomp=ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
    A <- ret$A
    B <- ret$B
    C <- ret$C

    ## Compute RD and SD with their cutoff values
    Xfit <- A %*% t(krp(C, B))
    rdsq <- apply((Xwide - Xfit)^2, 1, sum)        # RD_i = ||X_i-Xhat_i||F;     fit = sum(RD_i^2)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, crit=crit, robust=FALSE)
    sd <- .cutoff.sd(A, crit=crit, robust=FALSE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I, J, K)

    ## dimnames back
    nfac <- paste0("F",1:ncomp)
    dimnames(A) <- list(dn[[1]],nfac)
    dimnames(B) <- list(dn[[2]],nfac)
    dimnames(C) <- list(dn[[3]],nfac)
    dimnames(Xfit) <- dn
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]

    ret <- list(fit=ret$f, fp=ret$fp, ss=ret$ss, A=A, B=B, C=C, Xhat=Xfit, const=ret$const, iter=ret$iter,
        rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
        robust=FALSE, coda.transform="none")

    class(ret) <- "parafac"
    ret
}

## Robust PARAFAC
.Parafac.rob <- function (X, ncomp, center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
    scale=FALSE, scale.mode=c("B", "A", "C"),
    const="none", conv=1e-6, start="svd", maxit=10000,
    ncomp.rpca, alpha=0.75, robiter=100, crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)

    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]
    dn <- dimnames(X)

    ## If centering and/or scaling were requested, but no (robust) function
    ##  was specified, take by default median and mad
    if(is.logical(center) && center)
        center <- median
    if(is.logical(scale) && scale)
        scale <- mad

    X <- do3Scale(X, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    ssx <- sum(X^2)
    Xwide <- unfold(X)

    Ahat <- matrix(0, I, ncomp)

    ## define number of outliers the algorithm should resists
    h <- round(alpha*I)

    ## Step 1 RobPCA of XA
    if(trace)
        cat("\nStep 1. Perform robust PCA on the unfolded matrix.")

    outrobpca <- PcaHubert(Xwide, k=ncomp.rpca, kmax=ncol(Xwide), alpha=alpha, mcd=FALSE, trace=trace)
    Hset <- sort(sort(outrobpca@od, index.return=TRUE)$ix[1:h])
    Xhat <- Xwide[Hset,]
    fitprev <- 0
    changeFit <- 1 + conv
    iter <- 0

    while (changeFit > conv & iter <= robiter)
    {
        iter <- iter + 1

        ##  Step 2 - PARAFAC analysis
        ret <- cp_als(Xhat, h, J, K, ncomp=ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
        Ah <- ret$A
        Bh <- ret$B
        Ch <- ret$C

        ## Step 3 - Fit the model
        KR <- krp(Ch, Bh)                   # Khatri-Rao product
        for(i in 1:I) {
            vJKx1 <- matrix(X[i,,], 1, J*K)
            Ahat[i,] <- pracma::pinv(KR) %*% t(vJKx1)
        }

        ## The above is equivallent to the following:
        ##  Ahat <- Xwide %*% t(pracma::pinv(KR))

        Xfit <- Ahat %*% t(KR)

        ##  Step 4  - Computation of the residual distances
        rdsq <- apply((Xwide - Xfit)^2, 1, sum)
        rd <- sqrt(rdsq)
        Hset <- sort(sort(rdsq, index.return=TRUE)$ix[1:h])
        fit <- sum(rdsq[Hset])
        Xhat <- Xwide[Hset,]

        ##  Step 5  Fit of the model
        changeFit <- if(fitprev == 0) 1 + conv else abs(fit-fitprev)/fitprev
        fitprev <- fit
    }

    ## Reweighting
    cutoff.rd <- .cutoff.rd(rd, crit=crit, h)
    flag <- rd <= cutoff.rd
    Xflag <- X[flag,,]
    Xflag_wide <- unfold(Xflag)
    ret <- cp_als(Xflag_wide, nrow(Xflag_wide), J, K, ncomp=ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
    Arew <- ret$A
    Brew <- ret$B
    Crew <- ret$C

    KRrew <- krp(Crew, Brew)            # Khatri-Rao product

    Arew <- matrix(0, I, ncomp)
    for(i in 1:I) {
        vJKx1 <- matrix(X[i,,], 1, J*K)
        Arew[i,] <- pracma::pinv(KRrew) %*% t(vJKx1)
    }

    ## The above is equivallent to the following:
    ##  Arew <- Xwide %*% t(pracma::pinv(KRrew))

    Xfit <- Arew %*%t (KRrew)
    rdsq <- apply((Xwide - Xfit)^2, 1, sum)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, crit=crit, h)
    fit <- sum(rdsq[rd <= cutoff.rd])
    fp <- 100*(1-fit/ssx)

    for(i in 1:ncomp) {
        Arew[,i] <- Arew[,i]*norm(as.matrix(Brew[,i]),type="F")*norm(as.matrix(Crew[,i]),type="F")
        Brew[,i] <- Brew[,i]/norm(as.matrix(Brew[,i]),type="F")
        Crew[,i] <- Crew[,i]/norm(as.matrix(Crew[,i]),type="F")
    }

    sd <- .cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I, J, K)

    ## dimnames back
    nfac <- paste0("F", 1:ncomp)
    dimnames(Arew) <- list(dn[[1]], nfac)
    dimnames(Brew) <- list(dn[[2]], nfac)
    dimnames(Crew) <- list(dn[[3]], nfac)
    dimnames(Xfit) <- dn
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]

    res <- list(fit=fit, fp=fp, ss=ssx, A=Arew, B=Brew, C=Crew, Xhat=Xfit, const=ret$const,
                flag=flag, Hset=Hset, iter=iter, alpha=alpha,
                rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
                pcaobj=outrobpca, robust=TRUE, coda.transform="none")

    class(res) <- "parafac"
    res
}

## Classical PARAFAC for compositional data
.Parafac.ilr <- function (X, ncomp, center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
        scale=FALSE, scale.mode=c("B", "A", "C"),
        const="none", conv=1e-6, start="svd", maxit=10000, coda.transform=c("ilr", "clr"),
        crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    coda.transform <- match.arg(coda.transform)

    ## ncomp is the number of components
    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]
    dn <- dimnames(X)

    Xilr <- if(coda.transform == "ilr") ilrArray(X)
                 else if(coda.transform == "clr") clrArray(X)
                 else NULL

    if(coda.transform == "ilr")
        J <- J - 1

    ## centering the compositions
    Xilr <- do3Scale(Xilr, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwide <- unfold(Xilr)

    ret <- cp_als(Xwide, I, J, K, ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
    A <- ret$A
    B <- ret$B
    C <- ret$C


    ## Compute RD and SD with their cutoff values
    Xfit <- A %*% t(krp(C, B))
    rdsq <- apply((Xwide - Xfit)^2, 1, sum)
    fit <- sum(rdsq)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, crit=crit, robust=FALSE)
    sd <- .cutoff.sd(A, crit=crit, robust=FALSE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I, J, K)

    ## Back-transformation of loadings to clr
    if(coda.transform == "clr")     #do nothing
        Bclr <- B
    else
    {
        V <- matrix(0, nrow = J+1, ncol = J)
        for(i in seq_len(ncol(V)))
        {
            V[1:i, i] <- 1/i
            V[i + 1, i] <- (-1)
            V[, i] <- V[, i] * sqrt(i/(i + 1))
        }
        Bclr <- V %*% B
    }

    ## dimnames back
    znames <- paste0("Z", 1:dim(B)[1])
    nfac <- paste0("F", 1:ncomp)
    dimnames(A) <- list(dn[[1]], nfac)
    dimnames(B) <- if(coda.transform == "clr") list(dn[[2]], nfac) else list(znames, nfac)
    dimnames(Bclr) <- list(dn[[2]], nfac)
    dimnames(C) <- list(dn[[3]], nfac)
    dimnames(Xfit) <- list(dn[[1]], znames, dn[[3]])
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]

    res <- list(fit=ret$f, fp=ret$fp, ss=ret$ss, A=A, B=B, Bclr=Bclr, C=C, Xhat=Xfit, const=ret$const, iter=ret$iter,
        rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
        robust=FALSE, coda.transform=coda.transform)

    class(res) <- "parafac"
    res
}

## Robust PARAFAC for compositional data
.Parafac.rob.ilr <- function (X, ncomp, center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
        scale=FALSE, scale.mode=c("B", "A", "C"),
        const="none", conv=1e-6, start="svd", maxit=10000, coda.transform=c("ilr"),
        ncomp.rpca, alpha=0.75, robiter=100, crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    coda.transform <- match.arg(coda.transform)

    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]
    dn <- dimnames(X)

    Xilr <- ilrArray(X)
    J <- J - 1

    ## If centering and/or scaling were requested, but no (robust) function
    ##  was specified, take by default median and mad
    if(is.logical(center) && center)
        center <- median
    if(is.logical(scale) && scale)
        scale <- mad

    ## centering the compositions
    Xilr <- do3Scale(Xilr, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    ssx <- sum(Xilr^2)
    Xwide <- unfold(Xilr)

    Ahat <- matrix(0,dim(X)[1],ncomp)

    ## define number of outliers the algorithm should resists
    h <- round(alpha * dim(X)[1])

    ## Step 1 RobPCA XA
    outrobpca <- PcaHubert(Xwide, k=ncomp.rpca, kmax=ncol(Xwide), alpha=alpha, mcd=FALSE)
    Hset <- sort(sort(outrobpca@od, index.return=TRUE)$ix[1:h])
    Xhat <-Xwide[Hset,]
    fitprev <- 0
    changeFit <- 1 + conv
    iter <- 0
    while (changeFit > conv & iter <= robiter) {
        iter <- iter+1

        ## Step 2 - PARAFAC analysis
        ret <- cp_als(Xhat, h, J, K, ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
        Ah <- ret$A
        Bh <- ret$B
        Ch <- ret$C

        KR <- krp(Ch, Bh)  # Khatri-Rao product
        for(i in 1:dim(X)[1]) {
            vJKx1 <- matrix(Xilr[i, ,], 1, J*K)
            Ahat[i,] <- pracma::pinv(KR) %*% t(vJKx1)
        }
        Xfit <- Ahat %*%t (KR)

        ## Step 4  Computation of the rd
        rdsq <- apply((Xwide - Xfit)^2, 1, sum)
        rd <- sqrt(rdsq)
        Hset <- sort(sort(rdsq, index.return=TRUE)$ix[1:h])
        fit <- sum(rdsq[Hset])
        Xhat <- Xwide[Hset,]

        ## Step 5  Fit of the model
        changeFit <- if(fitprev == 0) 1 + conv else abs(fit-fitprev)/fitprev
        fitprev <- fit
    }

    ## Reweighting
    cutoff.rd <- .cutoff.rd(rd, crit=crit, h)
    flag <- rd <= cutoff.rd
    Xflag <- Xilr[flag,,]
    Xflag_wide <- unfold(Xflag)
    ret <- cp_als(Xflag_wide, nrow(Xflag_wide), J, K, ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
    Arew <- ret$A
    Brew <- ret$B
    Crew <- ret$C

    KR <- krp(Crew, Brew)  # Khatri-Rao product
    Arew <- matrix(0, I, ncomp)
    for(i in 1:I) {
        vJKx1 <- matrix(Xilr[i, ,], 1, J*K)
        Arew[i,] <- pracma::pinv(KR) %*%t (vJKx1)
    }
    Xfit <- Arew %*% t(KR)
    rdsq <- apply((Xwide - Xfit)^2, 1, sum)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, crit=crit, h)
    fit <- sum(rdsq[rd <= cutoff.rd])
    fp <- 100*(1-fit/ssx)

    for(i in 1:ncomp) {
        Arew[,i] <- Arew[,i] * norm(as.matrix(Brew[,i]), type="F") * norm(as.matrix(Crew[,i]), type="F")
        Brew[,i] <- Brew[,i] / norm(as.matrix(Brew[,i]), type="F")
        Crew[,i] <- Crew[,i] / norm(as.matrix(Crew[,i]), type="F")
    }

    ## Back-transformation of loadings to clr
    V <- matrix(0, nrow = J+1, ncol = J)
    for (i in seq_len(ncol(V))) {
      V[1:i, i] <- 1/i
      V[i + 1, i] <- (-1)
      V[, i] <- V[, i] * sqrt(i/(i + 1))
    }
    Bclr <- V %*% Brew

    sd <- .cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I, J, K)

    ## dimnames back
    nfac <- paste0("F", 1:ncomp)
    znames <- paste0("Z", 1:dim(Brew)[1])
    dimnames(Arew) <- list(dn[[1]], nfac)
    dimnames(Brew) <- list(znames, nfac)
    dimnames(Bclr) <- list(dn[[2]], nfac)
    dimnames(Crew) <- list(dn[[3]], nfac)
    dimnames(Xfit) <- list(dn[[1]], znames, dn[[3]])
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]

    res <- list(fit=fit, fp=fp, ss=ssx, A=Arew, B=Brew, Bclr=Bclr, C=Crew, Xhat=Xfit, const=ret$const,
            flag=flag, Hset=Hset, iter=iter, alpha=alpha,
            rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
            pcaobj=outrobpca, robust=TRUE, coda.transform=coda.transform)

    class(res) <- "parafac"
    res
}

## - dd     = distance-distance plot
## - comp   = paired component plot for a single mode
## - percomp= per component plot
## - allcomp= all components plot
##
plot.parafac <- function(x, which=c("dd", "comp", "percomp", "allcomp", "all"), ask = (which=="all" && dev.interactive(TRUE)), id.n, ...)
{
    which <- match.arg(which)
    op <- if(ask) par(ask = TRUE) else list()
    on.exit(par(op))

    if((which == "all" || which == "dd")) {
        ret <- .ddplot(x, id.n=id.n, ...)           # distance-distance plot
    }

    if((which == "all" || which == "comp")) {
        ret <- .compplot.parafac(x, ...)            # paired components plot
    }

    if((which == "all" || which == "percomp")) {
        ret <- .percompplot.parafac(x, ...)         # per-component plot
    }

    if((which == "all" || which == "allcomp")) {
        ret <- .allcompplot(x, ...)                 # all-component plot
    }

    invisible(ret)
}

print.parafac <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    ncomp <- dim(x$A)[2]

    cat("\nPARAFAC analysis with ", ncomp, " components.\nFit value:", x$fit,
        "\nFit percentage:", round(x$fp,2), "%\n")
    msg <- ""
    if(x$robust)
        msg <- paste0(msg, "Robust")
    if(x$coda.transform != "none"){
        if(nchar(msg) > 0)
            msg <- paste0(msg, ", ")
        tr <- if(x$coda.transform == "clr") "clr-transformed" else "ilr-transformed"
        msg <- paste0(msg, tr, "\n")
    }
    cat(msg, "\n")

    invisible(x)
}
