Tucker3 <- function(X, P=2, Q=2, R=2,
    center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
    scale=FALSE, scale.mode=c("B", "A", "C"),
    conv=1e-6, start="svd",
    robust=FALSE, coda.transform=c("none", "ilr", "clr"),
    ncomp.rpca=0, alpha=0.75, robiter=100, crit=0.975, trace=FALSE)
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

    stopifnot(alpha <=1 & alpha >= 0.5)

    ret <- if(robust & ilr) .Tucker3.rob.ilr(X=X, P=P, Q=Q, R=R, conv=conv,
            start=start, center=center, center.mode=center.mode,
            scale=scale, scale.mode=scale.mode, coda.transform=coda.transform,
            ncomp.rpca=ncomp.rpca, alpha=alpha, robiter=robiter, crit=crit,
            trace=trace)
           else if(!robust & !ilr) .Tucker3(X=X, P=P, Q=Q, R=R, conv=conv,
            start=start, center=center, center.mode=center.mode,
            scale=scale, scale.mode=scale.mode, crit=crit, trace=trace)
           else if(!robust & ilr) .Tucker3.ilr(X=X, P=P, Q=Q, R=R, conv=conv,
            start=start, center=center, center.mode=center.mode,
            scale=scale, scale.mode=scale.mode, coda.transform=coda.transform,
            crit=crit, trace=trace)
           else if(robust & !ilr) .Tucker3.rob(X=X, P=P, Q=Q, R=R, conv=conv,
            start=start, center=center, center.mode=center.mode,
            scale=scale, scale.mode=scale.mode, ncomp.rpca=ncomp.rpca,
            alpha=alpha, robiter=robiter, crit=crit, trace=trace)

    ## Total sum of squares, TUCKER3 fit and fit percentage:
    ## ret$ss <- sum(X^2)
    ## ret$fit will be ||X_A - A G_A kron(C',B')||^2 where X_A and G_A denote the matricized (frontal slices) data array and core array
    ## ret$fp is equal to: 100*(ss-ret$fit)/ss

    ret$call <- call
    ret
}

##
## Classical (non-robust) Tucker3 (non-compositional data)
##
.Tucker3 <- function(X, P=2, Q=2, R=2, conv, start, center=FALSE,
        center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale=FALSE, scale.mode=c("B", "A", "C"),
        crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]
    dn <- dimnames(X)

    ## center and scale
    X <- do3Scale(X, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwide <- unfold(X)

    ret <- t3_als(Xwide, I, J, K, P=P, Q=Q, R=R, start=start, conv=conv)
    A <- ret$A
    B <- ret$B
    C <- ret$C
    GA <- ret$GA

    ## Compute RD and SD with their cutoff values; flag observations as outliers
    Xfit <- A %*% GA %*% t(kronecker(C, B))
    rdsq <- apply((Xwide - Xfit)^2, 1, sum)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, crit=crit, robust=FALSE)   # (mean(rd^(2/3)) + sd(rd^(2/3)) * qnorm(crit))^(3/2)
    sd <- .cutoff.sd(A, crit=crit, robust=FALSE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I, J, K)

    ## compute "intrinsic eigenvalues" eigenvalues for A-mode:
    La <- GA %*% t(GA)
    Y <- permute(GA, P, Q, R)
    Lb <- Y %*% t(Y)
    Y <- permute(Y, Q, R, P)
    Lc <- Y %*% t(Y)

    ## dimnames back
    dimnames(A) <- list(dn[[1]], paste("F", 1:P, sep=""))
    dimnames(B) <- list(dn[[2]], paste("F", 1:Q, sep=""))
    dimnames(C) <- list(dn[[3]], paste("F", 1:R, sep=""))
    dimnames(GA) <- list(paste("F", 1:P, sep=""), paste("F", 1:(Q*R), sep=""))
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]
    dimnames(La) <- list(paste("F", 1:P, sep=""), paste("F", 1:P, sep=""))
    dimnames(Lb) <- list(paste("F", 1:Q, sep=""), paste("F", 1:Q, sep=""))
    dimnames(Lc) <- list(paste("F", 1:R, sep=""), paste("F", 1:R, sep=""))

    ret <- list(fit=ret$fit, fp=ret$fp, ss=ret$ss, A=A, B=B, C=C, GA=GA,
        La=La, Lb=Lb, Lc=Lc, Xhat=Xfit, flag=flag, const=ret$const, iter=ret$iter,
        rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
        robust=FALSE, coda.transform="none")

    class(ret) <- "tucker3"
    ret
}

## Robust, no ilr transformation
.Tucker3.rob <- function(X, P=2, Q=2, R=2, conv, start, center=FALSE,
        center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale=FALSE, scale.mode=c("B", "A", "C"),
        ncomp.rpca, alpha=alpha, robiter, crit=0.975, trace=FALSE)
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

    ## define num. of outliers the algorithm should resists
    h <- round(alpha * I)

    ## Step 1 RobPCA of XA
    if(trace)
        cat("\nStep 1. Perform robust PCA on the unfolded matrix.")

    outrobpca <- PcaHubert(Xwide, k=ncomp.rpca, kmax=ncol(Xwide), alpha=alpha, mcd=FALSE, trace=trace)
    Hset <- sort(sort(outrobpca@od, index.return=TRUE)$ix[1:h])
    Xhat <- Xwide[Hset,]
    fitprev <- 0
    changeFit <- 1 + conv
    iter <- 0

    while(changeFit > conv & iter <= robiter)
    {
        iter <- iter+1

        ## Step 2 - TUCKER3 analysis
        ret <- t3_als(Xhat, h, J, K, P=P, Q=Q, R=R, start=start, conv=conv)

        Ah <- ret$A        # hxP
        Bh <- ret$B        # JxQ
        Ch <- ret$C        # KxR
        G  <- ret$GA       # PxQ*R

        Ahat <- Xwide %*% pracma::pinv(G %*% t(kronecker(Ch, Bh)))
        Xfit <- Ahat %*% G %*% t(kronecker(Ch, Bh))

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
    cutoff.rd <- .cutoff.rd(rd, h, crit=crit)
    flag <- (rd <= cutoff.rd)
    Xflag <- X[flag,,]
    Xflag_wide <- unfold(Xflag)

    ## run Tucker on weighted dataset
    ret <- t3_als(Xflag_wide, nrow(Xflag_wide), J, K, P, Q, R, start=start, conv=1e-10)
    Arew <- ret$A
    Brew <- ret$B
    Crew <- ret$C
    Grew <- ret$GA

    Arew <- Xwide %*% pracma::pinv(Grew %*% t(kronecker(Crew, Brew)))
    Xfit <- Arew %*% Grew %*% t(kronecker(Crew, Brew))

    rdsq <- apply((Xwide - Xfit)^2, 1, sum)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, h, crit=crit)
    fit <- sum(rdsq[rd <= cutoff.rd])
    fp <- 100*(1-fit/ssx)

    sd <- .cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I,J,K)

    ## compute "intrinsic eigenvalues" eigenvalues for A-mode:
    La <- Grew %*% t(Grew)
    Y <- permute(Grew, P, Q, R)
    Lb <- Y %*% t(Y)
    Y <- permute(Y, Q, R, P)
    Lc <- Y %*% t(Y)

    ## dimnames back
    dimnames(Arew) <- list(dn[[1]], paste0("F", 1:P))
    dimnames(Brew) <- list(dn[[2]], paste0("F", 1:Q))
    dimnames(Crew) <- list(dn[[3]], paste0("F", 1:R))
    dimnames(Grew) <- list(paste0("F", 1:P), paste0("F", 1:(Q*R)))
    dimnames(Xfit) <- dn
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]
    dimnames(La) <- list(paste("F", 1:P, sep=""), paste("F", 1:P, sep=""))
    dimnames(Lb) <- list(paste("F", 1:Q, sep=""), paste("F", 1:Q, sep=""))
    dimnames(Lc) <- list(paste("F", 1:R, sep=""), paste("F", 1:R, sep=""))

    ret <- list(fit=fit, fp=fp, ss=ssx, A=Arew, B=Brew, C=Crew, GA=Grew,
       La=La, Lb=Lb, Lc=Lc, Xhat=Xfit,
       flag=flag, Hset=Hset, iter=iter, alpha=alpha,
       rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
       pcaobj=outrobpca, robust=TRUE, coda.transform="none")

    class(ret) <- "tucker3"

    ret
}

##
## Classical (non-robust) Tucker3 for compositional data
##
.Tucker3.ilr <- function(X, P=2, Q=2, R=2, conv, start, center=FALSE,
    center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale=FALSE, scale.mode=c("B", "A", "C"),
    coda.transform=c("ilr", "clr"),
    crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    coda.transform <- match.arg(coda.transform)

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

    ## centering and scaling the compositions
    Xilr <- do3Scale(Xilr, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwide <- unfold(Xilr)

    ret <- t3_als(Xwide, I, J, K, P=P, Q=Q, R=R, start=start, conv=conv)

    A <- ret$A
    B <- ret$B
    C <- ret$C
    GA <- ret$GA

    ## Compute RD and SD with their cutoff values; flag observations as outliers
    Xfit <- A %*% GA %*% t(kronecker(C, B))
    rdsq <- apply((Xwide - Xfit)^2, 1, sum)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, crit=crit, robust=FALSE)

    sd <- .cutoff.sd(A, crit=crit, robust=FALSE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I, J, K)

    ## compute "intrinsic eigenvalues" eigenvalues for A-mode:
    La <- GA %*% t(GA)
    Y <- permute(GA, P, Q, R)
    Lb <- Y %*% t(Y)
    Y <- permute(Y, Q, R, P)
    Lc <- Y %*% t(Y)

    ## Back-transformation of loadings to clr
    if(coda.transform == "clr")     # do nothing
        Bclr <- B
    else
    {
        V <- matrix(0, nrow=J+1, ncol=J)
        for (i in seq_len(ncol(V)))
        {
            V[1:i, i] <- 1/i
            V[i + 1, i] <- (-1)
            V[, i] <- V[, i] * sqrt(i/(i + 1))
        }
        Bclr <- V %*% B
    }

    ## dimnames back
    znames <- paste0("Z", 1:dim(B)[1])
    dimnames(A) <- list(dn[[1]], paste0("F", 1:P))
    dimnames(B) <- if(coda.transform == "clr") list(dn[[2]], paste0("F", 1:Q)) else list(znames, paste0("F", 1:Q))
    dimnames(Bclr) <- list(dn[[2]], paste0("F", 1:Q))
    dimnames(C) <- list(dn[[3]], paste0("F", 1:R))
    dimnames(GA) <- list(paste0("F", 1:P), paste0("F", 1:(Q*R)))
    dimnames(Xfit) <- list(dn[[1]], znames, dn[[3]])
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]
    dimnames(La) <- list(paste("F", 1:P, sep=""), paste("F", 1:P, sep=""))
    dimnames(Lb) <- list(paste("F", 1:Q, sep=""), paste("F", 1:Q, sep=""))
    dimnames(Lc) <- list(paste("F", 1:R, sep=""), paste("F", 1:R, sep=""))

    ret <- list(fit=ret$fit, fp=ret$fp, ss=ret$ss,
            A=A, B=B, Bclr=Bclr, C=C, GA=GA,
            La=La, Lb=Lb, Lc=Lc, Xhat=Xfit, flag=flag, iter=ret$iter,
            rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
            robust=FALSE, coda.transform=coda.transform)

    class(ret) <- "tucker3"

    ret
}

.Tucker3.rob.ilr <- function(X, P=2, Q=2, R=2, conv, start, center=FALSE,
        center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale=FALSE, scale.mode=c("B", "A", "C"),
        coda.transform=c("ilr"),
        ncomp.rpca, alpha=alpha, robiter, crit=0.975, trace=FALSE)
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

    ## define num of outliers the algorithm should resists
    h <- round(alpha * I)

    ## Step 1 RobPCA XA
    if(trace)
        cat("\nStep 1. Perform robust PCA on the unfolded matrix.")

    outrobpca <- PcaHubert(Xwide, k=ncomp.rpca, kmax=ncol(Xwide), alpha=alpha, mcd=FALSE)
    Hset <- sort(sort(outrobpca@od, index.return=TRUE)$ix[1:h])
    Xhat <- Xwide[Hset,] #Xunf
    fitprev <- 0
    changeFit <- 1 + conv
    iter <- 0

    while(changeFit > conv & iter <= robiter)
    {
        iter <- iter + 1

        ## Step 2 - PARAFAC analysis
        ret <- t3_als(Xhat, h, J, K, P=P, Q=Q, R=R, start=start, conv=conv)

        Ah <- ret$A        # hxP
        Bh <- ret$B        # JxQ
        Ch <- ret$C        # KxR
        G  <- ret$GA       # PxQ*R

        Ahat <- Xwide %*% pracma::pinv(G %*% t(kronecker(Ch, Bh)))
        Xfit <- Ahat %*% G %*% t(kronecker(Ch, Bh))

        ## Step 4  Computation of the RD
        rdsq <- apply((Xwide - Xfit)^2, 1, sum)
        rd <- sqrt(rdsq)
        Hset <- sort(sort(rdsq, index.return=TRUE)$ix[1:h])
        fit <- sum(rdsq[Hset])
        Xhat <- Xwide[Hset,]

        ## Step 5  Fit of the model
        changeFit <- if(fitprev == 0) 1 + conv else abs(fit-fitprev)/fitprev
        fitprev <- fit
    }

    ## reweighting
    cutoff.rd <- .cutoff.rd(rd, h, crit=crit)
    flag <- rd <= cutoff.rd
    Xflag <- Xilr[flag,,]
    Xflag_wide <- unfold(Xflag)

    ## run Tucker on weighted dataset
    ret <- t3_als(Xflag_wide, nrow(Xflag_wide), J, K, P, Q, R, start=start, conv=1e-10)
    Arew <- ret$A
    Brew <- ret$B
    Crew <- ret$C
    Grew <- ret$GA

    Arew <- Xwide %*% pracma::pinv(Grew %*% t(kronecker(Crew, Brew)))
    Xfit <- Arew %*% Grew %*% t(kronecker(Crew, Brew))

    rdsq<- apply((Xwide - Xfit)^2, 1, sum)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, h, crit=crit)
    fit <- sum(rdsq[rd <= cutoff.rd])
    fp <- 100*(1-fit/ssx)

    sd <- .cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I, J, K)

    ## compute "intrinsic eigenvalues" eigenvalues for A-mode:
    La <- Grew %*% t(Grew)
    Y <- permute(Grew, P, Q, R)
    Lb <- Y %*% t(Y)
    Y <- permute(Y, Q, R, P)
    Lc <- Y %*% t(Y)

    ## Back-transformation of loadings to clr
    V <- matrix(0, nrow = J+1, ncol = J)
    for (i in seq_len(ncol(V)))
    {
        V[1:i, i] <- 1/i
        V[i + 1, i] <- (-1)
        V[, i] <- V[, i] * sqrt(i/(i + 1))
    }
    Bclr <- V %*% Brew


    ## dimnames back
    znames <- paste0("Z", 1:dim(Brew)[1])
    dimnames(Arew) <- list(dn[[1]], paste0("F", 1:P))
    dimnames(Brew) <- list(znames, paste0("F", 1:Q))
    dimnames(Bclr) <- list(dn[[2]], paste0("F", 1:Q))
    dimnames(Crew) <- list(dn[[3]], paste0("F", 1:R))
    dimnames(Grew) <- list(paste0("F", 1:P), paste0("F", 1:(Q*R)))
    dimnames(Xfit) <- list(dn[[1]], znames, dn[[3]])
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]
    dimnames(La) <- list(paste("F", 1:P, sep=""), paste("F", 1:P, sep=""))
    dimnames(Lb) <- list(paste("F", 1:Q, sep=""), paste("F", 1:Q, sep=""))
    dimnames(Lc) <- list(paste("F", 1:R, sep=""), paste("F", 1:R, sep=""))

    ret <- list(fit=fit, fp=fp, ss=ssx, A=Arew, B=Brew, Bclr=Bclr, C=Crew, GA=Grew,
       La=La, Lb=Lb, Lc=Lc, Xhat=Xfit,
       flag=flag, Hset=Hset, iter=iter, alpha=alpha,
       rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
       pcaobj=outrobpca, robust=TRUE, coda.transform=coda.transform)

    class(ret) <- "tucker3"

    ret
}

## - dd     = distance-distance plot
## - comp   = paired component plot for a single mode
## - jbplot = joint biplot
## - tjplot = trajectory plot
##
#' @export
plot.tucker3 <- function(x, which=c("dd", "comp", "allcomp", "jbplot", "tjplot", "all"), ask = (which=="all" && dev.interactive(TRUE)), id.n, ...)
{
    which <- match.arg(which)
    op <- if(ask) par(ask = TRUE) else list()
    on.exit(par(op))

    if((which == "all" || which == "dd")) {
        ret <- .ddplot(x, id.n=id.n, ...)           # distance-distance plot
    }

    if((which == "all" || which == "comp")) {
        ret <- .compplot.tucker3(x, ...)            # paired components plot
    }

    if((which == "all" || which == "allcomp")) {
        ret <- .allcompplot(x, ...)                 # all components plot
    }

    if((which == "all" || which == "jbplot")) {
        ret <- .JBPlot(x, ...)                      # Joint biplot
    }

    if((which == "all" || which == "tjplot")) {
        ret <- .TJPlot(x, ...)                      # Trajectory biplot
    }

    invisible(ret)
}

#' @export
print.tucker3 <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    P <- dim(x$A)[2]
    Q <- dim(x$B)[2]
    R <- dim(x$C)[2]

    cat("\nTucker3 analysis with ",P,"x",Q,"x",R," components.\nFit value:", x$fit, "\nFit percentage:", round(x$fp,2), "%\n")
    msg <- ""
    if(x$robust)
        msg <- paste0(msg, "Robust")
    if(x$coda.transform != "none"){
        if(nchar(msg) > 0)
            msg <- paste0(msg, ", ")

        tr <- if(x$coda.transform == "clr") "clr-transformed" else "ilr-transformed"
        msg <- paste0(msg, tr, "\n")
    }
    cat(msg)

    invisible(x)
}
