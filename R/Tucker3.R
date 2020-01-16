Tucker3 <- function(X, P=2, Q=2, R=2,
    conv=1e-6,
    center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
    scale=FALSE, scale.mode=c("B", "A", "C"),
    robust=FALSE, coda.transform=c("none", "ilr"),
    ncomp.rpca=0, alpha=0.75, maxiter=100, crit=0.975, trace=FALSE)
{
    call = match.call()
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    coda.transform <- match.arg(coda.transform)
    ilr <- coda.transform != "none"

    stopifnot(alpha <=1 & alpha >= 0.5)

    if(robust & ilr)
    {
        ret <- .Tucker3.rob.ilr(X=X, P=P, Q=Q, R=R, conv=conv,
            center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode, coda.transform=coda.transform, ncomp.rpca=ncomp.rpca, alpha=alpha, maxiter=maxiter, crit=crit, trace=trace)
    }
    else if(!robust & !ilr)
    {
        ret <- .Tucker3(X=X, P=P, Q=Q, R=R, conv=conv, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode, crit=crit, trace=trace)
    }
    else if(!robust & ilr)                  # classical for compositional data
    {
        ret <- .Tucker3.ilr(X=X, P=P, Q=Q, R=R, conv=conv, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode, coda.transform=coda.transform, crit=crit, trace=trace)
    }
    else if(robust & !ilr)                  # robust, non-compositional data
    {
        ret <- .Tucker3.rob(X=X, P=P, Q=Q, R=R, conv=conv, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode, ncomp.rpca=ncomp.rpca, alpha=alpha, maxiter=maxiter, crit=crit, trace=trace)
    }
    else
        stop("Not yet implemented!")

    ## Store in the output object the total sum of squares
    ret$ss <- sum(X^2)

    ##
    ## ret$fit will be ||X_A - A G_A kron(C',B')||^2 where X_A and G_A denote the matricized (frontal slices) data array and core array
    ## ret$fp is equal to: 100*(ss-ret$fit)/ss

    ret$call <- call
    ret
}

## ILR transformation (see package 'chemometrics'
.ilrV <- function(x)
{
    x.ilr = matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
    for (i in 1:ncol(x.ilr))
    {
        x.ilr[, i] = sqrt((i)/(i + 1)) * log(((apply(as.matrix(x[,1:i]), 1, prod))^(1/i))/(x[, i + 1]))
    }

    if (is.data.frame(x))
        x.ilr <- data.frame(x.ilr)
    return(x.ilr)
}

## Geometric mean
.gm <- function (x)
{
    if (!is.numeric(x))
        stop("x has to be a vector of class numeric")
    if(any(na.omit(x == 0))) 0 else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
}


##
## Classical (non-robust) Tucker3 (non-compositional data)
##
.Tucker3 <- function(X, P=2, Q=2, R=2, conv=1e-10, center=FALSE,
        center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale=FALSE, scale.mode=c("B", "A", "C"),
        crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    dn <- dimnames(X)
    di <- dim(X)

    ## center and scale
    X <- do3Scale(X, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)

    ret <- .TK(X=X, P=P, Q=Q, R=R, conv=conv, crit=crit)
    ret$robust <- FALSE
    ret$ilr <- FALSE

    out.sd <- .cutoff.sd(ret$A, crit=crit, robust=FALSE)
    ret$sd <- out.sd$sd
    ret$cutoff.sd <- out.sd$cutoff.sd

    class(ret) <- "tucker3"
    ret
}

##
## Classical (non-robust) Tucker3 for compositional data
##
.Tucker3.ilr <- function(X, P=2, Q=2, R=2, conv=1e-10, center=FALSE,
    center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale=FALSE, scale.mode=c("B", "A", "C"),
    coda.transform=c("ilr"),
    crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    coda.transform <- match.arg(coda.transform)

    dn <- dimnames(X)
    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]

    Xtall <- tallArray(X)
    Xilr <- .ilrV(Xtall)
    Xwideilr <- tall2wide(Xilr, I, J, K)

    J <- J - 1
    Xarrayilr <- toArray(Xwideilr, I, J, K)
    Xarrayilr <- do3Scale(Xarrayilr, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwideilr <- unfold(Xarrayilr)

    TUCKER <- .TK(Xarrayilr, P=P, Q=Q, R=R, conv=conv)

    Arew <- TUCKER$A
    Brew <- TUCKER$B
    Crew <- TUCKER$C
    Grew <- TUCKER$GA

    Arew <- Xwideilr %*% .pinv(Grew %*% t(kronecker(Crew, Brew)))
    Xfitrew <- Arew %*% Grew %*% t(kronecker(Crew, Brew))
    out.Xhat.rew <- array(Xfitrew, c(I, J, K))

    odsqrew <- apply((Xwideilr - Xfitrew)^2, 1, sum)
    RD <- sqrt(odsqrew)
    cutoffOD <- .cutoff.rd(RD, robust=FALSE)
    flag <- RD <= cutoffOD

    ## Back-transformation of loadings to clr
    V <- matrix(0, nrow=J+1, ncol=J)
    for (i in 1:ncol(V))
    {
        V[1:i, i] <- 1/i
        V[i + 1, i] <- (-1)
        V[, i] <- V[, i] * sqrt(i/(i + 1))
    }
    Bclr <- V %*% Brew

    out.sd <- .cutoff.sd(Arew, crit=crit, robust=FALSE)

    ## dimnames back
    znames <- paste("Z", 1:dim(Brew)[1], sep="")
    dimnames(Arew) <- list(dn[[1]], paste("F", 1:P, sep=""))
    dimnames(Brew) <- list(znames, paste("F", 1:Q, sep=""))
    dimnames(Bclr) <- list(dn[[2]], paste("F", 1:Q, sep=""))
    dimnames(Crew) <- list(dn[[3]], paste("F", 1:R, sep=""))
    dimnames(Grew) <- list(paste("F", 1:P, sep=""), paste("F", 1:(Q*R), sep=""))
    dimnames(out.Xhat.rew) <- list(dn[[1]], znames, dn[[3]])

    names(flag) <- dn[[1]]
    names(RD) <- dn[[1]]

    ret <- list(fit=TUCKER$fit, fp=TUCKER$fp,
            A=Arew, B=Brew, Bclr=Bclr, C=Crew, GA=Grew,
            Zhat=out.Xhat.rew, flag=flag, iter=TUCKER$iter,
            rd=RD, cutoff.rd=cutoffOD, sd=out.sd$sd, cutoff.sd=out.sd$cutoff.sd,
            robust=FALSE, ilr=TRUE)

    class(ret) <- "tucker3"

    ret
}

## Robust, no ilr transformation
.Tucker3.rob <- function(X, P=2, Q=2, R=2, conv=1e-6, center=FALSE,
        center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale=FALSE, scale.mode=c("B", "A", "C"),
        ncomp.rpca, alpha=alpha, maxiter=100, crit=0.975, trace=FALSE)
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

    ## define num. of outliers the algorithm should resists
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

    while(changeFit > conv & iter <= maxiter)
    {
        iter <- iter+1
        Xnew <- array(Xhat, c(dim(Xhat)[1], J, K))

        ## Step 2 - PARAFAC analysis
        TUCKER <- .TK(Xnew, P=P, Q=Q, R=R, conv=conv)

        Ah <- TUCKER$A        # hxP
        Bh <- TUCKER$B        # JxQ
        Ch <- TUCKER$C        # KxR
        G  <- TUCKER$GA       # PxQ*R

        Ahat <- Xwide %*% .pinv(G %*% t(kronecker(Ch, Bh)))
        Xfit <- Ahat %*% G %*% t(kronecker(Ch, Bh))
        out.Xhat <- array(Xfit, c(I, J, K))

        ## Step 4  Computation of the od
        odsq <- apply((Xwide-Xfit)^2, 1, sum)    #Xunf
        od <- sqrt(odsq)
        Hset <- sort(sort(odsq,index.return=TRUE)$ix[1:h])
        fit <- sum(odsq[Hset])
        Xhat <- Xwide[Hset,]                     #Xunf

        ## Step 5  Fit of the model
        changeFit <- if(fitprev == 0) 1 + conv else abs(fit-fitprev)/fitprev
        fitprev <- fit
    }

    ## reweighting
    cutoffOD <- .cutoff.rd(od, h)

    out.rd <- od
    flag <- (out.rd <= cutoffOD)
    Xflag <- X[flag,,]      #X
    dim <- dim(Xflag)[1]
    Xflag_unf <- matrix(Xflag,dim,J*K)
    Xflag_unf <- array(Xflag_unf, dim=c(dim(Xflag_unf)[1],J,K))

    ## run Tucker on weighted dataset
    modelTKrew <- .TK(Xflag_unf, P, Q, R)
    Arew <- modelTKrew$A
    Brew <- modelTKrew$B
    Crew <- modelTKrew$C
    Grew <- modelTKrew$GA

    Arew <- Xwide %*% .pinv(Grew %*% t(kronecker(Crew, Brew)))
    Xfitrew <- Arew %*% Grew %*% t(kronecker(Crew, Brew))
    out.Xhat.rew <- array(Xfitrew, c(I,J,K))

    odsqrew <- apply((Xwide - Xfitrew)^2, 1, sum)
    out.rd <- sqrt(odsqrew)
    cutoffOD <- .cutoff.rd(out.rd, h)
    flag <- (out.rd <= cutoffOD)
    out.sd <- .cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)


    ## dimnames back
    znames <- paste("Z", 1:dim(Brew)[1], sep="")
    dimnames(Arew) <- list(dn[[1]], paste("F", 1:P, sep=""))
    dimnames(Brew) <- list(dn[[2]], paste("F", 1:Q, sep=""))
    dimnames(Crew) <- list(dn[[3]], paste("F", 1:R, sep=""))
    dimnames(Grew) <- list(paste("F", 1:P, sep=""), paste("F", 1:(Q*R), sep=""))
    dimnames(out.Xhat.rew) <- list(dn[[1]], znames, dn[[3]])

    names(flag) <- dn[[1]]
    names(out.rd) <- dn[[1]]

    fp <- TUCKER$fp

    ret <- list(fit=fit, fp=fp, A=Arew, B=Brew, C=Crew, GA=Grew, Zhat=out.Xhat.rew,
       flag=flag, Hset=Hset, iter=iter, alpha=alpha,
       rd=out.rd, cutoff.rd=cutoffOD, sd=out.sd$sd, cutoff.sd=out.sd$cutoff.sd,
       robust=TRUE, ilr=FALSE)

    class(ret) <- "tucker3"

    ret
}

.Tucker3.rob.ilr <- function(X, P=2, Q=2, R=2, conv=1e-10, center=FALSE,
        center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"), scale=FALSE, scale.mode=c("B", "A", "C"),
        coda.transform=c("ilr"),
        ncomp.rpca, alpha=alpha, maxiter=100, crit=0.975, trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    coda.transform <- match.arg(coda.transform)

    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]

    dn <- dimnames(X)

    Xtall <- tallArray(X)
    Xilr <- .ilrV(Xtall)
    Xwideilr <- tall2wide(Xilr, I, J, K)

    J <- J-1
    Xarrayilr <- toArray(Xwideilr, I, J, K)
    Xarrayilr <- do3Scale(Xarrayilr, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwideilr <- unfold(Xarrayilr)

    ## define num. of outliers the algorithm should resists
    h <- round(alpha * I)

    ## Step 1 RobPCA XA
    if(trace)
        cat("\nStep 1. Perform robust PCA on the unfolded matrix.")

    outrobpca <- PcaHubert(Xwideilr, k=ncomp.rpca, kmax=ncol(Xwideilr), alpha=alpha, mcd=FALSE)
    Hset <- sort(sort(outrobpca@od, index.return=TRUE)$ix[1:h])
    Xhat <- Xwideilr[Hset,] #Xunf
    fitprev <- 0
    changeFit <- 1 + conv
    iter <- 0

    while(changeFit > conv & iter <= maxiter)
    {
        iter <- iter+1
        Xnew <- array(Xhat, c(dim(Xhat)[1], J, K))

        ## Step 2 - PARAFAC analysis
        TUCKER <- .TK(Xnew, P=P, Q=Q, R=R, conv=conv)

        Ah <- TUCKER$A        # hxP
        Bh <- TUCKER$B        # JxQ
        Ch <- TUCKER$C        # KxR
        G  <- TUCKER$GA       # PxQ*R

        Ahat <- Xwideilr %*% .pinv(G %*% t(kronecker(Ch, Bh)))
        Xfit <- Ahat %*% G %*% t(kronecker(Ch, Bh))
        out.Xhat <- array(Xfit, c(I, J, K))

        ## Step 4  Computation of the od
        odsq <- apply((Xwideilr-Xfit)^2, 1, sum)    #Xunf
        od <- sqrt(odsq)
        Hset <- sort(sort(odsq,index.return=TRUE)$ix[1:h])
        fit <- sum(odsq[Hset])
        Xhat <- Xwideilr[Hset,]                     #Xunf

        ## Step 5  Fit of the model
        changeFit <- if(fitprev == 0) 1 + conv else abs(fit-fitprev)/fitprev
        fitprev <- fit
    }

    ## reweighting
    cutoffOD <- .cutoff.rd(od, h)
    out.rd <- od
    flag <- (out.rd <= cutoffOD)
    Xflag <- Xarrayilr[flag,,]      #X
    dim <- dim(Xflag)[1]
    Xflag_unf <- matrix(Xflag,dim,J*K)
    Xflag_unf <- array(Xflag_unf, dim=c(dim(Xflag_unf)[1],J,K))

    ## run Tucker on weighted dataset
    modelTKrew <- .TK(Xflag_unf, P, Q, R)
    Arew <- modelTKrew$A
    Brew <- modelTKrew$B
    Crew <- modelTKrew$C
    Grew <- modelTKrew$GA

    Arew <- Xwideilr %*% .pinv(Grew %*% t(kronecker(Crew, Brew)))
    Xfitrew <- Arew %*% Grew %*% t(kronecker(Crew, Brew))
    out.Xhat.rew <- array(Xfitrew, c(I, J, K))

    odsqrew <- apply((Xwideilr - Xfitrew)^2, 1, sum) #Xunf
    out.rd <- sqrt(odsqrew)
    cutoffOD <- .cutoff.rd(out.rd, h)

    flag <- (out.rd <= cutoffOD)

    ## Back-transformation of loadings to clr
    V <- matrix(0, nrow = J+1, ncol = J)
    for (i in 1:ncol(V))
    {
        V[1:i, i] <- 1/i
        V[i + 1, i] <- (-1)
        V[, i] <- V[, i] * sqrt(i/(i + 1))
    }
    Bclr <- V %*% Brew
    out.sd <- .cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)

    ## dimnames back
    znames <- paste("Z", 1:dim(Brew)[1], sep="")
    dimnames(Arew) <- list(dn[[1]], paste("F", 1:P, sep=""))
    dimnames(Brew) <- list(znames, paste("F", 1:Q, sep=""))
    dimnames(Bclr) <- list(dn[[2]], paste("F", 1:Q, sep=""))
    dimnames(Crew) <- list(dn[[3]], paste("F", 1:R, sep=""))
    dimnames(Grew) <- list(paste("F", 1:P, sep=""), paste("F", 1:(Q*R), sep=""))
    dimnames(out.Xhat.rew) <- list(dn[[1]], znames, dn[[3]])

    names(flag) <- dn[[1]]
    names(out.rd) <- dn[[1]]

    fp <- TUCKER$fp

    ret <- list(fit=fit, fp=fp, A=Arew, B=Brew, Bclr=Bclr, C=Crew, GA=Grew, Zhat=out.Xhat.rew,
       flag=flag, Hset=Hset, iter=iter, alpha=alpha,
       rd=out.rd, cutoff.rd=cutoffOD, sd=out.sd$sd, cutoff.sd=out.sd$cutoff.sd,
       robust=TRUE, ilr=TRUE)

    class(ret) <- "tucker3"

    ret
}

## - dd     = distance-distance plot
## - comp   = paired component plot for a single mode
## - jbplot = joint biplot
## - tjplot = trajectory plot
##
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
        ret <- .allcompplot(x, ...)         # all components plot
    }

    if((which == "all" || which == "jbplot")) {
        ret <- .JBPlot(x, ...)                      # Joint biplot
    }

    if((which == "all" || which == "tjplot")) {
        ret <- .TJPlot(x, ...)                      # Trajectory biplot
    }

    invisible(ret)
}

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

    cat("\nTucker3 analysis with ",P,"x",Q,"x",R," components.\nFit value:", round(x$fp,2), "%\n")
    msg <- ""
    if(x$robust)
        msg <- paste(msg, "Robust", sep="")
    if(x$ilr){
        if(nchar(msg) > 0)
            msg <- paste(msg, ", ", sep="")
        msg <- paste(msg, "ilr-transformed", "\n", sep="")
    }
    cat(msg)

    invisible(x)
}
