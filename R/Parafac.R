##  Return:
##  .Parafac:       list(fit=fit, A=A, B=B, C=C, Xhat=Xfit, RD=RD)
##  .Parafac.rob:   list(fit=fit, A=Arew, B=Brew, C=Crew, Xhat=Xhat.rew,
##                      flag=flag, Hset=Hset, iter=iter, RD=out.rd)
##  .Parafac.ilr:   list(fit=fit, A=A, B=B, Bclr=Bclr, C=C, Zhat=Zfit, RD=RD)
##  .Parafac.rob.ilr: list(fit=fit, A=Arew, B=Brew, Bclr=Bclr, C=Crew,
##                      Zhat=Xhat.rew, flag=flag, Hset=Hset,
##                      iter=iter, RD=out.rd)
###############################
##
Parafac <- function(X, ncomp=2,
    center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
    scale=FALSE, scale.mode=c("B", "A", "C"),
    const="none", conv=1e-6, start="svd", maxit=10000,
    robust=FALSE, coda.transform=c("none", "ilr"),
    ncomp.rpca=0, alpha=0.75, robiter=100, crit=0.975,      # arguments for the robust parafac
    trace=FALSE)
{
    center.mode <- match.arg(center.mode)
    scale.mode <- match.arg(scale.mode)
    coda.transform <- match.arg(coda.transform)
    ilr <- coda.transform != "none"

    call <- match.call()
    stopifnot(alpha <=1 & alpha >= 0.5)

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

    ## Store in the output object the total sum of squares
    ret$ss <- sum(X^2)

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
        pc <- PcaHubert(A, alpha=alpha, k=ncol(A), mcd=FALSE)
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
##  - const: constraints 9defaultes to "none"), orth=orthogonality constraints,
##      nonneg=nonnegativity, zerocor=zero correlation constraints
.Parafac <- function (X, ncomp, center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
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

    modelPar <- cp_als(X, ncomp=ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
    A <- modelPar$A
    B <- modelPar$B
    C <- modelPar$C

    KR <- krp(C, B)  # Khatri-Rao product
    Xfitw <- A %*% t(KR)
    Xfit <- array(Xfitw, c(I,J,K))

    rdsq <- apply((Xwide-Xfitw)^2,1,sum)
    fit <- sum(rdsq)
    rd <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(rd, robust=FALSE)

    out.sd <- .cutoff.sd(A, crit=crit, robust=FALSE)

    ## dimnames back
    nfac <- paste("F",1:ncomp,sep="")
    dimnames(A) <- list(dn[[1]],nfac)
    dimnames(B) <- list(dn[[2]],nfac)
    dimnames(C) <- list(dn[[3]],nfac)
    dimnames(Xfit) <- list(dn[[1]],dn[[2]],dn[[3]])
    names(rd) <- dn[[1]]

    ret <- list(fit=modelPar$f, fp=modelPar$fp, A=A, B=B, C=C, Xhat=Xfit, const=modelPar$const,
        rd=rd, cutoff.rd=cutoff.rd, sd=out.sd$sd, cutoff.sd=out.sd$cutoff.sd, robust=FALSE, ilr=FALSE)
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

    ## ncomp is the number of components
    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]

    dn <- dimnames(X)

    X <- do3Scale(X, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwide <- unfold(X)

    Ahat <- matrix(0, dim(X)[1], ncomp)

    ## define num. of outliers the algorithm should resists
    h <- round(alpha*dim(X)[1])

    ## Step 1 RobPCA XA
    outrobpca <- PcaHubert(Xwide, k=ncomp.rpca, kmax=ncol(Xwide), alpha=alpha, mcd=FALSE)
    Hset <- sort(sort(outrobpca@od, index.return=TRUE)$ix[1:h])
    Xhat <-Xwide[Hset,]
    fitprev <- 0
    changeFit <- 1 + conv
    iter <- 0
    while (changeFit > conv & iter <= robiter)
    {
        iter <- iter+1

        ##  Step 2 - PARAFAC analysis
        modelPar <- cp_als(Xhat, h, J, K, ncomp=ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
        Ah<-modelPar$A
        Bh<-modelPar$B
        Ch<-modelPar$C

        ## Step 3 - Fit the model
        KR <- krp(Ch, Bh)  # Khatri-Rao product
        for(i in 1:dim(X)[1])
        {
            vJKx1 <- matrix(X[i,,], 1, J*K)
            Ahat[i,] <- .pinv(KR) %*% t(vJKx1)
        }
        Xfit <- Ahat %*% t(KR)

        ##  Step 4  - Computation of the residual distances
        rdsq <- apply((Xwide-Xfit)^2, 1, sum)
        rd <- sqrt(rdsq)
        Hset <- sort(sort(rdsq, index.return=TRUE)$ix[1:h])
        fit <- sum(rdsq[Hset])
        Xhat <- Xwide[Hset,]#Xunf

        ##  Step 5  Fit of the model
        changeFit <- if(fitprev == 0) 1 + conv else abs(fit-fitprev)/fitprev
        fitprev <- fit
    }

    ## reweighting
    cutoffOD <- .cutoff.rd(rd, h)
    flag <- rd <= cutoffOD
    Xflag <- X[flag,,]
    dim <- dim(Xflag)[1]
    Xflag_unf <- matrix(Xflag, dim, J*K)
    modelParrew <- cp_als(Xflag_unf, dim, J, K, ncomp=ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
    Arew <- modelParrew$A
    Brew <- modelParrew$B
    Crew <- modelParrew$C

    KRrew <- krp(Crew, Brew)  # Khatri-Rao product
    Arew <- matrix(0,dim(X)[1],ncomp)
    for(i in 1:dim(X)[1])
    {
        vJKx1 <- matrix(X[i,,],1,J*K)#X
        Arew[i,] <- .pinv(KRrew)%*%t(vJKx1)
    }
    Xfitrew <- Arew%*%t(KRrew)
    odsqrew <- apply((Xwide-Xfitrew)^2, 1, sum) #Xunf
    out.rd <- odrew <- sqrt(odsqrew)
    cutoffOD <- .cutoff.rd(odrew, h)
    flag <- odrew <= cutoffOD
    Xhat.rew <- array(Xfitrew, c(I,J,K))

    for(i in 1:ncomp)
    {
        Arew[,i]<-Arew[,i]*norm(as.matrix(Brew[,i]),type="F")*norm(as.matrix(Crew[,i]),type="F")
        Brew[,i]<-Brew[,i]/norm(as.matrix(Brew[,i]),type="F")
        Crew[,i]<-Crew[,i]/norm(as.matrix(Crew[,i]),type="F")
    }

    out.sd <- .cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)

    ## dimnames back
    nfac <- paste("F",1:ncomp,sep="")
    dimnames(Arew) <- list(dn[[1]],nfac)
    dimnames(Brew) <- list(dn[[2]],nfac)
    dimnames(Crew) <- list(dn[[3]],nfac)
    dimnames(Xhat.rew) <- list(dn[[1]],dn[[2]],dn[[3]])
    names(flag) <- dn[[1]]
    names(out.rd) <- dn[[1]]

    res <- list(fit=fit, fp=modelPar$fp, A=Arew, B=Brew, C=Crew, Xhat=Xhat.rew, const=modelParrew$const,
                flag=flag, Hset=Hset, iter=iter, alpha=alpha,
                rd=out.rd, cutoff.rd=cutoffOD, sd=out.sd$sd, cutoff.sd=out.sd$cutoff.sd,
                pcaobj=outrobpca,
                robust=TRUE, ilr=FALSE)

    class(res) <- "parafac"
    res
}

## Classical PARAFAC for compositional data
.Parafac.ilr <- function (X, ncomp, center=FALSE, center.mode=c("A", "B", "C", "AB", "AC", "BC", "ABC"),
        scale=FALSE, scale.mode=c("B", "A", "C"),
        const="none", conv=1e-6, start="svd", maxit=10000, coda.transform=c("ilr"),
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

    Xtall <- tallArray(X)
    Xilr <- .ilrV(Xtall)
    Xwideilr <- tall2wide(Xilr, I, J, K)

    #############################################################
    J <- J-1
    Xarrayilr <- array(Xwideilr, c(I, J, K))

    ## centering the compositions
    Xarrayilr <- do3Scale(Xarrayilr, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwideilr <- unfold(Xarrayilr)

    modelPar <- cp_als(Xwideilr, I, J, K, ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
    A <- modelPar$A
    B <- modelPar$B
    C <- modelPar$C

    KR <- krp(C, B)  # Khatri-Rao product
    Zfitw <- A%*%t(KR)
    Zfit <- array(Zfitw,c(I,J,K))

    rdsq <- apply((Xwideilr-Zfitw)^2,1,sum)
    fit <- sum(rdsq)
    RD <- sqrt(rdsq)
    cutoff.rd <- .cutoff.rd(RD, robust=FALSE)

    ## Back-transformation of loadings to clr
    V <- matrix(0, nrow = J+1, ncol = J)
    for(i in 1:ncol(V))
    {
        V[1:i, i] <- 1/i
        V[i + 1, i] <- (-1)
        V[, i] <- V[, i] * sqrt(i/(i + 1))
    }
    Bclr <- V %*% B

    out.sd <- .cutoff.sd(A, crit=crit, robust=FALSE)

    ## dimnames back
    nfac <- paste("F",1:ncomp,sep="")
    dimnames(A) <- list(dn[[1]],nfac)
    dimnames(B) <- list(NULL,nfac)
    dimnames(Bclr) <- list(dn[[2]],nfac)
    dimnames(C) <- list(dn[[3]],nfac)
    dimnames(Zfit) <- list(dn[[1]],NULL,dn[[3]])
    names(RD) <- dn[[1]]

    res <- list(fit=fit, fp=modelPar$fp, A=A, B=B, Bclr=Bclr, C=C, Zhat=Zfit, const=modelPar$const,
        rd=RD, cutoff.rd=cutoff.rd, sd=out.sd$sd, cutoff.sd=out.sd$cutoff.sd, robust=FALSE, ilr=TRUE)
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

    ## ncomp is the number of components
    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]

    dn <- dimnames(X)

    Xtall <- tallArray(X)
    Xilr <- .ilrV(Xtall)
    Xwideilr <- tall2wide(Xilr, I, J, K)

    J <- J-1
    Xarrayilr <- array(Xwideilr, c(I, J, K))

    ## centering the compositions
    Xarrayilr <- do3Scale(Xarrayilr, center=center, center.mode=center.mode, scale=scale, scale.mode=scale.mode)
    Xwideilr <- unfold(Xarrayilr)

    Ahat<-matrix(0,dim(X)[1],ncomp)

    #define num. of outliers the algorithm should resists
    h <- round(alpha * dim(X)[1])

    ## Step 1 RobPCA XA
    outrobpca <- PcaHubert(Xwideilr, k=ncomp.rpca, kmax=ncol(Xwideilr), alpha=alpha, mcd=FALSE)
    Hset <- sort(sort(outrobpca@od, index.return=TRUE)$ix[1:h])
    Xhat <-Xwideilr[Hset,]
    fitprev <- 0
    changeFit <- 1 + conv
    iter <- 0
    while (changeFit > conv & iter <= robiter) {
        iter <- iter+1

        ## Step 2 - PARAFAC analysis
        modelPar <- cp_als(Xhat, h, J, K, ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
        Ah <- modelPar$A
        Bh <- modelPar$B
        Ch <- modelPar$C

        KR <- krp(Ch, Bh)  # Khatri-Rao product
        for(i in 1:dim(X)[1]) {
            vJKx1 <- matrix(Xarrayilr[i, ,], 1, J*K)
            Ahat[i,] <- .pinv(KR) %*% t(vJKx1)
        }
        Xfit <- Ahat %*%t (KR)

        ## Step 4  Computation of the rd
        rdsq <- apply((Xwideilr-Xfit)^2, 1, sum)
        rd <- sqrt(rdsq)
        Hset <- sort(sort(rdsq, index.return=TRUE)$ix[1:h])
        fit <- sum(rdsq[Hset])
        Xhat <- Xwideilr[Hset,]#Xunf

        ## Step 5  Fit of the model
        if (fitprev == 0) {
            changeFit <- 1 + conv
        } else {
            changeFit <- abs(fit-fitprev)/fitprev
        }
        fitprev <- fit
    }

    ## Reweighting
    cutoffOD <- .cutoff.rd(rd, h)
    out.rd <- rd
    flag <- (out.rd <= cutoffOD)
    Xflag <- Xarrayilr[flag,,] #X
    dim <- dim(Xflag)[1]
    Xflag_unf <- matrix(Xflag,dim,J*K)
    modelParrew <- cp_als(Xflag_unf, dim, J, K, ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
    Arew <- modelParrew$A
    Brew <- modelParrew$B
    Crew <- modelParrew$C

    KRrew <- krp(Crew, Brew)  # Khatri-Rao product
    Arew <- matrix(0, dim(X)[1], ncomp)
    for(i in 1:dim(X)[1])
    {
        vJKx1 <- matrix(Xarrayilr[i, ,], 1, J*K)#X
        Arew[i,] <- .pinv(KRrew) %*%t (vJKx1)
    }
    Xfitrew <- Arew %*% t(KRrew)
    odsqrew <- apply((Xwideilr-Xfitrew)^2, 1, sum) #Xunf
    out.rd <- sqrt(odsqrew)
    cutoffOD <- .cutoff.rd(out.rd, h)
    flag <- out.rd <= cutoffOD
    Xhat.rew <- array(Xfitrew, c(I, J, K))

    for(i in 1:ncomp) {
        Arew[,i] <- Arew[,i] * norm(as.matrix(Brew[,i]), type="F") * norm(as.matrix(Crew[,i]), type="F")
        Brew[,i] <- Brew[,i] / norm(as.matrix(Brew[,i]), type="F")
        Crew[,i] <- Crew[,i] / norm(as.matrix(Crew[,i]), type="F")
    }

    ## Back-transformation of loadings to clr
    V <- matrix(0, nrow = J+1, ncol = J)
    for (i in 1:ncol(V)) {
      V[1:i, i] <- 1/i
      V[i + 1, i] <- (-1)
      V[, i] <- V[, i] * sqrt(i/(i + 1))
    }
    Bclr <- V %*% Brew

    out.sd <- .cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)

    ## dimnames back
    nfac <- paste("F", 1:ncomp, sep="")
    dimnames(Arew) <- list(dn[[1]], nfac)
    dimnames(Brew) <- list(NULL, nfac)
    dimnames(Bclr) <- list(dn[[2]], nfac)
    dimnames(Crew) <- list(dn[[3]], nfac)
    dimnames(Xhat.rew) <- list(dn[[1]], NULL, dn[[3]])
    names(flag) <- dn[[1]]
    names(out.rd) <- dn[[1]]

    res <- list(fit=fit, fp=modelPar$fp, A=Arew, B=Brew, Bclr=Bclr, C=Crew, Zhat=Xhat.rew, const=modelParrew$const,
            flag=flag, Hset=Hset, iter=iter, alpha=alpha,
            rd=out.rd, cutoff.rd=cutoffOD, sd=out.sd$sd, cutoff.sd=out.sd$cutoff.sd,
            pcaobj=outrobpca,
            robust=TRUE, ilr=TRUE)

    class(res) <- "parafac"
    res
}

## - dd     = distance-distance plot
## - comp   = paired component plot for a single mode
##
plot.parafac <- function(x, which=c("dd", "comp", "percomp", "allcomp", "all"), ask = (which=="all" && dev.interactive(TRUE)), id.n, ...)
{
    which <- match.arg(which)
    op <- if(ask) par(ask = TRUE) else list()
    on.exit(par(op))

    if((which == "all" || which == "dd")) {
        ret <- .ddplot(x, id.n=id.n, ...)       # distance-distance plot
    }

    if((which == "all" || which == "comp")) {
        ret <- .compplot.parafac(x, ...)        # paired components plot
    }

    if((which == "all" || which == "percomp")) {
        ret <- .percompplot.parafac(x, ...)       # per-component plot
    }

    if((which == "all" || which == "allcomp")) {
        ret <- .allcompplot(x, ...)       # per-component plot
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

    cat("\nPARAFAC analysis with ", ncomp, " components.\nFit value:", round(x$fp,2), "%\n")
    msg <- ""
    if(x$robust)
        msg <- paste(msg, "Robust", sep="")
    if(x$ilr){
        if(nchar(msg) > 0)
            msg <- paste(msg, ", ", sep="")
        msg <- paste(msg, "ilr-transformed", "\n", sep="")
    }
    cat(msg, "\n")

    invisible(x)
}

##
## Adapted from package ThreeWay
##
##  n, m and p - Number of A-, B- and C-mode entities respectively
##  r - Number of extracted components
##
##  ort1, ort2 and ort3 - Type of constraints on A, B and C respectively
##      (1 for no constraints, 2 for orthogonality constraints, 3 for zero correlations constraints)
##
.CPfunc <- function (X, n, m, p, r, ort1, ort2, ort3, start, conv, maxit, A, B, C, trace=FALSE)
{
    X <- as.matrix(X)

    ftiter <- matrix(0, maxit/10, 2)
    mintripcos <- 0
    ssx <- sum(X^2)

    if(start == 0) {
        if(n >= r) {
            AUT <- eigen(X %*% t(X))
            A <- AUT$vectors[, 1:r]
        } else {
            A <- orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
            A <- A[1:n, ]
        }

        Z <- permute(X, n, m, p)
        if (m >= r) {
            AUT <- eigen(Z %*% t(Z))
            B <- AUT$vectors[, 1:r]
        } else {
            B <- orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
            B <- B[1:m, ]
        }

        Z <- permute(Z, m, p, n)
        if (p >= r) {
            AUT <- eigen(Z %*% t(Z))
            C <- AUT$vectors[, 1:r]
        } else {
            C <- orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
            C <- C[1:p, ]
        }
    }

    if(start == 1) {
        if (n >= r) {
            A <- orth(matrix(runif(n * r, 0, 1), nrow = n) - 0.5)
        } else {
            A <- orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
            A <- A[1:n, ]
        }
        if (m >= r) {
            B <- orth(matrix(runif(m * r, 0, 1), nrow = m) - 0.5)
        } else {
            B <- orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
            B <- B[1:m, ]
        }
        if (p >= r) {
            C <- orth(matrix(runif(p * r, 0, 1), nrow = p) - 0.5)
        } else {
            C <- orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
            C <- C[1:p, ]
        }
    }

    H <- matrix(0, r, r^2)
    for(ii in 1:r) {
        H[ii, (ii - 1) * r + ii] <- 1
    }
    H1 <- permute(H, r, r, r)
    H1 <- permute(B %*% H1, m, r, r)
    H1 <- permute(C %*% H1, p, r, m)
    f <- sum((X - A %*% H1)^2)
    if(trace)
        cat(paste("Candecomp/Parafac function value at Start is ", f, sep = " "), fill = TRUE)

    fold <- f + 2 * conv * f
    iter <- 0
    BB <- t(B) %*% B
    CC <- t(C) %*% C

    while((fold - f > conv * f | iter < 2) & f > conv^2 & iter < maxit) {
        fold <- f
        Z1 <- permute(X, n, m, p)
        Z1 <- permute(t(B) %*% Z1, r, p, n)
        Z1 <- permute(t(C) %*% Z1, r, n, r)
        XF <- Z1 %*% t(H)
        if (ort1 == 1) {
            FF <- BB * CC
            A <- XF %*% solve(FF)
        }
        if (ort1 == 2) {
            SVD <- svd(XF)
            A <- SVD$u %*% t(SVD$v)
        }
        if (ort1 == 3) {
            FF <- BB * CC
            SVD <- svd(XF - matrix(1, n, 1) %*% apply(XF, 2, mean))
            A <- SVD$u %*% t(SVD$v) + matrix(1, n, 1) %*% apply(XF, 2, mean) %*% solve(FF)
        }

        AA <- t(A) %*% A
        Z <- permute(X, n, m, p)
        Z1 <- permute(Z, m, p, n)
        Z1 <- permute(t(C) %*% Z1, r, n, m)
        Z1 <- permute(t(A) %*% Z1, r, m, r)
        XF <- Z1 %*% t(H)
        if (ort2 == 1) {
            FF <- AA * CC
            B <- XF %*% solve(FF)
        }
        if (ort2 == 2) {
            SVD <- svd(XF)
            B <- SVD$u %*% t(SVD$v)
        }

        if (ort2 == 3) {
            FF <- AA * CC
            SVD <- svd(XF - matrix(1, m, 1) %*% apply(XF, 2, mean))
            B <- SVD$u %*% t(SVD$v) + matrix(1, m, 1) %*% apply(XF, 2, mean) %*% solve(FF)
        }

        BB <- t(B) %*% B
        Z <- permute(Z, m, p, n)
        Z1 <- permute(Z, p, n, m)
        Z1 <- permute(t(A) %*% Z1, r, m, p)
        Z1 <- permute(t(B) %*% Z1, r, p, r)
        XF <- Z1 %*% t(H)
        if (ort3 == 1) {
            FF <- AA * BB
            C <- XF %*% solve(FF)
        }
        if (ort3 == 2)  {
            SVD <- svd(XF)
            C <- SVD$u %*% t(SVD$v)
        }
        if (ort3 == 3) {
            FF <- AA * BB
            SVD <- svd(XF - matrix(1, p, 1) %*% apply(XF, 2, mean))
            C <- SVD$u %*% t(SVD$v) + matrix(1, p, 1) %*% apply(XF, 2, mean) %*% solve(FF)
        }

        CC <- t(C) %*% C
        if (ort3 == 1) {
            f <- ssx - mtrace(CC %*% FF)
        } else {
            H1 <- permute(H, r, r, r)
            H1 <- permute(B %*% H1, m, r, r)
            H1 <- permute(C %*% H1, p, r, m)
            f <- sum((X - A %*% H1)^2)
        }

        iter <- iter + 1
        if ((iter%%10) == 0) {
            tripcos <- min(congruence(A, A) * congruence(B, B) * congruence(C, C))
            if (iter == 10)
                mintripcos <- tripcos
            if (tripcos < mintripcos)
                mintripcos <-  tripcos
            if ((iter%%1000) == 0 & trace)
                cat(paste("Minimal Triple cosine =", tripcos), fill = TRUE)

            ftiter[iter/10, ] <- c(f, tripcos)
        }

        if ((iter%%50) == 0 & trace) {
            cat(paste("f=", f, "after", iter, "iters; diff.=", (fold - f), sep = " "), fill = TRUE)
        }
    }

    ftiter <- ftiter[1:iter/10, ]           # take only the first iter/10 rows
    fp <- 100 - 100 * f/ssx
    tripcos <- min(congruence(A, A) * congruence(B, B) * congruence(C, C))
    names(tripcos) <- c("Minimal triple cosine")

    if (iter < 10) {
        mintripcos <- tripcos
    }

    if(trace) {
        cat(paste("Candecomp/Parafac function value is", f, "after",
            iter, "iterations", sep = " "), fill = TRUE)
        cat(paste("Fit percentage is", fp, "%", sep = " "), fill = TRUE)
    }

    out <- list(A=A, B=B, C=C, f=f, fp=fp, iter=iter, tripcos=tripcos, mintripcos=mintripcos, ftiter=ftiter)
    out
}
