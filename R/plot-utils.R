## Plot RD against SD (outlier map) for a paradac or a tucker 3 object
##
## it is assumed that the object has list elements A (A-mode loadings),
##     RD (residual distances) and SD (score distances).
##
## - robust=TRUE generates the robust outlier map, otherwise classical
## - crit is the quantile of the reference distribution
##
##  Distance-Distance Plot:
##  Plot the vector y=RD (residual distances) against
##  x=SD (score distances). Identify by a label the id.n
##  observations with largest RD. If id.n is not supplied, calculate
##  it as the number of observations larger than cutoff. Use cutoff
##  to draw a horisontal and a vertical line. Draw also a dotted line
##  with a slope 1.
.ddplot <- function(obj, crit=0.975, id.n, labs, ...)
{
    .label <- function(x, y, id.n=3, labs=seq_len(length(x)), ...) {
        if(id.n > 0) {
            xrange <- par("usr")
            xrange <- xrange[2] - xrange[1]
            n <- length(y)
            ind <- sort(y, index.return=TRUE)$ix
            ind <- ind[(n-id.n+1):n]
            text(x[ind] + xrange/30, y[ind], labs[ind], ...)
        }
    }


    RD <- obj$rd
    critRD <- obj$cutoff.rd
    SD <- obj$sd
    critSD <- obj$cutoff.sd
    A <- obj$A
    n <- length(RD)

    if(obj$robust)
    {
        xl <- "SD robust"
        yl <- "RD robust"
    }
    else
    {
        xl <- "SD non-robust"
        yl <- "RD non-robust"
    }

    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    } else
        id.n <- length(which(RD > critRD))


    if(missing(labs))
        labs <- rownames(A)
    if(is.null(labs))
        labs <- seq_len(nrow(A))
    xlim <- c(0, max(SD, critSD))
    xlim[2] <- xlim[2] + 0.1 * xlim[2]
    plot(SD, RD, xlab=xl, ylab=yl, xlim=xlim, ylim=c(0,max(RD, critRD)), cex.lab=1, type="p", ...)
    .label(SD, RD, id.n=id.n, labs=labs, cex=0.8, ...)
    abline(h=critRD, lty=2)
    abline(v=critSD, lty=2)

    return(invisible(list(SD=SD, RD=RD, critSD=critSD, critRD=critRD)))
}

##
##  Loadings plot for modes A, B and C of Tucker3 model
##
.compplot.tucker3 <- function (x, mode=c("A", "B", "C"), choices=1L:2L, xlim, ylim, arrows=TRUE, 
    cex.text.A=0.8, cex.text.B=0.8, x.off, y.off, ...)
{

    mode <- match.arg(mode)

    A <- x$A
    B <- if(x$coda.transform == "ilr") x$Bclr else x$B
    C <- x$C
    GA<- x$GA
    P <- dim(A)[2]
    Q <- dim(B)[2]
    R <- dim(C)[2]

    p <- if(mode=="A") P else if(mode=="B") Q else R

    if(length(choices) != 2 || min(choices) < 1 || max(choices) > p || choices[1] == choices[2]) {
        warning("Wrong components choosen! Components 1 and 2 will be used.")
        choices <- 1L:2L
    }
    comp1 <- choices[1]
    comp2 <- choices[2]

    ## Unfold the core array for B and C mode
    GG <- toArray(GA, P, Q, R)

    if(mode == "A") {
        GGG <- unfold(GG, mode="A")
        Fx <- kronecker(C,B) %*% t(GGG)
        qrB <- qr(Fx)                # The QR Decomposition of a Matrix B == Q %*% R
        Tx <- qr.R(qrB)
        tilde <- A %*% solve(Tx)
        lab <- rownames(A)
    } else if(mode == "B") {
        GGG <- unfold(GG, mode="B")
        Fx <- kronecker(C,A) %*% t(GGG)
        qrB <- qr(Fx)                # The QR Decomposition of a Matrix B == Q %*% R
        Tx <- qr.R(qrB)
        tilde <- B %*% solve(Tx)
        lab <- rownames(B)
    } else {
        GGG <- unfold(GG, mode="C")
        Fx <- kronecker(B,A) %*% t(GGG)
        qrB <- qr(Fx)                # The QR Decomposition of a Matrix B == Q %*% R
        Tx <- qr.R(qrB)
        tilde <- C %*% solve(Tx)
        lab <- rownames(C)
    }

    eps <- 0.1
    if(missing(xlim)) {
        xlim <- c(min(tilde[,comp1]), max(tilde[,comp1]))
        xlim[1] <- xlim[1] - abs(eps*xlim[1])
        xlim[2] <- xlim[2] + abs(eps*xlim[2])
        if(xlim[1] > 0) xlim[1] <- 0
        if(xlim[2] < 0) xlim[2] <- 0
    }
    if(missing(ylim)) {
        ylim <- c(min(tilde[,comp2]), max(tilde[,comp2]))
        ylim[1] <- ylim[1] - abs(eps*ylim[1])
        ylim[2] <- ylim[2] + abs(eps*ylim[2])
        if(ylim[1] > 0) ylim[1] <- 0
        if(ylim[2] < 0) ylim[2] <- 0
    }

    plot(tilde[,choices], type="n", xlab=paste("Axis", comp1), ylab=paste("Axis", comp2), xlim=xlim, ylim=ylim, cex=1.2, ...)
    abline(v=0, h=0, lty=2)
    
    tilde.text <- matrix(0, nrow=nrow(tilde), ncol=ncol(tilde))
    if(!missing(x.off))
        tilde.text[, 1] <- tilde.text[, 1] + x.off 
    if(!missing(y.off))
        tilde.text[, 2] <- tilde.text[, 2] + y.off 
        
    if(mode == "B" & arrows) {
        arrows(0, 0, tilde[,comp1], tilde[,comp2], code = 2, length = 0.09)
        text(tilde[,comp1] + tilde.text[, 1], tilde[,comp2] + tilde.text[, 2], lab, cex=cex.text.B, pos=1, ...)
    }else
        text(tilde[,comp1] + tilde.text[, 1], tilde[,comp2] + tilde.text[, 2], lab, cex=cex.text.A, col=4, ...)

    return(invisible(x))
}

##
##  Loadings plot for modes A, B and C of PARAFAC model
##
.compplot.parafac <- function (x, mode=c("A", "B", "C"), choices=1L:2L, xlim, ylim, arrows=TRUE, ...)
{
    mode <- match.arg(mode)

    A <- x$A
    B <- if(x$coda.transform=="ilr") x$Bclr else x$B
    C <- x$C
    ncomp <- x$ncomp

    if(length(choices) != 2 || min(choices) < 1 || max(choices) > ncomp || choices[1] == choices[2]) {
        warning("Wrong components choosen! Components 1 and 2 will be used.")
        choices <- 1L:2L
    }
    comp1 <- choices[1]
    comp2 <- choices[2]


    ## create a long, matricized array of diagonal matrices of size 'ncomp'
    D <- do.call(rbind, lapply(1:ncomp, function(x) diag(ncomp)))

    if(mode == "A") {
        ## orthonormalization MODE A
        WA <- kronecker(C, B) %*% D
        qrstr <- qr(WA)
        R <- qr.R(qrstr)                # T transformation matrix
        tilde  <- (A %*% R)
    } else if(mode == "B") {
        ## orthonormalization MODE B
        WB <- kronecker(C, A) %*% D
        qrstr <- qr(WB)
        R <- qr.R(qrstr)
        tilde <- B %*% R
        obl_B <- R
        coordB <- as.vector(obl_B)      # oblique original axes
    } else {
        ## orthonormalization MODE C
        WC <- kronecker(B, A) %*% D
        qrstr <- qr(WC)
        R <- qr.R(qrstr)
        tilde <- C %*% R
    }

    eps <- 0.1
    if(missing(xlim)) {
        xlim <- c(min(tilde[,comp1]), max(tilde[,comp1]))
        xlim[1] <- xlim[1] - abs(eps*xlim[1])
        xlim[2] <- xlim[2] + abs(eps*xlim[2])
        if(xlim[1] > 0) xlim[1] <- 0
        if(xlim[2] < 0) xlim[2] <- 0
    }
    if(missing(ylim)) {
        ylim <- c(min(tilde[,comp2]), max(tilde[,comp2]))
        ylim[1] <- ylim[1] - abs(eps*ylim[1])
        ylim[2] <- ylim[2] + abs(eps*ylim[2])
        if(ylim[1] > 0) ylim[1] <- 0
        if(ylim[2] < 0) ylim[2] <- 0
    }

    plot(tilde[,choices], type="n", xlim=xlim, ylim=ylim, xlab="First component", ylab="Second component", cex=0.8, ...)
    abline(v=0, h=0, lty = 2)
    
    if(mode == "B" & arrows) {
        arrows(0, 0, tilde[,comp1], tilde[,comp2], code = 2, length = 0.09)
        text(tilde[, comp1], tilde[, comp2], labels=rownames(tilde), cex=0.8, pos=1, ...)
    }else
        text(tilde[, comp1], tilde[, comp2], labels=rownames(tilde), cex=0.8, ...)

    return(invisible(x))
}

.allcompplot <- function (x, mode=c("C", "B", "A"), xlim, ylim, xlab, ylab, legend.position="topleft", points=TRUE, ...)
{
    mode <- match.arg(mode)

    C <- x$C
    if(mode == "A")
        C <- x$A
    else if(mode=="B")
        C <- x$B

    if(missing(ylim))
        ylim <- c(min(C), max(C))
    if(missing(xlab))
        xlab <- "Time"
    if(missing(ylab))
        ylab <- "Component score"
    time <- seq_len(nrow(C))
    names(time) <- rownames(C)
    plot(time, C[,1], ylim=ylim, xlab=xlab, ylab=ylab, type="n", xaxt="n", ...)
    for(i in seq_len(ncol(C)))
    {
        if(points)
            points(time, C[,i], pch=i, col=i, bg=i)
        lines(time, C[,i], col=i, lty=i)
    }

    abline(h=0, lty="dotted")
    axis(side=1, at=time, labels=names(time))
    myseq <- seq_len(ncol(C))
    if(!is.null(legend.position) && legend.position != "none")
        legend(legend.position, pch=myseq, col=myseq, pt.bg=myseq,
            legend=paste("Component", myseq), lty=myseq)
}

.percompplot.parafac <- function (x, comp=1, ...)
{
    stopifnot(comp <= x$ncomp)

    Ah <- x$A
    Bh <- if(x$coda.transform=="ilr") x$Bclr else x$B
    Ch <- x$C
    ncomp <- x$ncomp
    Anames <- rownames(Ah)
    Bnames <- rownames(Bh)
    Cnames <- rownames(Ch)


    x1 <- Ah[,comp]/sqrt(sum(Ah[,comp]^2))
    x2 <- Bh[,comp]/sqrt(sum(Bh[,comp]^2))
    x3 <- Ch[,comp]/sqrt(sum(Ch[,comp]^2))
##    x1 <- Ah[,comp]/sqrt(nrow(Ah))
##    x2 <- Bh[,comp]/sqrt(nrow(Bh))
##    x3 <- Ch[,comp]/sqrt(nrow(Ch))

    eps <- 0.1
    ylim1 <- min(x1,x2,x3)
    ylim2 <- max(x1,x2,x3)
    ylim <- c(ylim1 - ylim1 * eps, ylim2 + ylim2 * eps)

    plot(c(-2,0,2), c(ylim1,0,ylim2), type="n", ylab=paste("Component", comp), xlab= "", xaxt="n")
    text(1, x1, Anames,cex=0.5)
    text(-0.25, x3, Cnames,  cex=0.5)
    text(-1, x2, Bnames, cex=0.5)

    arrows(0, x1, 0.9, x1, code = 2, length = 0.09)
    arrows(0, x2, -0.9, x2,  code = 2, length = 0.09)
    abline(v=0)
    axis(1, at = c(-1, 0, 1),
    labels = c("Second Mode","Third Mode", "First Mode "), cex.axis=0.7, tick=FALSE)
##    axis(2, at = c(-0.5,0,0.5),labels = c(-0.5,0,0.5),cex.axis=0.9,tick=FALSE)

    return(invisible(x))
}


## Joint biplot (for Tucker 3)
.JBPlot <- function (x, alfa=.5, comp=1, cex.text.A=1.0, cex.text.B=1.0, rotate=c(1,1,1,1), A.x.off, A.y.off, B.x.off, B.y.off, ...)
{
    ##  cex.text.A and cex.text.B - text size of the A-mode and B-mode labels respectively
    ##  rotate <- c(A.x, A.y, B.x, B.y) - a vector with 4 elements which are 1 or -1
    ##  A.x.off, A.y.off: x- and y-offset of the labels of mode A
    ##  B.x.off, B.y.off: x- and y-offset of the labels of mode B (labels of the arrows)
     
    ## A is loadings matrix for first mode
    ## Bclr is loadings matrix for second mode only two components (clr transformation)
    ## C is loadings matrix for third mode
    ## G is wide unfolded core array

    ## alfa is a number [0,1]
    ## comp is the frontal slice of core array (1 is first frontal
    ##  slice of core array, 2 second , r-th ...)
    ##
    ## Warning if all elements of first C loading are negative we
    ##  need to change the sign of all the elements of first B loading

    A <- x$A
    B <- if(x$coda.transform == "ilr") x$Bclr else x$B
    C <- x$C
    GA<- x$GA

    I <- dim(A)[1]
    J <- dim(B)[1]
    K <- dim(C)[1]
    P <- dim(A)[2]
    Q <- dim(B)[2]
    R <- dim(C)[2]

    ## Select the frontal slice of core array
    GG <- toArray(GA, P, Q, R)
    GGG <- GG[,, comp]

    ## stretching or shrinking parameters
    k <- alfa
    ssa <- (I/J)^(.25)
    ssb <- (J/I)^(.25)

    ## coordinates
    SVDG <- svd(GGG)
    singv <- diag(SVDG$d)

    Atilde <- ssa*(A %*% SVDG$u %*% singv^(k))
    Btilde <- ssb*(B %*% SVDG$v %*% singv^((1-k)))

    Atilde[, 1] <- rotate[1] * Atilde[, 1]
    Atilde[, 2] <- rotate[2] * Atilde[, 2]
    Btilde[, 1] <- rotate[3] * Btilde[, 1]
    Btilde[, 2] <- rotate[4] * Btilde[, 2]

    if(!missing(A.x.off))
        Atilde[, 1] <- Atilde[, 1] + A.x.off 
    if(!missing(A.y.off))
        Atilde[, 2] <- Atilde[, 2] + A.y.off 

    Btilde.text <- matrix(0, nrow=nrow(Btilde), ncol=ncol(Btilde))
    if(!missing(B.x.off))
        Btilde.text[, 1] <- Btilde.text[, 1] + B.x.off 
    if(!missing(B.y.off))
        Btilde.text[, 2] <- Btilde.text[, 2] + B.y.off 

    ## plot
    ## Warning
    ## VT::09.01.2024 - reduce Btilde[,2] a bit, because the text of the arrows will be written below
    aa <- sum(sign(C[,1]))
    if(K == -aa) {
        plot(c(min(Atilde[,1],-Btilde[,1]), max(Atilde[,1],-Btilde[,1])),
             c(min(Atilde[,2], Btilde[,2]), max(Atilde[,2], Btilde[,2])),
                type="n", xlab="First axis", ylab="Second axis", ...)
        abline(v=0, h=0, lty=1)
        text(-Btilde[,1] + Btilde.text[,1], Btilde[,2] + Btilde.text[, 2], rownames(B), col=1, cex=cex.text.B, pos=1)
        arrows(0, 0, -Btilde[,1], Btilde[,2], code=2, length=0.09)
        text(Atilde[,1], Atilde[,2], rownames(A), col=4, cex=cex.text.A)
    } else {
        plot(c(min(Atilde[,1], Btilde[,1]), max(Atilde[,1], Btilde[,1])),
             c(min(Atilde[,2], Btilde[,2]-0.1), max(Atilde[,2],Btilde[,2])),
                type="n", xlab="First axis", ylab="Second axis", ...)
        abline(v=0, h=0, lty=1)
        text(Btilde[,1] + Btilde.text[,1], Btilde[,2] + Btilde.text[, 2], rownames(B), col=1, cex=cex.text.B, pos=1)
        arrows(0, 0, Btilde[,1], Btilde[,2], code=2, length=0.09)
        text(Atilde[,1], Atilde[,2], rownames(A), col=4, cex=cex.text.A)
    }

    return(invisible(x))
}

## Trajectory biplot (for Tucker 3)
.TJPlot <- function (x, choices, arrows=TRUE, longnames=TRUE, 
    cex.line=0.5, cex.text=0.6, cex.text.arrow=0.8, ...)
{
    
    ##  cex.line         <- 0.5      # trajectory line symbols
    ##  cex.text         <- 0.6      # trajectory text
    ##  cex.text.arrow   <- 0.8      # arrow text
    
    ## A is loadings matrix for first mode .
    ## Bclr is loadings matrix for second mode only two components (clr transformation).
    ## C is loadings matrix for third mode
    ## G is wide unfolded core array.

    A <- x$A
    B <- if(x$coda.transform == "ilr") x$Bclr else x$B
    C <- x$C
    GA<- x$GA

    I <- dim(A)[1]
    J <- dim(B)[1]
    K <- dim(C)[1]
    P <- dim(A)[2]
    Q <- dim(B)[2]
    R <- dim(C)[2]

    if(missing(choices))
        choices <- 1:I

    ssa <- (I/J)^(.25)
    ssb <- (J/I)^(.25)

    ## Make unfolded core array for B and C mode
    GG <- toArray(GA, P, Q, R)
    GB <- unfold(GG, mode="B")

    ## Make unfolded A x C label
    labA <- matrix(0, I, K)
    for(k in 1:K) {
        for(i in 1:I) {
            if(longnames)
                labA[i, k] <- paste(rownames(A)[i], rownames(C)[k], sep="x")
            else
                labA[i, k] <- paste(i, k, sep="-")
        }
    }
    rowlabACco <- as.vector(labA)

    ## The QR Decomposition of a Matrix B == Q %*% R
    qrB <- qr(B)
    Bco <- qr.Q(qrB)
    Tx <- qr.R(qrB)

    ## to enforce positive diagonals of R, and thereby
    ## get a unique factorisation. Is this correct?
    D <- diag(sign(diag(Tx)))

##  !!!!!!!!!!!!!!!!!!
##    Bco <- Bco %*% D
    Bco <- Bco %*% D/2
    Tx <- D %*% Tx

    ACco <- kronecker(C,A) %*% t(GB) %*% solve(Tx)

    ## Here I try to rescale the arrows - otherwise they are not visible
    Bco <- 0.8 * max(ACco)/max(Bco) * Bco
    
    ## collect the selected points for the choices objects
    cx <- c()
    for(j in choices) cx <- c(cx, seq(j, I*K, I))

    #PLOT
    
    ## VT::09.01.2024 - reduce Bco[,2] a bit, because the text of the arrows will be written below
    if(arrows) {
        x <- c(min(0, ACco[cx, 1], Bco[,1]), max(0, ACco[cx, 1], Bco[,1]))
        y <- c(min(0, ACco[cx, 2], Bco[,2]-0.15), max(0, ACco[cx, 2], Bco[,2]))
    } else {
        x <- c(min(ACco[cx, 1]), max(ACco[cx, 1]))
        y <- c(min(ACco[cx, 2]), max(ACco[cx, 2]))
    }
    plot(x, y, xlab="First trajectory axis", ylab="Second trajectory axis", type="n", ...)
    abline(v=0, h=0, lty = 2)

    if(arrows) {
        arrows(0, 0, Bco[,1], Bco[,2], code = 2, length = 0.09)
        text(Bco[,1], Bco[,2], rownames(B), col=1, cex=cex.text.arrow, pos=1)
    }

    for(i in choices) {
        aa <- seq(i, I*K, I)
##        print(i)
##        print(aa)
        lines(c(ACco[aa,1]), c(ACco[aa,2]), col=2, lty = 3, type = "o", cex=cex.line)
        text(ACco[aa[c(1,K)],1], ACco[aa[c(1,K)],2], rowlabACco[aa[c(1,K)]], col=1, cex=cex.text)
    }

    return(invisible(x))
}
