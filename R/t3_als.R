##
##  X is an n x m x p array
##
##  - the function ALS() expects a n x m*p matrix
##
.TK <- function(X, P=2, Q=2, R=2, conv=1e-10, crit=0.975)
{
    di <- dim(X)
    n <- di[1]
    m <- di[2]
    p <- di[3]
    dn <- dimnames(X)

    Y <- unfold(X)     # default mode is "A", returns n x m*p matrix
    tt <- .ALS(Y, n, m, p, P, Q, R, conv)

    Xfit <- tt$A %*% tt$H %*% t(kronecker(tt$C, tt$B))
    out.Xhat <- array(Xfit, c(n, m, p))
    odsq <- apply((Y-Xfit)^2, 1, sum)
    RD <- sqrt(odsq)

##    critRD <- (mean(RD^(2/3)) + sd(RD^(2/3)) * qnorm(crit))^(3/2)
    critRD <- .cutoff.rd(RD, crit=crit, robust=FALSE)

    flag <- (RD <= critRD)

    ## dimnames back
    dimnames(tt$A) <- list(dn[[1]], paste("F", 1:P, sep=""))
    dimnames(tt$B) <- list(dn[[2]], paste("F", 1:Q, sep=""))
    dimnames(tt$C) <- list(dn[[3]], paste("F", 1:R, sep=""))
    dimnames(tt$H) <- list(paste("F", 1:P, sep=""), paste("F", 1:(Q*R), sep=""))
    names(RD) <- dn[[1]]
    names(flag) <- dn[[1]]

    ret <- list(fit=tt$f, fp=tt$fp,
                A=tt$A, B=tt$B, C=tt$C, GA=tt$H,
                La=tt$La, Lb=tt$Lb, Lc=tt$Lc,
                iter=tt$iter, flag=flag, rd=RD, cutoff.rd=critRD)
}

##
##  Tucker3
##  Adapted from T3funcrep() in package ThreeWay
##
##  Alternating least squares for Tucker3 model
##
.ALS <- function(X, n, m, p, r1, r2, r3, conv)
{
    X <- as.matrix(X)

	## initialize A, B and C
	ss <- sum(X^2)
	dys <- 0

	## rational starts via eigendecompositions
	EIG <- eigen(X %*% t(X))
	A <- EIG$vectors[, 1:r1]
	Z <- permute(X, n, m, p)		# yields m x p x n array: n,m,p ==> m,p,n
	EIG <- eigen(Z %*% t(Z))
	B <- EIG$vectors[, 1:r2]
	Z <- permute(Z, m, p, n)		# yields p x n x m array: m,p,n ==> p,n,m
	EIG <- eigen(Z %*% t(Z))
	C <- EIG$vectors[, 1:r3]

	## Update Core
	Z <- permute(t(A) %*% X, r1, m, p)
	Z <- permute(t(B) %*% Z, r2, p, r1)
	H <- permute(t(C) %*% Z, r3, r1, r2)

	## Evaluate f
    f <- ss - sum(H^2)


	iter <- 0
	fold <- f + 2*conv*f
	while(fold-f > f*conv)
    {
		iter <- iter+1
		fold <- f

		## update A (Z = X*C'x B' - GS Z*Z'*A)
		Z <- permute(X,n,m,p)
		Z <- permute(t(B)%*%Z,r2,p,n)
		Z <- permute(t(C)%*%Z,r3,n,r2)			 # yields n x r2 x r3 array
		A <- qr.Q(qr(Z%*%(t(Z)%*%A)),complete=FALSE)
	
		# update B
		Z <- permute(X, n, m, p)
		Z <- permute(Z, m, p, n)
		Z <- permute(t(C) %*% Z, r3, n, m)
		Z <- permute(t(A) %*% Z, r1, m, r3)		 # yields m x r3 x r1 array
		B <- qr.Q(qr(Z %*% (t(Z) %*% B)), complete=FALSE)

		## update C
		Z <- permute(t(A) %*% X, r1, m, p)
		Z <- permute(t(B) %*% Z, r2, p, r1)			 # yields p x r1 x r2 array
		C <- qr.Q(qr(Z %*% (t(Z) %*% C)), complete=FALSE)

		## Update Core
		Z <- permute(t(A) %*% X, r1, m, p)
		Z <- permute(t(B) %*% Z, r2, p, r1)
		H <- permute(t(C) %*% Z, r3, r1, r2)

		## Evaluate f
		f <- ss - sum(H^2)
	}

    ss <- sum(X^2)
    fp <- 100*(ss-f)/ss

    ## compute "intrinsic eigenvalues"
    ## eigenvalues for A-mode:
    La <- H %*% t(H)
    Y <- permute(H, r1, r2, r3)
    Lb <- Y %*% t(Y)
    Y <- permute(Y, r2, r3, r1)
    Lc <- Y %*% t(Y)

    list(A=A, B=B, C=C, H=H,
         f=f, fp=fp,
         iter=iter,
         La=La, Lb=Lb, Lc=Lc)
}
