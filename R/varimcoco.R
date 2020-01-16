## Addapted from ThreeWay
##
##  - A, B, C and H are the component matrices and the core array from a Tucker3 solution
##  - wa_rel, wb_rel and wc_rel are the requested relative weights for the three modes.
##      Default to all zero
##  - rot1, rot2 and rot3 tell on which mode to rotate (0/1)
##  - nanal is the number of soultuins to search for, not used currently by the caller
#3
.varimcoco <- function(A, B, C, H, wa_rel, wb_rel, wc_rel, rot1=1, rot2=1, rot3=1, nanal=5, trace=FALSE)
{
    conv <- 1e-6

    n <- nrow(A); r1 <- ncol(A)
    m <- nrow(B); r2 <- ncol(B)
    p <- nrow(C); r3 <- ncol(C)

    if(trace)
        cat(paste("Results from ", nanal, " random starts"), fill=TRUE)

    ## CHECK INPUT SPECIFICATION
    sz <- dim(H)
    if(sz[1] != r1 | sz[2] != r2*r3) {
	   stop("Error: sizes of A,B and C incompatible with H!")
    }

    gam1 <- gam2 <- gam3 <- gama <- gamb <- gamc <- 1

    ## Compute Natural Weights for core rotation:
    ## ensure natural proportion for three core terms
    w1=1/r2/r3
    w2=1/r1/r3
    w3=1/r1/r2
    ww=(w1+w2+w3)
    w1=w1/ww						# normalizes weights to sum of 1
    w2=w2/ww
    w3=w3/ww
    ssh <- sum(H^2)^2 / r1/r2/r3		# mean sq of squared core elements
    								    # makes contribution of qmax on uniform H equal to 1
    w1=w1/ssh
    w2=w2/ssh
    w3=w3/ssh
    v1=gam1*w1/(w1+w2+w3)
    v2=gam2*w2/(w1+w2+w3)
    v3=gam3*w3/(w1+w2+w3)

    # Compute WEIGHTS for columnwise orthonormal component matrices
    # relative weights times natural weights (these are "natural" given columnwise orthonormality!)
    wa=n/r1
    wb=m/r2
    wc=p/r3
    wa=wa_rel*wa
    wb=wb_rel*wb
    wc=wc_rel*wc

    # START LOOP FOR .. RUNS OF ALGORITHM
    Ssup <- matrix(NA, nrow=r1*r1, ncol=nanal)
    Tsup <- matrix(NA, nrow=r2*r2, ncol=nanal)
    Usup <- matrix(NA, nrow=r3*r3, ncol=nanal)

    func <- 1:nanal
    for(ana in 1:nanal) {
        if(rot1 == 1){
            set.seed(ana^2, kind=NULL, normal.kind=NULL)
            S <- orth(matrix(runif(r1*r1, 0, 1), r1, r1) - 0.5)
        } else {
            S <- diag(1, nrow=r1)
        }

        if (rot2==1){
            set.seed(ana^2, kind=NULL, normal.kind=NULL)
            T <- orth(matrix(runif(r2*r2, 0, 1), r2, r2) - 0.5)
    	} else{
    		T <- diag(1, nrow=r2)
    	}
        if(rot3 == 1){
            set.seed(ana^2, kind=NULL, normal.kind=NULL)
            U <- orth(matrix(runif(r3*r3, 0, 1), r3, r3) - 0.5)
        } else{
            U <- diag(1, nrow=r3)
        }

        ## Initialize Y, AS, BT and CU
        ## multiply by weights ^ 0.25
        Y  <- (w1 + w2 + w3) ^ 0.25 * S %*% H %*% kronecker(t(U), t(T))		
        AS <- wa ^ 0.25 * A %*% t(S)								
        BT <- wb ^ 0.25 * B %*% t(T)								
        CU <- wc ^ 0.25 * C %*% t(U)			

    ## EVALUATE FUNCTION VALUE AT START
    ## three-way orthomax part	
    	f1=sum((t(Y)*t(Y))^2)-v1/(r2*r3)*sum((colSums(t(Y)*t(Y)))^2)
    	Y <- permute(Y, r1, r2, r3)						# Y of order r2 x r3r
    	f1=f1-v2/(r1*r3)*sum((colSums(t(Y)*t(Y)))^2)
    	Y <- permute(Y, r2, r3, r1)						# Y of order r3 x r1r2
    	f1=f1-v3/(r1*r2)*sum((colSums(t(Y)*t(Y)))^2)
    	Y=permute(Y, r3, r1, r2)						# Y of order r1 x r2r3

    	# orthomax for component matrices
    	L=AS
    	L=L*L
    	f2a=sum(L^2)-gama/n*sum((colSums(L))^2)
    	L=BT
    	L=L*L
    	f2b=sum(L^2)-gamb/m*sum((colSums(L))^2)
    	L=CU
    	L=L*L
    	f2c=sum(L^2)-gamc/p*sum((colSums(L))^2)
    	f2=f2a+f2b+f2c
    	f=f1+f2
		
        ## START OF ITERATIVE PART
    	iter=0
    	fold=f-2*conv*abs(f)
    	while ((abs(f-fold))>(conv*abs(f))){  			
    		iter=iter+1
    		fold=f
    		
    		# UPDATING S
    		if(rot1 == 1){
    			# orthomax of     Y' = (w1+w2+w3)^.25 H_F ' with gamma = v1
    			# + orthomax of   AS = wa^.25A 		    with gamma = gama
           		ORTHMAXs=orthmax2(t(Y),AS,v1,gama,conv)			
    			AS=ORTHMAXs$B2
    			K=ORTHMAXs$B1
    			S=t(ORTHMAXs$T)%*%S
    		}
    		if(rot1 == 0){
    			K=t(Y)
    		}
    		Y <- permute(t(K[1:(r2*r3), ]), r1, r2, r3)				# Y of order r2 x r3r1

    		# UPDATING T
    		if (rot2==1){
    			# orthom of   Y' = (w1+w2+w3)^.25 H_F ' with gamma=v2
    			#+ orthom of BT = wb^.25B              with gamma=gamb
                ORTHMAXt=orthmax2(t(Y),BT,v2,gamb,conv)			
    			BT=ORTHMAXt$B2
    			K=ORTHMAXt$B1
    			T=t(ORTHMAXt$T)%*%T
    		}
    		if (rot2==0){
    			K=t(Y)
    		}
    		Y <- permute(t(K[1:(r1*r3), ]), r2, r3, r1)		 		# Y of order r3 x r1r2

    		# UPDATING U
    		if (rot3==1){
    			# orthom of   Y' = (w1+w2+w3)^.25 H_F ' with gamma=v3
    			# + orthom of CU = wc^.25C              with gamma=gamc
    			ORTHMAXu=orthmax2(t(Y),CU,v3,gamc,conv)	
    			CU=ORTHMAXu$B2
    			K=ORTHMAXu$B1
    			U=t(ORTHMAXu$T)%*%U
    		}
    		if (rot3==0){
    			K=t(Y)
    		}
    		Y <- permute(t(K[1:(r1*r2), ]), r3, r1, r2)	 		# Y of order r1 x r2r3

    		# EVALUATE ORTHCOCO FUNCTION:
    		# three-way orthomax part
    		f1=sum((t(Y)*t(Y))^2)-v1/(r2*r3)*sum((colSums(t(Y)*t(Y)))^2)
    		
    		Y=permute(Y, r1, r2, r3)							# Y of order r2 x r3r1
    		f1=f1-v2/(r1*r3)*sum((colSums(t(Y)*t(Y)))^2)
    		Y=permute(Y, r2, r3, r1)							# Y of order r3 x r1r2
    		f1=f1-v3/(r1*r2)*sum((colSums(t(Y)*t(Y)))^2)
    		Y=permute(Y, r3, r1, r2)							# Y of order r1 x r2r3
    		# orthomax for component matrices
    		L=AS
    		L=L*L
    		f2a=sum(L^2)-gama/n*sum((colSums(L))^2)
    		L=BT
    		L=L*L
    		f2b=sum(L^2)-gamb/m*sum((colSums(L))^2)
    		L=CU
    		L=L*L
    		f2c=sum(L^2)-gamc/p*sum((colSums(L))^2)
    		f2=f2a+f2b+f2c
    		f=f1+f2
    	}

        if(trace)
    	   cat(paste("Run no. ", ana, " f = ", round(f, digits=3), "(core: ", round(f1, digits=3), "; A: ", round(f2a, digits=3), "; B: ", round(f2b, digits=3), "; C: ", round(f2c, digits=3), "), ", iter, " iters"), fill=TRUE)

    	func[ana] <- f
    	Ssup[, ana] <- S
    	Tsup[, ana] <- T
    	Usup[, ana] <- U
    }

    ## FINISH LOOP FOR ANALYSES

    ## REPORT BEST SOLUTION
    f <- max(func)
    maxi <- which.max(func)					

    ## Rotation matrices
    S <- matrix(Ssup[, maxi], nrow(S), ncol(S))
    T <- matrix(Tsup[, maxi], nrow(T), ncol(T))
    U <- matrix(Usup[, maxi], nrow(U), ncol(U))

    if(trace)
        cat(paste("Three-way orthomax function value for best solution is ", round(f, digits=4)), fill=TRUE)

    ##  Rotated core:
    AS <- A %*% t(S)
    BT <- B %*% t(T)
    CU <- C %*% t(U)

    # reflect to positive sums in A,B and C
    sg=sign(colSums(AS))
    for (i in 1:length(sg)){
    	if(sg[i]==0){
    		sg[i]=1
    	}
    }

    S=diag(sg,nrow=r1)%*%S
    AS=AS%*%diag(sg,nrow=r1)

    sg=sign(colSums(BT))
    for (i in 1:length(sg)){
    	if (sg[i]==0){
    		sg[i]=1
    	}
    }
    T=diag(sg,nrow=r2)%*%T
    BT=BT%*%diag(sg,nrow=r2)

    sg=sign(colSums(CU))
    for (i in 1:length(sg)){
    	if (sg[i]==0){
    		sg[i]=1
    	}
    }
    U=diag(sg,nrow=r3)%*%U
    CU=CU%*%diag(sg,nrow=r3)

    ## compute core
    K=S%*%H%*%kronecker(t(U),t(T))

    ## recompute final orthomax values
    L=AS
    L=L*L
    f2a=(sum(L^2)-gama/n*sum((colSums(L))^2))*n/r1
    L=BT
    L=L*L
    f2b=(sum(L^2)-gamb/m*sum((colSums(L))^2))*m/r2
    L=CU
    L=L*L
    f2c=(sum(L^2)-gamc/p*sum((colSums(L))^2))*p/r3

    ## three-way orthomax part
    Y=((w1+w2+w3)^.25)*K
    f1=sum((t(Y)*t(Y))^2)-v1/(r2*r3)*sum((colSums(t(Y)*t(Y)))^2)
    Y <- permute(Y, r1, r2, r3)						# Y of order r2 x r3r1
    f1=f1-v2/(r1*r3)*sum((colSums(t(Y)*t(Y)))^2)
    Y <- permute(Y, r2, r3, r1)						# Y of order r3 x r1r2
    f1=f1-v3/(r1*r2)*sum((colSums(t(Y)*t(Y)))^2)
    Y <- permute(Y, r3, r1, r2)						# Y of order r1 x r2r3

    if(trace)
    {
        cat("Varimax values of core and AS, BT and CU (all based on natural weights)", fill=TRUE)
        cat("                Core             A               B               C ",fill=TRUE)
        cat(paste("              ",round(f1,digits=3),"          ",round(f2a,digits=3),"         ",round(f2b,digits=3),"          ",round(f2c,digits=3)),fill=TRUE)
        cat("These simplicity values are based on 'natural' weights and therefore comparable across matrices",fill=TRUE)
        cat("When multiplied by the relative weights, they give the contribution to the overall simplicity value",fill=TRUE)
        cat("(they are I^2/p, J^2/q or K^2/r, resp., times the sum of the variances of squared values)",fill=TRUE)		
    }

    out <- list(AS=AS, BT=BT, CU=CU, K=K, S=S, T=T, U=U, f=f, f1=f1, f2a=f2a, f2b=f2b, f2c=f2c, func=func)

    return(out)	
}
