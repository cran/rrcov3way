t3_als <- function(X, I, J, K, P, Q, R, start, conv=1e-10)
{

    ssx <- sum(X^2)

    if(!is.list(start) && length(start) != 1)
        stop("'start' must be either a list with elements A, B, C and GA, or a single character - one of 'random' or 'svd'!")

    if(!is.list(start) && !(start %in% c("random", "svd")))
        stop("'start' must be either a list with elements A, B, C and GA, or one of 'random' or 'svd'!")

    if(is.list(start)) {
        if(length(start) != 4 | sum(names(start) %in% c("A", "B", "C", "GA")) != 4)
            stop("'start' must be a list with elements A, B, C and GA!")
        if(!all(unlist(lapply(start, is.numeric))))
            stop("'start' must be a list containing 4 numeric matrices!")
    }

    tt <- if(is.list(start))
              ThreeWay::T3funcrep(X, I, J, K, P, Q, R, start=2, conv=conv, A=start$A, B=start$B, C=start$C, H=start$GA)
          else
              ThreeWay::T3funcrep(X, I, J, K, P, Q, R, start=if(start=="svd") 0 else 1, conv=conv)

    ret <- list(f=tt$f, fp=tt$fp, ss=ssx,
                A=tt$A, B=tt$B, C=tt$C, GA=tt$H,
                La=tt$La, Lb=tt$Lb, Lc=tt$Lc,
                iter=tt$iter)
}
