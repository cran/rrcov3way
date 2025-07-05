unimcd <- function(y, quan){
    out <- list()
    ncas <- length(y)
    len <- ncas-quan+1

    if(len == 1){
        out$tmcd <- mean(y)
        out$smcd <- sqrt(var(y))
    } else {
        ay <- c()
        I <- order(y)
        y <- y[I]
        ay[1] <- sum(y[1:quan])
        for(samp in 2:len){
            ay[samp]<-ay[samp-1]-y[samp-1]+y[samp+quan-1]
        }
        ay2<-ay^2/quan
        sq<-c()
        sq[1]<-sum(y[1:quan]^2)-ay2[1]
        for(samp in 2:len){
            sq[samp]<-sq[samp-1]-y[samp-1]^2+y[samp+quan-1]^2-ay2[samp]+ay2[samp-1]
        }
        sqmin<-min(sq)
        Isq<-order(sq)
        ndup<-sum(sq == sqmin)
        ii<-Isq[1:ndup]
        slutn<-c()
        slutn[1:ndup]<-ay[ii]
        initmean<-slutn[floor((ndup+1)/2)]/quan
        initcov<-sqmin/(quan-1)
        res<-(y-initmean)^2/initcov
        sortres<-sort(res)
        factor<-sortres[quan]/qchisq(quan/ncas,1)
        initcov<-factor*initcov
        res<-(y-initmean)^2/initcov
        quantile<-qchisq(0.975,1)
        out$weights<-(res<quantile)
        out$tmcd<-sum(y*out$weights)/sum(out$weights)
        out$smcd<-sqrt(sum((y-out$tmcd)^2*out$weights)/(sum(out$weights)-1))
        Iinv<-order(I)
        out$weights<-out$weights[Iinv]
    }
    return(out)
}
