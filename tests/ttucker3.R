## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(rrcov3way))

## Example with rotation
data(elind)
t3 <- Tucker3(elind, 3, 2, 2)
xout <- do3Rotate(t3, c(3, 3, 3), rotate=c("A", "B", "C"))
xout$vvalue

w <- c(NA, 0, 0.5, 1, 2.5, 3, 3.5, 4, 5, 10, Inf)
res <- matrix(NA, nrow=length(w), ncol=7)
for(i in 1:length(w))
{
    res[i, 1] <- res[i, 2] <- res[i, 3] <- w[i]
    if(is.na(w[i]))
        x <- do3Rotate(t3, rotate=c())      ## no rotation
    else if(is.finite(w[i]))
        x <- do3Rotate(t3, rep(w[i], 3))
    else
        x <- do3Rotate(t3, rep(1e18, 3))

    res[i, 4:7] <- round(x$vvalue,3)
}
colnames(res) <- c("A", "B", "C", "Core", "A", "B", "C")
rownames(res) <- 1:nrow(res)
res

## Example with the UNIDO Manufacturing value added data
data(va3way)
dim(va3way)

## Treat quickly and dirty the zeros in the data set (if any)
va3way[va3way==0] <- 0.001

(res <- Tucker3(va3way))
print(res$fit)
print(res$A)
print(res$B)
print(res$C)
print(res$rd)
print(res$cutoff.rd)

(res.r <- Tucker3(va3way, robust=TRUE))
print(res.r$fit)
print(res.r$A)
print(res.r$B)
print(res.r$C)
print(res.r$rd)
print(res$cutoff.rd)
print(res.r$Hset)
print(res.r$flag)

(res.c <- Tucker3(va3way, coda.transform="ilr"))
print(res.c$fit)
print(res.c$A)
print(res.c$B)
print(res.c$Bclr)
print(res.c$C)
print(res.c$rd)
print(res$cutoff.rd)

(res.rc <- Tucker3(va3way, robust=TRUE, coda.transform="ilr"))
print(res.rc$fit)
print(res.rc$A)
print(res.rc$B)
print(res.rc$Bclr)
print(res.rc$C)
print(res.rc$rd)
print(res$cutoff.rd)
print(res.rc$Hset)
print(res.rc$flag)
