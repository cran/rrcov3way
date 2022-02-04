## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(rrcov3way))

## t3_als() ======================================================
data(elind)
a <- unfold(elind)
try(rrcov3way:::t3_als(a, 23, 6, 7, 2, 2, 2, start=c(1, 2, 3)))  # not a list, length != 1
try(rrcov3way:::t3_als(a, 23, 6, 7, 2, 2, 2, start=0))           # not a list, must be "random" or "svd"
try(rrcov3way:::t3_als(a, 23, 6, 7, 2, 2, 2, start="rand"))      # not a list, must be "random" or "svd"
try(rrcov3way:::t3_als(a, 23, 6, 7, 2, 2, 2,                     # if list, elements must be in (A, B, C, GA)
    start=list(A="A", B="B", C="C", D="GA")))
try(rrcov3way:::t3_als(a, 23, 6, 7, 2, 2, 2,                     # if list, elements must be numeric matrices
    start=list(A="A", B="B", C="C", GA="GA")))

set.seed(98765)

## Example with rotation
data(elind)
t3 <- Tucker3(elind, 3, 2, 2)

## IGNORE_RDIFF_BEGIN

xout <- do3Rotate(t3, c(3, 3, 3), rotate=c("A", "B", "C"))
xout$vvalue

w <- c(NA, 0, 0.5, 1, 2.5, 3, 3.5, 4, 5, 10, Inf)
res <- matrix(NA, nrow=length(w), ncol=7)
for(i in seq_along(w))
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

## IGNORE_RDIFF_END

colnames(res) <- c("A", "B", "C", "Core", "A", "B", "C")
rownames(res) <- seq_len(nrow(res))
res

## Example with the UNIDO Manufacturing value added data
data(va3way)
dim(va3way)

## Treat quickly and dirty the zeros in the data set (if any)
va3way[va3way==0] <- 0.001

set.seed(123456)

## Using robustness with clr transformation
try(res <- Tucker3(va3way, robust=TRUE, coda.transform="clr"))

## Rejected values of parameter 'crit'
try(res <- Tucker3(va3way, crit=c(1:10)))       # length different than 1
try(res <- Tucker3(va3way, crit=-1))            # crit non-positive
try(res <- Tucker3(va3way, crit=2))             # crit >= 1
res <- Tucker3(va3way, crit=0.2)                # crit < 0.5 --> crit=1-crit

##  Standard Tucker 3
(res <- Tucker3(va3way, center=TRUE, scale=TRUE))
print(res$fit)
print(res$A)
print(res$B)
print(res$C)
print(res$rd)
print(res$cutoff.rd)

## Robust Tucker 3
(res.r <- Tucker3(va3way, robust=TRUE, center=TRUE,
    scale=TRUE, trace=TRUE))
print(res.r$fit)
print(res.r$A)
print(res.r$B)
print(res.r$C)
print(res.r$rd)
print(res$cutoff.rd)
print(res.r$Hset)
print(res.r$flag)

## Compositional  Tucker 3
(res.c <- Tucker3(va3way, coda.transform="ilr"))
print(res.c$fit)
print(res.c$A)
print(res.c$B)
print(res.c$Bclr)
print(res.c$C)
print(res.c$rd)
print(res$cutoff.rd)

## Robust, compositional Tucker 3
(res.rc <- Tucker3(va3way, robust=TRUE, coda.transform="ilr",
    center=TRUE, scale=TRUE, trace=TRUE))
print(res.rc$fit)
print(res.rc$A)
print(res.rc$B)
print(res.rc$Bclr)
print(res.rc$C)
print(res.rc$rd)
print(res$cutoff.rd)
print(res.rc$Hset)
print(res.rc$flag)

##================================================================
##
##  Example with the TV data from ThreeWay: ratings of 15 American television
##  shows on 16 bipolar scales made by 30 students: 15x16x30 array
##  in A-mode, i.e. 15x450.
data(TV, package="ThreeWay")

## Transform to 3-way array, to pass to Parafac() and set the dimnames
tv <- toArray(TV[[1]], 30, 16, 15, mode="B")
dimnames(tv)[[1]] <- TV[[4]]
dimnames(tv)[[2]] <- TV[[2]]
dimnames(tv)[[3]] <- TV[[3]]

(tvt3 <- Tucker3(tv, 2, 3, 4))
tvt3$A; tvt3$B; tvt3$C; tvt3$GA
tvt3$cutoff.rd
sort(tvt3$rd)
tvt3$cutoff.sd
sort(tvt3$sd)

(rtvt3 <- Tucker3(tv, 2, 3, 4, robust=TRUE))
rtvt3$A; rtvt3$B; rtvt3$C; rtvt3$GA
rtvt3$cutoff.rd
sort(rtvt3$rd)
rtvt3$cutoff.sd
sort(rtvt3$sd)

## ===============================================================
##
##  Compositional data and robustness

data(ulabor)

(res0 <- Tucker3(ulabor))
res0$A; res0$B; res0$C; res0$GA
res0$cutoff.rd
sort(res0$rd)
res0$cutoff.sd
sort(res0$sd)

(res <- Tucker3(ulabor, robust=TRUE, coda.transform="ilr"))
res$A; res$B; res$C; res$GA
res$cutoff.rd
sort(res$rd)
res$cutoff.sd
sort(res$sd)

(res1 <- Tucker3(ulabor, coda.transform="clr"))
res1$A; res1$B; res1$C; res1$GA
res1$cutoff.rd
sort(res1$rd)
res1$cutoff.sd
sort(res1$sd)
