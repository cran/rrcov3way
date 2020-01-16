## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(rrcov3way))

## VT::14.01.2020 - added ## IGNORE_RDIFF_BEGIN/END
##  when printing the A, B and C matrices - because this could cause
##  differences on some platform, e.g. see for example sign indeterminacy
##

## Example with the UNIDO Manufacturing value added data
data(va3way)
dim(va3way)

## Treat quickly and dirty the zeros in the data set (if any)
va3way[va3way==0] <- 0.001

(res <- Parafac(va3way))
print(res$fit)
## IGNORE_RDIFF_BEGIN
print(res$A)
print(res$B)
print(res$C)
## IGNORE_RDIFF_END
print(res$rd)
print(res$cutoff.rd)

(res.r <- Parafac(va3way, robust=TRUE))
print(res.r$fit)
## IGNORE_RDIFF_BEGIN
print(res.r$A)
print(res.r$B)
print(res.r$C)
## IGNORE_RDIFF_END
print(res.r$rd)
print(res$cutoff.rd)
print(res.r$Hset)
print(res.r$flag)

(res.c <- Parafac(va3way, coda.transform="ilr"))
print(res.c$fit)
## IGNORE_RDIFF_BEGIN
print(res.c$A)
print(res.c$B)
print(res.c$Bclr)
print(res.c$C)
## IGNORE_RDIFF_END
print(res.c$rd)
print(res$cutoff.rd)

(res.rc <- Parafac(va3way, robust=TRUE, coda.transform="ilr"))
print(res.rc$fit)
## IGNORE_RDIFF_BEGIN
print(res.rc$A)
print(res.rc$B)
print(res.rc$Bclr)
print(res.rc$C)
## IGNORE_RDIFF_END
print(res.rc$rd)
print(res$cutoff.rd)
print(res.rc$Hset)
print(res.rc$flag)
