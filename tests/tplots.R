## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(rrcov3way))
set.seed(123456)

#############
##
## Example with the UNIDO Manufacturing value added data

data(va3way)
dim(va3way)

## Treat quickly and dirty the zeros in the data set (if any)
va3way[va3way==0] <- 0.001

##  Tucker 3 =====================================================
##

## IGNORE_RDIFF_BEGIN
(res <- Tucker3(va3way))
print(res$fit)
print(res$A)

res1 <- reflect(res, mode="C", rsign=-1)
print(res1$C)
## IGNORE_RDIFF_END

## Distance-distance plot
plot(res, which="dd", main="Distance-distance plot")

## Paired component plot, mode A
plot(res, which="comp", main="Paired component plot (mode A)")

## choices must be of length 2 (warning)
plot(res, which="comp", choices=1:3, main="Paired component plot (mode A)")

## Paired component plot, mode B
plot(res, which="comp", mode="B", main="Paired component plot (mode B)")

## Paired component plot, mode C
plot(res, which="comp", mode="C", main="Paired component plot (mode C)")

## All component plot - mode C (default)
plot(res, which="allcomp", main="All component plot")

## All component plot - mode A
plot(res, which="allcomp", mode="A", main="All component plot, mode A")

## All component plot - mode B
plot(res, which="allcomp", mode="B", main="All component plot, mode B")

## Joint biplot
plot(res, which="jbplot", main="Joint biplot")
plot(res1, which="jbplot", main="Joint biplot")             # reflected C

## Trajectory
plot(res, which="tjplot", main="Trajectory biplot")
plot(res, which="tjplot", choices=c(1:4), arrows=FALSE, main="Trajectory biplot")
plot(res1, which="tjplot", main="Trajectory biplot")        # reflected C

## Parafac =======================================================
##

## IGNORE_RDIFF_BEGIN
(res <- Parafac(va3way))
print(res$fit)
print(res$A)
## IGNORE_RDIFF_END

## Distance-distance plot
plot(res, which="dd", main="Distance-distance plot")

## Plot Orthonormalized A-mode component plot
plot(res, which="comp", mode="A", main="Component plot, A-mode")

## choices must be of length 2 (warning)
plot(res, which="comp", choices=1:3, main="Component plot, A-mode")

## Plot Orthonormalized B-mode component plot
plot(res, which="comp", mode="B", main="Component plot, B-mode")

## Plot Orthonormalized C-mode component plot
plot(res, which="comp", mode="C", main="Component plot, C-mode")

## Per component plot
plot(res, which="percomp", main="Per component plot")

## All component plot
plot(res, which="allcomp", main="All component plot")

## .ddplot =======================================================
data(elind)
cp <- Parafac(elind)
rcp <- Parafac(elind, robust=TRUE)
rrcov3way:::.ddplot(cp)                     # non-robust
rrcov3way:::.ddplot(rcp)                    # robust

rrcov3way:::.ddplot(rcp, labs=NULL)         # the labels will be 1:n
rrcov3way:::.ddplot(cp, id.n=5)             # id.n is specified
try(rrcov3way:::.ddplot(cp, id.n=-1))       # id.n must be between 1 and n
