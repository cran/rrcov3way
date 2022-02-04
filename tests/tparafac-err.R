## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(rrcov3way))

set.seed(123456)

## Example with the UNIDO Manufacturing value added data
data(va3way)
dim(va3way)

## Treat quickly and dirty the zeros in the data set (if any)
va3way[va3way==0] <- 0.001

## IGNORE_RDIFF_BEGIN
res <- Parafac(va3way, trace=TRUE)              # tracing
## IGNORE_RDIFF_END

## Using robustness with clr transformation
try(res <- Parafac(va3way, robust=TRUE, coda.transform="clr"))

## Rejected values of parameter 'crit'
try(res <- Parafac(va3way, crit=c(1:10)))       # length different than 1
try(res <- Parafac(va3way, crit=-1))            # crit non-positive
try(res <- Parafac(va3way, crit=2))             # crit >= 1

res <- Parafac(va3way, crit=0.2)                # crit < 0.5 --> crit=1-crit

## Test cp_als(): the input array
try(rrcov3way:::cp_als(va3way))                 # missing ncomp

set.seed(98765)
rrcov3way:::cp_als(va3way, ncomp=2)             # OK, 3-way array
rrcov3way:::cp_als(unfold(va3way), ncomp=2,
    n=49, m=5, p=14)                            # OK, unfolded 3-way array

try(rrcov3way:::cp_als("abc", ncomp=2))         # error, not an array or matrix

try(rrcov3way:::cp_als(unfold(va3way), ncomp=2))# missing dimensions
try(rrcov3way:::cp_als(unfold(va3way), ncomp=2,
    n=50, m=5, p=14))                           # n != dim(Xa)[1]
try(rrcov3way:::cp_als(unfold(va3way), ncomp=2,
    n=49, m=1, p=14))                           # m*p != dim(Xa)[2]

## Test cp_als(): the constraints
try(Parafac(va3way, const="abc"))               # wrong constraint
res <- Parafac(va3way, const=c("none", "none")) # length of const < 3
res$const

## Test cp_als(): the initial values
try(Parafac(va3way, start=c(1:2)))      # wrong start
try(Parafac(va3way, start="abc"))       # wrong start

Parafac(va3way, start="svd")
Parafac(va3way, const="nonneg", start="svd")
Parafac(va3way, const="orth", start="svd")
Parafac(va3way, const="zerocor", start="svd")

set.seed(12345)
n <- 49
m <- 5
p <- 14
r <- 2

A <- matrix(runif(max(n,r) * r), max(n,r))[1:n, , drop=FALSE]
B <- matrix(runif(max(m,r) * r), max(m,r))[1:m, , drop=FALSE]
C <- matrix(runif(max(p,r) * r), max(p,r))[1:p, , drop=FALSE]
Parafac(va3way, const="nonneg", start=list(A=A, B=B, C=C))

A <- pracma::orth(matrix(rnorm(max(n,r) * r), max(n,r)))[1:n, , drop=FALSE]
B <- pracma::orth(matrix(rnorm(max(m,r) * r), max(m,r)))[1:m, , drop=FALSE]
C <- pracma::orth(matrix(rnorm(max(p,r) * r), max(p,r)))[1:p, , drop=FALSE]
Parafac(va3way, start=list(A=A, B=B, C=C))
try(Parafac(va3way, const="nonneg", start=list(A=A, B=B, C=C)))
