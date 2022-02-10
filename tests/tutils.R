## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(rrcov3way))
set.seed(123456)

## unfold() ======================================================
data(elind)
try(unfold(elind[,,1]))
(a <- unfold(elind, mode="A"))
(b <- unfold(elind, mode="B"))
(c <- unfold(elind, mode="C"))

## tallArray(), wideArray() ======================================
all.equal(elind[,,1], rrcov3way:::tallArray(elind)[1:23,], check.attributes=FALSE)
all.equal(elind[,,7], rrcov3way:::tallArray(elind)[(6*23+1):(7*23),], check.attributes=FALSE)
all.equal(a, rrcov3way:::wideArray(elind))

## ilr, clr, gm ==================================================
data(ulabor)
rrcov3way:::.ilrV(ulabor[,,1])
is.data.frame(rrcov3way:::.ilrV(as.data.frame(ulabor[,,1])))

set.seed(1234)
x <- rnorm(100)
x[100] <- 0
rrcov3way:::.gm(x)
try(rrcov3way:::.gm(c(1, 5, 6, "a")))

## orthogonal matrices ===========================================
v <- matrix(rnorm(3), 3)
v <- v / norm(v, 'F')
x <- diag(3) - 2 * (v %*% t(v))
rrcov3way:::is.orthogonal(x)
rrcov3way:::is.orthonormal(x)
zapsmall(x %*% t(x))

## reorder, do3Postprocess, do3Scale =============================
data(elind)

cp <- Parafac(elind, const="orth")
cp_r1 <- reorder(cp)                    # order =TRUE (default)
cp_r2 <- reorder(cp, order=FALSE)       # order = FALSE
cp_r3 <- reorder(cp, order=c(2,1))      # order is a vector
try(cp_r4 <- reorder(cp, order=c(1,1))) # order is wrong
cpx <- do3Postprocess(cp, -1, -1, -1, reorder=TRUE)
cpy_A <- do3Scale(cp, renorm.mode="A")
cpy_B <- do3Scale(cp, renorm.mode="B")
cpy_C <- do3Scale(cp, renorm.mode="C")

t3 <- Tucker3(elind)
t3_r1 <- reorder(t3)                    # order =TRUE (default), mode=A
t3_r2 <- reorder(t3, order=FALSE)       # order = FALSE
t3_r3 <- reorder(t3, order=c(2,1))      # order is a vector
try(t3_r4 <- reorder(t3, order=c(1,1))) # order is wrong
t3_r5 <- reorder(t3, mode="B")          # order =TRUE (default), mode=B
t3_r5 <- reorder(t3, mode="C")          # order =TRUE (default), mode=C
t3x <- do3Postprocess(t3, -1, -1, -1, TRUE, TRUE, TRUE)
t3y_A <- do3Scale(t3, renorm.mode="A")
t3y_B <- do3Scale(t3, renorm.mode="B")
t3y_C <- do3Scale(t3, renorm.mode="C")

##  coordinates(), weights() - for Parafac and Tucker 3 ==========
data(elind)

cp <- Parafac(elind)
try(coordinates(cp))
try(weights(cp))

cp_a <- Parafac(elind, const=c("orth", "none", "none"))
## IGNORE_RDIFF_BEGIN
coordinates(cp_a, mode="A")
coordinates(cp_a, mode="A", type="unit")
coordinates(cp_a, mode="A", type="principal")
weights(cp_a, mode="A")
## IGNORE_RDIFF_END

cp_b <- Parafac(elind, const=c("none", "orth", "none"))
## IGNORE_RDIFF_BEGIN
coordinates(cp_b, mode="B")
coordinates(cp_b, mode="B", type="unit")
coordinates(cp_b, mode="B", type="principal")
weights(cp_b, mode="B")
## IGNORE_RDIFF_END

cp_c <- Parafac(elind, const=c("none", "none", "orth"))
## IGNORE_RDIFF_BEGIN
coordinates(cp_c, mode="C")
coordinates(cp_c, mode="C", type="unit")
coordinates(cp_c, mode="C", type="principal")
weights(cp_c, mode="C")
## IGNORE_RDIFF_END

t3 <- Tucker3(elind)
coordinates(t3, mode="A")
coordinates(t3, mode="A", type="unit")
coordinates(t3, mode="A", type="principal")
weights(t3, mode="A")

coordinates(t3, mode="B")
coordinates(t3, mode="B", type="unit")
coordinates(t3, mode="B", type="principal")
weights(t3, mode="B")

coordinates(t3, mode="C")
coordinates(t3, mode="C", type="unit")
coordinates(t3, mode="C", type="principal")
weights(t3, mode="C")

## mtrace() ======================================================
(a <- matrix(c(5,2,3, 4,-3,7, 4,1,2), ncol=3))
(b <- matrix(c(1,0,1, 0,1,2, 1,0,3), ncol=3))

mtrace(a)
mtrace(b)

## tr(A+B)=tr(A)+tr(B)
all.equal(mtrace(a) + mtrace(b), mtrace(a+b))

## tr(A)=tr(A')
all.equal(mtrace(a), mtrace(t(a)))

## tr(alphA)=alphatr(A)
alpha <- 0.5
all.equal(mtrace(alpha*a), alpha*mtrace(a))

##  tr(AB)=tr(BA)
all.equal(mtrace(a %*% b), mtrace(b %*% a))

##  tr(A)=tr(BAB-1)
all.equal(mtrace(a), mtrace(b %*% a %*% solve(b)))

try(mtrace(c(1, 2, 3, 4)))                  # not a matrix
try(mtrace(matrix(1:20, nrow=5)))           # not a square matrix
try(mtrace(matrix(rep("A", 25), nrow=5)))   # not a numeric matrix


## krp() =========================================================
a <- matrix(1:12, 3, 4)
b <- diag(1:4)
c <- diag(1:3)
krp(a, b)
krp(b, a)
try(krp(a, c))

## congruence ====================================================
## Ignore warning messages of rrcov built version
## IGNORE_RDIFF_BEGIN
require(rrcov)
## IGNORE_RDIFF_END

data(delivery, package="robustbase")
X <- getLoadings(PcaClassic(delivery))
Y <- getLoadings(PcaHubert(delivery, k=3))
round(congruence(X,Y),3)
try(congruence(X,Y[1:2,]))
