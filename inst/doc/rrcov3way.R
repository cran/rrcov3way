### R code from vignette source 'rrcov3way.Rnw'

###################################################
### code chunk number 1: rrcov3way.Rnw:444-445
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: intro
###################################################
library("rrcov3way")
data("elind")


###################################################
### code chunk number 3: intro-elind1
###################################################
dim(elind)


###################################################
### code chunk number 4: intro-elind2
###################################################
rownames(elind[, , 1])
colnames(elind[, , 1])


###################################################
### code chunk number 5: intro-elind-3
###################################################
rownames(elind[, 1, ])
colnames(elind[, 1, ])


###################################################
### code chunk number 6: intro-elind-preprocess
###################################################
elind <- do3Scale(elind, center = TRUE, scale = TRUE)


###################################################
### code chunk number 7: intro-parafac
###################################################
res <- Parafac(elind, ncomp = 3)
res


###################################################
### code chunk number 8: intro-parafac-1
###################################################
head(res$A)
res$B


###################################################
### code chunk number 9: intro-parafac-rob
###################################################
(resr <- Parafac(elind, ncomp = 3, robust = TRUE))


###################################################
### code chunk number 10: intro-parafac-dd
###################################################
oldpar <- par(mfrow = c(1,2))
plot(res, main = "Classical distance-distance plot")
plot(resr, main = "Robust distance-distance plot")
par(oldpar)


###################################################
### code chunk number 11: intro-parafac-comp
###################################################
oldpar <- par(mfrow = c(1,2))
plot(resr, which = "comp", main = "Paired component plot (mode A)")
plot(resr, which = "comp", mode = "B",
    main = "Paired component plot (mode B)")
par(oldpar)


###################################################
### code chunk number 12: structure-matricization
###################################################
x <- c(1, 0, 1, -1, 2, 1, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 1,
    0, 2, 0, 0, 0, 0, 1)
X <- array(x, dim = c(4, 3, 2))
dimnames(X) <- list(1:4, 1:3, 1:2)
X
(Xa <- unfold(X))
(Xb <- unfold(X, mode = "B"))
(Xc <- unfold(X, mode = "C"))


###################################################
### code chunk number 13: robust-pca
###################################################
data("elind")
(o <- Parafac(elind, robust = TRUE, ncomp.rpca = 11))
rrcov::screeplot(o$pcaobj, main = "Screeplot: elind data")
o1 <- Parafac(elind, robust = TRUE)
cat("\n Selected number of components: ", o1$pcaobj$k, "\n")


###################################################
### code chunk number 14: elind-center
###################################################
elind.cA <- do3Scale(elind, center = TRUE, center.mode = "A")
round(colMeans(elind.cA[, , 1]), 10)


###################################################
### code chunk number 15: elind-center2
###################################################
elind.cAB <- do3Scale(elind.cA, center = TRUE, center.mode = "B")
round(colMeans(elind.cAB[, , 1]), 10)


###################################################
### code chunk number 16: elind-center3
###################################################
elind.cAsB <- do3Scale(elind, center = TRUE, scale = TRUE)


###################################################
### code chunk number 17: elind-center4
###################################################
elind.cAsB <- do3Scale(elind, center = TRUE, scale = mad)


###################################################
### code chunk number 18: center
###################################################
A <- B <- C <- matrix(c(1, 2, 3, 4)) + 10
X0 <- A %*% t(krp(C, B))
X <- toArray(X0, 4, 4, 4)
Parafac(X, ncomp = 1)$fp
Parafac(do3Scale(X, center = TRUE), ncomp = 1)$fp
Parafac(X - mean(X), ncomp = 1)$fp
Parafac(do3Scale(X, center = TRUE, scale = TRUE), ncomp = 1)$fp


###################################################
### code chunk number 19: renormalize
###################################################
t3 <- Tucker3(elind, 3, 2, 2)
t3.norm <- do3Scale(t3, renorm.mode = "A")
rowSums(t3.norm$GA ^ 2)


###################################################
### code chunk number 20: renormalize-2
###################################################
cp <- Parafac(elind, ncomp = 3)
cp.norm <- do3Scale(cp, renorm.mode = "A")
colSums(cp.norm$B ^ 2)
colSums(cp.norm$C ^ 2)


###################################################
### code chunk number 21: rrcov3way.Rnw:867-871
###################################################
data("elind")
t3 <- Tucker3(elind, 3, 2, 2)
xout <- do3Rotate(t3, c(3, 3, 3), rotate = c("A", "B", "C"))
xout$vvalue


###################################################
### code chunk number 22: rrcov3way.Rnw:874-875
###################################################
xout$vvalue


###################################################
### code chunk number 23: rrcov3way.Rnw:878-894
###################################################
w <- c(NA, 0, 0.5, 1, 2.5, 3, 3.5, 4, 5, 10, Inf)
res <- matrix(NA, nrow = length(w), ncol = 7)
for(i in 1:length(w))
{
    res[i, 1] <- res[i, 2] <- res[i, 3] <- w[i]
    if(is.na(w[i]))
        x <- do3Rotate(t3, rotate = c())
    else if(is.finite(w[i]))
        x <- do3Rotate(t3, rep(w[i], 3))
    else
        x <- do3Rotate(t3, rep(1e18, 3))

    res[i, 4:7] <- round(x$vvalue,3)
}
colnames(res) <- c("w(A)", "w(B)", "w(C)", "Core", "A", "B", "C")
rownames(res) <- 1:nrow(res)


###################################################
### code chunk number 24: rrcov3way.Rnw:904-909
###################################################
library("xtable")
print(xtable(res, label="tab:varimax", caption = "Varimax values for
the core array and the component matrices for different joint
varimax rotations to the Tucker3 solution with 3x2x2 components
of the OECD data"), table.placement = "H")


###################################################
### code chunk number 25: rrcov3way.Rnw:916-917
###################################################
x2=do3Rotate(t3, rep(3, 3), rotate = c("A", "B", "C"))


###################################################
### code chunk number 26: rrcov3way.Rnw:919-923
###################################################
print(xtable(cbind(t3$A, x2$x$A), label = "tab:varimax-A",
    caption = "Component matrix A from the unrotated and
    rotated solution with relative weights (3, 3, 3)"),
    table.placement = "H")


###################################################
### code chunk number 27: rrcov3way.Rnw:926-930
###################################################
print(xtable(cbind(t3$B, x2$x$B), label = "tab:varimax-B",
    caption = "Component matrix B from the unrotated and
    rotated solution with relative weights (3, 3, 3)"),
    table.placement = "H")


###################################################
### code chunk number 28: rrcov3way.Rnw:933-937
###################################################
print(xtable(cbind(t3$C, x2$x$C), label = "tab:varimax-C",
    caption = "Component matrix C from the unrotated and
    rotated solution with relative weights (3, 3, 3)"),
    table.placement = "H")


###################################################
### code chunk number 29: rrcov3way.Rnw:940-945
###################################################
print(xtable(rbind(t3$GA, rep(NA, 4), x2$x$GA),
    label = "tab:varimax-core",
    caption = "Core array from the unrotated and rotated
    solution with relative weights (3, 3, 3)"),
    table.placement = "H")


###################################################
### code chunk number 30: oecd-reflect
###################################################
res <- Parafac(elind)
head(res$A)
head(res$B)
res1 <- do3Postprocess(res, reflectA = c(1, -1))
head(res1$A)
head(res1$B)


###################################################
### code chunk number 31: girls
###################################################
data("girls")
dim(girls)
head(girls[, , 1])
sum(girls ^ 2)


###################################################
### code chunk number 32: girls-preprocess
###################################################
X <- do3Scale(girls, center = TRUE, scale = TRUE, only.data = FALSE)
center <- X$center
X <- X$x

average.girl <- as.data.frame(matrix(center, ncol = 8, byrow = TRUE))
dimnames(average.girl) <- list(dimnames(X)[[3]], dimnames(X)[[2]])
average.girl$weight <- average.girl$weight / 10
average.girl$length <- average.girl$length / 10
average.girl$crrump <- average.girl$crrump / 10

sum(X ^ 2)


###################################################
### code chunk number 33: girls-average-plot
###################################################
p <- ncol(average.girl)
plot(rownames(average.girl), average.girl[,1],
    ylim = c(50, 1200),
    type = "n", xlab = "Age", ylab = "")
for(i in 1: p)
{
    lines(rownames(average.girl), average.girl[, i], lty = i, col = i)
    points(rownames(average.girl), average.girl[, i], pch = i, col = i)
}
legend <- colnames(average.girl)
legend[1] <- paste0(legend[1], "*")
legend[2] <- paste0(legend[3], "*")
legend[3] <- paste0(legend[4], "*")
legend("topleft", legend = legend, col = 1:p, lty = 1:p, pch = 1:p)


###################################################
### code chunk number 34: girls-models
###################################################
(t3 <- Tucker3(X, 3, 3, 2))
(cp <- Parafac(X, ncomp = 3))


###################################################
### code chunk number 35: girls-models-robust
###################################################
(t3r <- Tucker3(X, 3, 3, 2, robust = TRUE, ncomp.rpca = 3))
(cpr <- Parafac(X, ncomp = 3, robust = TRUE, ncomp.rpca = 3))


###################################################
### code chunk number 36: girls-ddplot
###################################################
oldpar <- par(mfrow = c(1, 2))
plot(cp, main = "Classical outlier map")
plot(cpr, main = "Robust outlier map")
par(oldpar)


###################################################
### code chunk number 37: girls-paired
###################################################
oldpar <- par(mfrow = c(1, 2))
plot(t3, which = "comp", choices = 1L:2L)
plot(t3, which = "comp", choices = c(1L, 3L))
par(oldpar)


###################################################
### code chunk number 38: girls-paired-B
###################################################
oldpar <- par(mfrow = c(1, 2))
plot(t3, which = "comp", choices = 1L:2L, mode = "B")
plot(t3, which = "comp", choices = c(1L,3L), mode = "B")
par(oldpar)


###################################################
### code chunk number 39: girls-percomp
###################################################
cp <- Parafac(X, ncomp = 2)
cp <- do3Postprocess(cp, reflectA = -1, reflectC = -1)
plot(cp, which = "percomp")


###################################################
### code chunk number 40: girls-allcomp
###################################################
plot(t3, which="allcomp")


###################################################
### code chunk number 41: amino-allcomp
###################################################
data("amino")
amino.cp <- Parafac(amino, ncomp = 3, const = "nonneg")
plot(amino.cp, which = "allcomp",
    xlab = "Wavelength", ylab = "Intensity", mode = "C",
    points = FALSE, legend.position = NULL)


###################################################
### code chunk number 42: girls-jointbiplot
###################################################
t3x <- do3Postprocess(t3, reflectC = -1)
plot(t3x, which = "jbplot")


###################################################
### code chunk number 43: girls-tjplot
###################################################
t3x <- do3Postprocess(t3, reflectB = -1)
plot(t3x, which = "tjplot", arrows = TRUE)


###################################################
### code chunk number 44: girls-tjplot2
###################################################
plot(t3x, which = "tjplot", choices = c(9, 11, 25, 2), arrows = FALSE)


###################################################
### code chunk number 45: Kojima
###################################################
data("Kojima")
dim(Kojima.girls)
head(dimnames(Kojima.girls)[[1]])
dimnames(Kojima.girls)[[2]]
dimnames(Kojima.girls)[[3]]


###################################################
### code chunk number 46: Kojima-center
###################################################
X <- do3Scale(Kojima.girls, center = TRUE, scale = sd)


###################################################
### code chunk number 47: Kojima-Parafac
###################################################
cp <- Parafac(X, ncomp = 3, const = c("orth", "none", "none"), conv = 1e-10)
cp


###################################################
### code chunk number 48: Kojima-Parafac-postprocess
###################################################
cp.norm <- do3Scale(cp, mode = "A")
cp.norm <- do3Postprocess(cp.norm, reflectB = c(-1, -1, -1))
cp.norm <- do3Postprocess(cp.norm, reorder = c(2, 3, 1))
b.pc <- coordinates(cp.norm, mode = "B", type = "principal")


###################################################
### code chunk number 49: Kojima-Parafac-Table1-2
###################################################
round(b.pc, 2)
round(weights(cp.norm), 2)
round(coordinates(cp.norm, mode = "C"), 2)


###################################################
### code chunk number 50: Kojima-percomp
###################################################
plot(cp.norm, which = "percomp")


###################################################
### code chunk number 51: Kojima-girls-congruence
###################################################
cp1 <- Parafac(X, ncomp = 1, const = c("orth", "none", "none"),
    maxit = 10000, conv = 1e-10)
cp2 <- Parafac(X, ncomp = 2, const = c("orth", "none", "none"),
    maxit = 10000, conv = 1e-10)
cp2 <- do3Postprocess(cp2, reorder = TRUE)
cp3 <- Parafac(X, ncomp = 3, const = c("orth", "none", "none"),
    maxit = 10000, conv = 1e-10)
cp3 <- do3Postprocess(cp3, reorder = TRUE)

round(congruence(cp1$A, cp1$A), 2)
round(congruence(cp2$A, cp1$A), 2)
round(congruence(cp3$A, cp1$A), 2)

round(congruence(cp1$A, cp2$A), 2)
round(congruence(cp2$A, cp2$A), 2)
round(congruence(cp3$A, cp2$A), 2)

round(congruence(cp1$A, cp3$A), 2)
round(congruence(cp2$A, cp3$A), 2)
round(congruence(cp3$A, cp3$A), 2)


###################################################
### code chunk number 52: rrcov3way.Rnw:1542-1543
###################################################
options(width = 120)


###################################################
### code chunk number 53: Arno-data
###################################################
data("Arno")
Arno[, , 1]


###################################################
### code chunk number 54: rrcov3way.Rnw:1549-1550
###################################################
options(width = 60)


###################################################
### code chunk number 55: Arno-t3-comp
###################################################
(Arnot3 <- Tucker3(Arno, P = 2, Q = 2, R = 1,
    center = TRUE, center.mode = "AB", coda.transform = "ilr"))


###################################################
### code chunk number 56: Arno-t3-comp-modeA
###################################################
plot(Arnot3, which = "comp", main = "(a) Paired component plot (mode A)")


###################################################
### code chunk number 57: Arno-t3-comp-modeB
###################################################
plot(Arnot3, which = "comp", mode = "B",
    main = "(b) Paired component plot (mode B)")


###################################################
### code chunk number 58: Arno-t3-comp-joint
###################################################
plot(Arnot3, which = "jbplot", main = "Joint biplot")


###################################################
### code chunk number 59: Arno-t3-comp-trajectory
###################################################
plot(Arnot3, which = "tjplot", main = "Trajectory biplot")


###################################################
### code chunk number 60: amino-surface
###################################################
data("amino")
x <- as.numeric(dimnames(amino)[[2]])
y <- as.numeric(dimnames(amino)[[3]])
persp(x, y, amino[2, , ], ticktype = "detailed",
    theta = -15, phi = 30, col = "lightblue",
    xlab = "Emission wavelength", ylab = "Excitation wavelength",
    zlab = "Intensity")


###################################################
### code chunk number 61: amino-excitation
###################################################
library("ggplot2")
library("reshape2")
mamino1= melt(amino[, 1 ,], value.name = "Intensity")
mamino2= melt(amino[, 30 ,], value.name = "Intensity")
mamino <- rbind(cbind(id = dimnames(amino)[[2]][1], mamino1),
    cbind(id = dimnames(amino)[[2]][30], mamino2))
colnames(mamino)[2] <- "Sample"
colnames(mamino)[3] <- "Emission Wavelength"
p <- ggplot(data=mamino, aes(x = `Emission Wavelength`, y = Intensity)) +
    geom_line(aes(colour = Sample)) +
    facet_wrap(~id, nrow = 2, scales = "free") +
theme_bw()
p


###################################################
### code chunk number 62: amino-PARAFAC
###################################################
set.seed(1234)
(p1 <- Parafac(amino, 3, const = "nonneg", start = "random"))
(p2 <- Parafac(amino, 3, const = "none", start = "random"))


###################################################
### code chunk number 63: amino-PARAFAC-nonneg
###################################################
plot(p1, which = "allcomp", mode = "B", points = FALSE,
    legend.position = NULL, xlab = "Wavelength",
    main = "(b)")


###################################################
### code chunk number 64: amino-PARAFAC-none
###################################################
plot(p2, which = "allcomp", mode = "B", points = FALSE,
    legend.position = NULL, xlab = "Wavelength",
    main = "(a)")


