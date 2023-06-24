######
##  VT::10.01.2020
##
##
##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/R/rrcov3way", load_code=roxygen2:::load_installed)
##
#' Amino acids fluorescence data.
#'
#' A data set containing five simple laboratory-made samples where each sample contains
#'  different amounts of tyrosine, tryptophan and phenylalanine dissolved in phosphate
#'  buffered water. The samples were measured by fluorescence (excitation 240-300 nm,
#'  emission 250-450 nm, 1 nm intervals) on a PE LS50B spectrofluorometer.
#'
#' @name amino
#' @docType data
#' @usage data(amino)
#' @format A three-way array with dimension \code{5x201x61}.
#'  The first dimension refers to the 5 samples. The second dimension
#'  refers to the emission measurements (250-450nm, 1nm intervals).
#'  The third dimension refers to the excitation (240-300 nm, 1nm intervals).
#' @source
#'  \url{https://ucphchemometrics.com/datasets/}.
#' @references
#'  Bro, R, PARAFAC: Tutorial and applications, Chemometrics and Intelligent Laboratory Systems, 1997, 38, 149-171
#'  Bro, R, Multi-way Analysis in the Food Industry. Models, Algorithms, and Applications. 1998. Ph.D. Thesis,
#'      University of Amsterdam (NL) & Royal Veterinary and Agricultural University (DK).
#'  Kiers, H.A.L. (1998) A three-step algorithm for Candecomp/Parafac analysis of large data
#'      sets with multicollinearity, Journal of Chemometrics, 12, 155-171.
#' @examples
#'  \dontrun{
#'
#'  data(amino)
#'  ##  Plotting Emission spectra
#'  oldpar <- par(mfrow=c(2,1))
#'  matplot(t(amino[,,1]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence emission spectra")
#'  matplot(t(amino[,,5]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence emission spectra")
#'  par(oldpar)
#'
#'  ##  Plotting excitation spectra
#'  oldpar <- par(mfrow=c(2,1))
#'  matplot(t(amino[,1,]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence excitation spectra")
#'  matplot(t(amino[,30,]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence excitation spectra")
#'  par(oldpar)
#'
#' }
#' @keywords datasets
NULL
#' Dorrit fluorescence data.
#'
#' A data set with 27 synthetic samples containing different concentrations of four analytes
#'  (hydroquinone, tryptophan, phenylalanine and dopa) measured in a Perkin-Elmer
#'  LS50 B fluorescence spectrometer.
#'
#'  Each fluorescence landscape corresponding to each
#'  sample in the original data set consists of 233 emission wavelengths
#'  (250-482 nm) and 24 excitation wavelengths
#'  (200-315 nm taken each 5 nm). The fluorescence data is three-way.
#'  Ideally, the data is trilinear, the components
#'  of the modes corresponding to concentrations (27), emission
#'  spectra (233) and excitation spectra (24). Hence a three-way
#'  PARAFAC model should be capable of uniquely and meaningfully
#'  describing the variation in the data set.
#'
#'  The data set is modified as described in Engelen and Hubert (2011):
#'  the emission wavelengths are taken at 2 nm,
#'  noisy parts situated at the excitation wavelengths from 200 to 230 nm and
#'  at emission wavelengths below 250 nm are excluded. The severe Rayleigh
#'  scattering areas present in all samples are replaced by interpolated values.
#'  Thus we end up with a \code{(27 x 116 x 18)} data array.
#'
#' @name dorrit
#' @docType data
#' @usage data(dorrit)
#' @format A three-way array with dimension \code{27 x 116 x 18}.
#'  The first dimension refers to the 27 samples. The second dimension
#'  refers to the emission measurements (251-481nm, 2nm intervals).
#'  The third dimension refers to the excitation (230-315nm, 5nm intervals).
#' @source
#'  \url{https://ucphchemometrics.com/datasets/}.
#' @references
#'  Bro, R, Sidiropoulos, ND and Smilde, AK (2002). Maximum likelihood
#'  fitting using ordinary least squares algorithms.
#'  \emph{Journal of Chemometrics}, \bold{16}(8--10), 387--400.
#'
#'  Riu, J and Bro, R (2003) Jack-knife estimation of standard errors and outlier detection in PARAFAC models.
#'      \emph{Chemometrics and Intelligent Laboratory Systems}, 65(1), 35--49
#'
#'  Engelen, S and Hubert, M (2011) Detecting outlying samples in a parallel factor analysis model,
#'      \emph{Analytica Chemica Acta} 705 155--165.
#'
#'  Baunsgaard, D (1999) Factors Affecting 3-way Modelling (PARAFAC) of Fluorescence
#'  Landscapes, Royal Veterinary and Agricultural University,
#'  Department of Dairy and Food Science, Frederiksberg, Denmark.
#' @examples
#'  \dontrun{
#'
#'  data(dorrit)
#'  ##  Plotting Emission spectra
#'  oldpar <- par(mfrow=c(2,1))
#'  matplot(t(dorrit[,,1]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence emission spectra")
#'  matplot(t(dorrit[,,5]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence emission spectra")
#'  par(oldpar)
#'
#'  ##  Plotting excitation spectra
#'  oldpar <- par(mfrow=c(2,1))
#'  matplot(t(dorrit[,1,]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence excitation spectra")
#'  matplot(t(dorrit[,30,]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence excitation spectra")
#'  par(oldpar)
#'
#'  persp(as.numeric(dimnames(dorrit)[[2]]),
#'      as.numeric(dimnames(dorrit)[[3]]), dorrit[4,,],
#'      xlab="Emission", ylab="Excitation", zlab="Intensity",
#'      theta = 30, phi = 30, expand = 0.5, col = "lightblue",
#'      ticktype="detailed")
#'
#'  pp <- Parafac(dorrit, ncomp=4, robust=TRUE)
#'  plot(pp)
#'
#' }
#' @keywords datasets
NULL
