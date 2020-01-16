######
##  VT::10.01.2020
##
##
##  roxygen2::roxygenise("C:/projects/statproj/R/rrcov3way")
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
#' @format A three-way array with dimension 5x201x61.
#'  The first dimension refers to the 5 samples. The second dimension refers to the emission measurements (250-450nm, 1nm intervals).
#'  The third dimension refers to the excitation (240-300 nm, 1nm intervals).
#'
#'  @source
#'  \url{http://www.models.life.ku.dk/Amino_Acid_fluo}.
#'
#'  @references
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
#'  par <- oldpar
#'
#'  ##  Plotting excitation spectra
#'  oldpar <- par(mfrow=c(2,1))
#'  matplot(t(amino[,1,]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence excitation spectra")
#'  matplot(t(amino[,30,]), type="l",
#'      xlab="Wavelength/nm", ylab="Intensity",
#'      main="Fluorescence excitation spectra")
#'  par <- oldpar
#'
#' }
#' @keywords datasets
NULL
