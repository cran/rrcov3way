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
#' Fictious data set of Kiers and Mechelen (2001).
#'
#' A fictitious data presented by Kiers and Mechelen (2001) to illustrate Tucker3 
#'  model: A set of six persons with scores on five response variables for four 
#'  different situations. The response variables indicate to what extent each 
#'  individual displays an emotional, sensitive, caring, thorough, or accurate 
#'  behavior. The data set is represented as a 6 x 5 x 4 array. The data are 
#'  chosen such that they correspond to a Tucker3 model with \code{P=2, Q=2, R=2}.
#'
#' @name Kiers2001
#' @docType data
#' @usage data(Kiers2001)
#' @format A three-way array with dimension \code{6 x 5 x 4}.
#'  The first dimension refers to the 6 persons. The second dimension
#'  refers to the five responce variables.
#'  The third dimension refers to four different situations.
#' @source
#' Kiers HAL, Mechelen IV (2001) Three-way component analysis: Principles and
#' illustrative application. Psychological Methods 6(1):84--110 
#' @references
#' Todorov, V., Simonacci, V., Gallo, M., and Di Palma, M. (2025). Robust tools
#'  for three-way component analysis of compositional data: The R package
#'  rrcov3way. Behaviormetrika. In press.
#'
#' @examples
#'
#'  data(Kiers2001)
#'
#'  t3 <- Tucker3(Kiers2001, P=2, Q=2, R=2)
#'  t3
#'  t3$A
#'  t3$B
#'  t3$C
#'  t3$G
#'  plot(t3)
#'  plot(t3, which="jbplot", xlim=c(0, 2))
#'
#' @keywords datasets
NULL
#' Student satisfaction data
#'
#' The primary assessment tool for analyzing student satisfaction is an 
#'  annual questionnaire organized into five categories: degree program, course 
#'  characteristics, teaching, equipment, and overall satisfaction. 
#'  Each questionnaire comprises 18 standardized questions rated on a ten-point 
#'  Likert scale (1 = "Strongly Disagree" to 10 = "Strongly Agree"). The  
#'  original microdata were aggregated by 10 faculties and six academic years. 
#'  Thus, we have a threeway array with dimensions \code{10 x 18 x 6}.
#'
#' @name satisfaction
#' @docType data
#' @usage data(satisfaction)
#' @format A three-way array with dimension \code{10 x 18 x 6}.
#'  The first dimension refers to the 10 faculties. The second dimension
#'  refers to 18 standardized questions rated on a ten-point Likert scale 
#'  (1 = "Strongly Disagree" to 10 = "Strongly Agree").
#'  The third dimension refers to six consequtive academic years (2012--2017).
#' @source
#'  Italian universities are mandated by Law No. 370/99 to evaluate
#'  teaching quality through student opinion surveys. The National Agency for
#'  the Evaluation of Universities and Research (ANVUR), established in 2006,
#'  supervises the periodic assessment of academic quality. Following Presidential
#'  Decree 76/2010, ANVUR standardized methodologies for evaluating institutions
#'  and degree programs, with a strong focus on student involvement.
#'  The University of Florence provided microdata, which were later aggregated 
#'  into data for 10 faculties, 18 questions and 6 academic years.
#' 
#' Simonacci V, Gallo M (2017) Statistical tools for student evaluation of academic
#'  educational quality. Quality & Quantity 51(2):565--579 
#' @references
#'  Gallo M, Simonacci V, Todorov V (2021) A compositional three-way approach
#'  for student satisfaction analysis. In: Filzmoser P, Hron K, Martin--Fernandez
#'  JA, Palarea-Albaladejo J (eds) Advances in Compositional Data Analysis,
#'  Springer, Cham, pp 143--162
#'
#' Todorov, V., Simonacci, V., Gallo, M., and Di Palma, M. (2025). Robust tools
#'  for three-way component analysis of compositional data: The R package
#'  rrcov3way. Behaviormetrika. In press.
 
#' @examples
#'
#'  data(satisfaction)
#'
#'  t3 <- Tucker3(satisfaction, P=2, Q=2, R=1, coda.transform="clr")
#'  t3
#'  t3$A
#'  t3$B
#'  t3$C
#'  t3$G
#'  plot(t3)
#'  plot(t3, which="jbplot", xlim=c(-1, 1))
#'
#' @keywords datasets
NULL
