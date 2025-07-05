
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `rrcov3way`: Robust Methods for Multiway Data Analysis

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version/rrcov3way)](https://cran.r-project.org/package=rrcov3way)
[![R-CMD-check](https://github.com/valentint/rrcov3way/workflows/R-CMD-check/badge.svg)](https://github.com/valentint/rrcov3way/actions)
[![Codecov test
coverage](https://codecov.io/gh/valentint/rrcov3way/branch/master/graph/badge.svg)](https://app.codecov.io/gh/valentint/rrcov3way?branch=master)
[![downloads](https://cranlogs.r-pkg.org/badges/rrcov3way)](https://cran.r-project.org/package=rrcov3way)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/rrcov3way)](https://cran.r-project.org/package=rrcov3way)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

The package provides robust methods for multiway data analysis by means
of *Parafac* and *Tucker 3* models (Engelen and Hubert (2011)
[doi:10.1016/j.aca.2011.04.043](https://doi.org/10.1016/j.aca.2011.04.043)).
Robust versions for compositional data are also provided (Gallo (2015)
[doi:10.1080/03610926.2013.798664](https://doi.org/10.1080/03610926.2013.798664),
Di Palma et al. (2018)
[doi:10.1080/02664763.2017.1381669](https://doi.org/10.1080/02664763.2017.1381669)).

Several optimization methods alternative to ALS are available (Simonacci
and Gallo (2019)
[doi:10.1016/j.chemolab.2019.103822](https://doi.org/10.1016/j.chemolab.2019.103822),
Simonacci and Gallo (2020)
[doi:10.1007/s00500-019-04320-9](https://doi.org/10.1007/s00500-019-04320-9)).

## Installation

The `rrcov3way` package is on CRAN (The Comprehensive R Archive Network)
and the latest release can be easily installed using the command

    install.packages("rrcov3way")

## Building from source

To install the latest stable development version from GitHub, you can
pull this repository and install it using

    ## install.packages("remotes")
    remotes::install_github("valentint/rrcov3way")

Of course, if you have already installed `remotes`, you can skip the
first line (I have commented it out).

## Example

This is a simple example which shows you basic functions of the package
using the OECD `elind` data set. The data consist of specialization
indices of electronics industries of 23 European countries for the years
1973-1979. The specialization index is defined as the proportion of the
monetary value of an electronic industry compared to the total export
value of manufactured goods of a country compared to the similar
proportion for the world as a whole.

``` r

## Load the package 'rrcov3way' and the data set
##  to be used in the examples:
##
## OECD Electronics Industries Data - Kroonenberg PM (2008)
##  23 countries x 6 industries x 7 years
library(rrcov3way)
#> Robust Methods for Multiway Data Analysis, Applicable also for
#> Compositional Data (version 0.2-5)
#> 
#> Attaching package: 'rrcov3way'
#> The following object is masked from 'package:stats':
#> 
#>     reorder
data(elind)
dim(elind)
#> [1] 23  6  7

##  The frontal slices (mode C) are matrices with dimension 23x6, 
##  representing the data for all countries and all industries 
##  in one year. The labels are:
rownames(elind[,,1])
#>  [1] "CA" "US" "JP" "AS" "NZ" "BL" "DA" "FR" "RF" "GR" "IR" "IT" "PB" "RU" "AU"
#> [16] "FI" "NO" "PO" "SP" "SV" "CH" "TU" "YU"
colnames(elind[,,1])
#> [1] "INFO" "RADI" "TELE" "STRU" "ELET" "COMP"

##  The lateral slices (mode B) represent the data
##  for all countries and years for one industry
##  - matrices with dimension 23 x 7 and the labels are:
rownames(elind[,1,])
#>  [1] "CA" "US" "JP" "AS" "NZ" "BL" "DA" "FR" "RF" "GR" "IR" "IT" "PB" "RU" "AU"
#> [16] "FI" "NO" "PO" "SP" "SV" "CH" "TU" "YU"
colnames(elind[,1,])
#> [1] "78" "79" "80" "82" "83" "84" "85"

##  First of all we center and scale the data, using the default procedures 
##  for centering and scaling.

elind <- do3Scale(elind, center=TRUE, scale=TRUE)

##  Next we perform classical PARAFAC analysis 
##  with the default number of components (ncomp=2).

## default PARAFAC, non-robust, no ilr transformation,
##  extract 3 components
res <- Parafac(elind, ncomp=3)
res
#> Call:
#> Parafac(X = elind, ncomp = 3)
#> 
#> 
#> PARAFAC analysis with  3  components.
#> Fit value: 2.522052 
#> Fit percentage: 57.97 %
#> 
```

## Community guidelines

### Report issues and request features

If you experience any bugs or issues or if you have any suggestions for
additional features, please submit an issue via the
[*Issues*](https://github.com/valentint/rrcov3way/issues) tab of this
repository. Please have a look at existing issues first to see if your
problem or feature request has already been discussed.

### Contribute to the package

If you want to contribute to the package, you can fork this repository
and create a pull request after implementing the desired functionality.

### Ask for help

If you need help using the package, or if you are interested in
collaborations related to this project, please get in touch with the
package maintainer.
