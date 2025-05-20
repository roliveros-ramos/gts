
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gts <a href="https://roliveros-ramos.github.io/gts/"><img src="man/figures/logo_small.png" align="right" height="124" /></a>

<!-- badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/gts)](https://CRAN.R-project.org/package=gts)
![GitHub R package
version](https://img.shields.io/github/r-package/v/roliveros-ramos/gts?label=GitHub)
[![R-CMD-check](https://github.com/roliveros-ramos/gts/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/roliveros-ramos/gts/actions/workflows/R-CMD-check.yaml)
[![GitHub
issues](https://img.shields.io/github/issues/roliveros-ramos/gts)](https://github.com/roliveros-ramos/gts/issues)
[![OpenSSF Best
Practices](https://www.bestpractices.dev/projects/2132/badge)](https://www.bestpractices.dev/projects/2132)
[![](http://cranlogs.r-pkg.org/badges/gts)](https://CRAN.R-project.org/package=gts)
[![](http://cranlogs.r-pkg.org/badges/grand-total/gts)](https://CRAN.R-project.org/package=gts)
[![codecov](https://codecov.io/gh/roliveros-ramos/gts/graph/badge.svg?token=HELOL3WS4G)](https://app.codecov.io/gh/roliveros-ramos/gts)
<!-- badges: end -->

### Gridded Time-Series manipulation and analysis

Management of gridded time-series (arrays with a time dimension),
normally geographical (longitude, latitude, depth/altitude). The package
provides methods to read and write gridded time-series (GTS) objects and
manage them using R time-series methods.

See <https://roliveros-ramos.github.io/gts/> for more details.

### Installation

``` r
# The easiest way to get gts is to install it from CRAN:
 #install.packages("gts")

# Alternatively, install the stable development version from OSMOSE drat repository:
install.packages("gts", repo="https://osmose-model.github.io/drat/")

# Or the development version from GitHub:
# install.packages("remotes")
remotes::install_github("roliveros-ramos/gts")
```

### Usage

For a quick introduction, check the worked the examples available from
the package:

``` r
library(gts)
vignette("gts")
```

For a more detailed explanation of the package philosophy, you can read
the article [gts: an R package for fitting complex ecological
models](https://doi.org/10.1111/2041-210X.14452).

### Contributions

If you find any bug, have questions about the documentation or requests
for enhancements, please [open an
issue](https://github.com/roliveros-ramos/gts/issues).

Contributions are accepted as pull requests. Please note that the gts
package is released with a [Contributor Code of
Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
By contributing to this project, you agree to abide by its terms.
