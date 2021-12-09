
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ascifit

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/ascifit)](https://CRAN.R-project.org/package=ascifit) -->
<!-- badges: end -->

The goal of `ascifit` is to have a set of common utilities to reproduce
results from the paper *Isotonic Regression Under Adversarial Sign
Corruption* by Matey Neykov, Shamindra Shrotriya.

The `ascifit` package is developed and maintained by

-   [Matey Neykov](https://mateyneykov.com/)
-   [Shamindra Shrotriya](https://www.shamindras.com/)

## Installation

You can install the development version of `sce` from
[GitHub](https://github.com/) using `pak`:

``` r
# install.packages("pak") # Run this if you haven't installed pak
pak::pkg_install("shamindras/ascifit")
```

or alernatively using `devtools`:

``` r
# install.packages("devtools") # Run this if you haven't installed devtools
devtools::install_github("shamindras/ascifit")
```

## Reproducibility

To reproduce figures from the main paper you can run the following:

``` r
library(ascifit)
```
