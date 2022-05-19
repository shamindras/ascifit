
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ascifit

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/ascifit)](https://CRAN.R-project.org/package=ascifit) -->
<!-- badges: end -->

The goal of `ascifit` is to have a set of common utilities to reproduce
results from the paper *Adversarial Sign Corrupted Isotonic Regression*
by Shamindra Shrotriya and Matey Neykov.

The `ascifit` package is developed and maintained by

-   [Shamindra Shrotriya](https://www.shamindras.com/)
-   [Matey Neykov](https://mateyneykov.com/)

## Citation

If you are in `R` you can simply run the following command to get the
`BibTeX` citation for `ascifit`:

``` r
citation("ascifit")
```

Alternatively, please use the following `BibTeX` citation:

``` bibtex
@misc{shrotriya2022asci,
  title  = {{A}dversarial {S}ign-{C}orrupted {I}sotonic {R}egression},
  author = {Shamindra Shrotriya and Matey Neykov},
  year   = {2022},
  eprint = {arXiv:},
  note   = {R package version 0.0.1}
}
```

## Installation

### R package (required)

You can install the development version of `ascifit` from
[GitHub](https://github.com/) using `pak`:

``` r
install.packages("pak") # Run this if you haven't installed pak
pak::pkg_install("shamindras/ascifit")
```

or *alternatively* using `devtools`:

``` r
install.packages("devtools") # Run this if you haven't installed devtools
devtools::install_github("shamindras/ascifit")
```

### Python conda environment (optional)

For your interest, we have provided some `jupyter` `Sagemath` notebooks
used to check some of the symbolic calculations in the
`experiments/python/notebooks` directory in some. To run these
(optional) notebooks locally, we have created the need to first install
the conda environment `cb01`, which contains all of the relevant
packages.

<details>
<summary>
Python Sagemath installation
</summary>

This can be done via the following steps:

1.  Clone the `ascifit` repo and `cd` into the local cloned repo.

``` bash
git clone https://github.com/shamindras/ascifit.git
cd ascifit
```

2.  Install [Anaconda](https://www.anaconda.com/products/individual)
    locally.
3.  At the top level of the `ascifit` directory install the conda
    environment `cb01` as follows:

``` bash
make conda_cb01
```

4.  Activate the conda environment as follows:

``` bash
conda activate cb01
```

5.  You can then run the following command in your terminal, which will
    launch `jupyter lab` in your browser.

``` bash
jupyter lab
```

6.  Using `jupyter lab` in your browser, navigate to
    `experiments/python/notebooks` and start running these as required.

-   Perhaps make a copy of the notebook first to avoid any merge
    conflicts.

</details>

## Reproducibility

After performing the relevant `ascifit` `R` package installation, the
figures from the main paper can be reproduced as follows:

``` r
library(ascifit)
# TODO: Add details
```

## Code of Conduct

Please note that the `ascifit` project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
