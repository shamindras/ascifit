
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ascifit

The goal of `ascifit` is to have a set of common utilities to reproduce
results from the paper *Adversarial Sign Corrupted Isotonic Regression*.

This repository is developed and maintained by

- [Shamindra Shrotriya](https://www.shamindras.com/)
- [Matey Neykov](https://www.shamindras.com/)

## Citation

To reference this work, please cite the latest version of the paper

``` bib
@misc{shrotriya2022ascireg,
  title   = {Adversarial Sign-Corrupted Isotonic Regression},
  author  = {Shamindra Shrotriya and Matey Neykov},
  year    = 2022,
  eprint  = {arXiv:2207.07075}
}
```

## Installation

In order to reproduce the figures you need to first install the `renv`
package from CRAN in a separate `R` session.

``` r
install.packages("renv")
```

You can then clone the repository locally as usual, for example using
SSH

``` bash
git@github.com:shamindras/ascifit.git #SSH approach
```

## Reproducibility

1.  After locally cloning the repository, you can open the
    `ascifit.Rproj` project file in RStudio. This will create a separate
    RStudio session.

2.  Once `ascifit.Rproj` is opened in RStudio, you can open the file

    ``` bash
    experiments/R/paper_results/create-ascifit-plots.R
    ```

3.  Select all the contents, and run the entire file. It will produce
    both `png/TeX` versions of the three figures, as follows:

    ``` bash
    out_ascifit1_50_reps.png # -> Figure 1 in paper
    sim_plt_50_reps.png # -> Figure 2 in paper
    ```

    **Note:** In the paper we use the `.tex` versions for the figures
    not the `.png` to reduce the resulting paper file size. However both
    the `.tex` and corresponding named `.png` files produce equivalent
    plots. The `.png` formats are typically easier to view.
