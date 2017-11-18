
<!-- README.md is generated from README.Rmd. Please edit that file -->
pbr
===

Description
-----------

This is the companion R Package to the book *Plant Breeding in R*. It is primarily used to load all of the packages necessary for completing the tasks outlined in *Plant Breeding in R*. However, there are some functions that are provided within the package that may be useful.

Installation
------------

You can install this package using the `devtools` package:

``` r
devtools::install_github("neyhartj/pbr")
```

When you install `pbr`, a number of other packages are installed. The **core** set of packages - those that you might need for most analyses - include:

-   `dplyr`, `tidyr`, `readr`, for data management and modeling
-   `ggplot2`, for visualization
-   `lme4`, for fitting mixed-effect models
-   `agridat`, for various agricultural experiment data

Additionally, other packages are installed that may be useful for specialized analysis, including those for:

-   Managing genetic experiments and QTL mapping

    -   `qtl`, for QTL mapping and simulations
    -   `rrBLUP`, for GWAS and genomic selection

-   Statistics:

    -   `agricolae`, for statistical procedures in agriculture research
    -   `BGLR`, for Bayesian linear regression
    -   `lsmeans`, for calculating least-square means

Usage
-----

Loading the `pbr` package also loads the **core** set of packages (dplyr, tidyr, readr, ggplot, lme4, agridat):

``` r
library(pbr)
#> Loading pbr: dplyr
#> Loading pbr: readr
#> Loading pbr: lme4
#> Loading pbr: agridat
```

Packages that are not in the core set must be loaded individually using the `library()` function, for instance:

``` r
library(qtl)
```

Functions
---------

The `pbr` package includes endogenous functions to support the other packages. Those include

|    Function    |                                  Purpose                                 |
|:--------------:|:------------------------------------------------------------------------:|
|  `dist_env()`  |        Clustering of environments based on phenoypic observations        |
|    `herit()`   | Estimating heritability and variance components from fitted model object |
| `herit_boot()` | Calculating confidence intervals for heritability based on bootstrapping |

Support
-------

Please [open an issue](https://github.com/neyhartj/pbr/issues/new) for support or to comment.
