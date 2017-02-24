# pbr

## Description

This is the companion R Package to the book *Plant Breeding in R*. It is primarily used to load all of the packages necessary for completing the tasks outlined in *Plant Breeding in R*. However, there are some functions that are provided within the package that may be useful.

## Installation

You can install this package using the `devtools` package:

```
devtools::install_github("neyhartj/pbr")
```

## Purpose

The `pbr` package serves two purposes: 1) install and load necessary packages for analysis, and 2) provide additional functions.

### Installed Packages

When you install `pbr`, the following packages are also installed:

  Package         Purpose
-------------- --------------------
  dplyr           Data management
  tidyr           Data management
  readr           Reading data
  ggplot2         Plotting and graphics
  qtl             QTL mapping and simulations
  rrBLUP          GWAS and genomic selection
  agricolae       Statistical analysis
  lattice         Lattice graphics/data
  
Loading the `pbr` package also loads those packages:

```
library(pbr)
```

### Functions

The `pbr` package includes endogenous functions to support the other packages. Those include

  Function            Purpose
---------------- ------------------------------
  plot_AMMI()       Flexible plotting of AMMI models
  dist_env()        Clustering of environments based on phenoypic observations

## Support

Please [open an issue](https://github.com/neyhartj/pbr/issues/new) for support or to comment.
