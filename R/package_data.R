#' Genotypic data on a two-row barley population
#'
#' @description
#' A \code{matrix} of genomewide marker data. 1622 barley individuals were genotyped
#' with 2716 DNA markers. The data has been curated and filtered to remove markers
#' and individuals with a high amount of missing data.
#'
#' @format An object of class \code{matrix} with 1622 rows and 2716 columns. The row names
#' are the names of barley lines, and the column names are the names of markers.
#'
#' @source
#' This data was downloaded from the Triticeae Toolbox at https://triticeaetoolbox.org/barley/.
#'
#'
"tr_cap_genos_mat"

#' Genotypic data on a two-row barley population
#'
#' @description
#' A \code{data.frame} of genomewide marker data. 1622 barley individuals were genotyped
#' with 2716 DNA markers. The data has been curated and filtered to remove markers
#' and individuals with a high amount of missing data.
#'
#' @format An object of class \code{data.frame} with 2716 rows and 1626 columns. The first 4
#' columns contain the marker names, the allele states, the chromosome, and the position
#' on each chromosome. The remaining columns contain the genotype states for each barley line.
#'
#' @source
#' This data was downloaded from the Triticeae Toolbox at https://triticeaetoolbox.org/barley/.
#'
"tr_cap_genos_hmp"


#' Phenotypic data on a two-row barley population
#'
#' @description
#' A \code{data.frame} of phenotypic data. 1598 barley individuals were phenotyped in 6 trials
#' for grain yield and plant height.
#'
#' @format An object of class \code{data.frame} with 2120 rows and 4 columns. Columns include
#' the names of barley lines, the trial in which the phenotype was recorded, values for grain yield,
#' and values for plant height.
#'
#' @source
#' This data was downloaded from the Triticeae Toolbox at https://triticeaetoolbox.org/barley/.
#'
"tr_cap_phenos"

