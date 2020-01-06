#' Example data from a two-row barley population
#'
#' @description
#' A \code{list} of phenotypic and molecular marker data on a population of 1590
#' two-row barley lines.
#'
#' @format An object of class \code{list} with three elements: 1) \code{pheno} - a
#' \code{data.frame} with 2104 rows and 4 columns. The first column contains the
#' line name, the second column contains the trial name, the and the remaining two
#' columns contain phenotypic data on two quantitative traits; 2) \code{geno} - a
#' \code{data.frame} with 2187 rows and 1594 columns. The first 4 columns contain
#' information on the name, allels, chromosome, and genetic
#' map position (respectively) of 2187 molecular markers. Remaining columns
#' contain information on the genotypic states at each of those markers for all
#' of the barley lines; 3) \code{geno_mat} - a
#' \code{matrix} with 1590 rows and 2187 columns. The elements of the matrix contain
#' the genotypic states at each of 2187 markers for each of 1590 barley lines.
#'
#' @source
#' This data was downloaded from the Triticeae Toolbox at https://triticeaetoolbox.org/barley/.
#'
"cap.barley"



#' Example data from the \code{PopVar} package.
#'
#' @description
#' A \code{list} of phenotypic and molecular marker data on a population of 165
#' six-row barley lines.
#'
#' @format An object of class \code{list} with two elements: 1) \code{pheno} - a
#' \code{data.frame} with 165 rows and 5 columns. The first column contains the
#' line name, and the remaining four columns contain phenotypic data on four
#' quantitative traits; 2) \code{geno} - a \code{data.frame} with 742 rows and 168 columns.
#' The first 3 columns contain information on the name, chromosome, and genetic
#' map position (respectively) of 742 molecular markers. Remaining columns
#' contain information on the genotypic states at each of those markers for all
#' of the barley lines.
#'
#' @source
#' This data was obtained and lightly edited from the \code{PopVar} package.
#'
"popvar.barley"
