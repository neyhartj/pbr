#' Genomewide Assocation Analysis Through Linear Mixed-Models
#'
#' @description
#' Performs a genomewide association analysis for the main effect of QTL. Most
#' arguments are taken from the \code{\link[GWAS]{rrBLUP}} function.
#'
#' @param pheno A data.frame of phenotypic data.
#' @param geno A data.frame of marker names, positions, and genotypes.
#' @param fixed A character vector with the column names of fixed effects. All
#' other columns will be assumed traits.
#' @param K_g The covariance matrix of random genotype effects (i.e. effect of
#' background markers). If \code{NULL}, it is determined from the additive relationship
#' matrix of the \code{geno} input. If missing, the effect of background markers
#' is not included in the model.
#' @param n.PC The number of principal components from singular value decomposition
#' of the \code{K_g} matrix to use to correct for population structure. If \code{0},
#' no correction for population structure is performed.
#' @param P3D Logical. Population parameters previous determined.
#' @param n.core The number of cores to use when calculating marker scores.
#' @param impute.method The method by which to impute the marker genotypes.
#'
#' @return
#' An object of class \code{gwas}, with marker scores for each trait.
#'
#'


