#' Heritability confidence interval by bootstrap
#'
#' @param object A model object of class \code{lm}, \code{aov}, \code{glm}, or \code{lmer}.
#' @param geno.term A character vector giving the model term for genotypes.
#' @param env.term A character vector giving the model term for environments. If \code{NULL} (default),
#' heritability is calculated within one environment. Note that \code{env.term = NULL} will
#' ignore multiple environments, if present.
#' @param ge.term A character vector giving the model term for genotype-by-environment interaction.
#' @param boot.reps The number of bootsrap replicates.
#' @param alpha The significance level for the confidence interval.
#'
#' @import modelr
#' @importFrom purrr map
#' @import tidyr
#' @import dplyr
#'
#' @export
#'
herit_boot <- function(object, ...) {

  UseMethod(generic = "herit_boot", object)

}
