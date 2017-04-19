#' Calculate heritability from fitted model
#'
#' @description
#' Calculates the heritability (on an entry-mean basis) of a trait given a fitted model object
#'
#' @param object A model object of class \code{lm}, \code{aov}, \code{glm}, or \code{lmer}.
#' @param geno.term A character vector giving the model term for genotypes.
#' @param env.term A character vector giving the model term for environments. If \code{NULL} (default),
#' heritability is calculated within one environment. Note that \code{env.term = NULL} will
#' ignore multiple environments, if present.
#' @param ge.term A character vector giving the model term for genotype-by-environment interaction.
#'
#'
#' @examples
#'
#' # Use the gauch.soy dataset
#' data("gauch.soy")
#'
#' # Filter
#' gauch_soy1 <- gauch.soy %>%
#'   group_by(env) %>%
#'   filter(n_distinct(gen, rep) == 28)
#'
#' # Fit a linear model using lm
#' mod <- lm(yield ~ gen + env + gen:env, data = gauch_soy1)
#' # Calculate heritability
#' herit(object = mod, geno.term = "gen", env.term = "env", ge.term = "gen:env")
#'
#' # Fit a mixed-effects model with lmer
#' mod1 <- lmer(yield ~ (1|gen) + env + (1|gen:env), data = gauch_soy1)
#' # Calculate heritability
#' herit(object = mod1, geno.term = "gen", env.term = "env", ge.term = "gen:env")
#'
#' # Exclude the GE term
#' mod2 <- lm(yield ~ gen + env, data = gauch_soy1)
#' herit(object = mod2, geno.term = "gen", env.term = "env")
#'
#' # A different dataset with only one rep
#' data("cornelius.maize")
#' mod3 <- lm(yield ~ gen + env, data = cornelius.maize)
#' herit(object = mod3, geno.term = "gen", env.term = "env")
#'
#' # One environment and one rep
#' cornelius_maize <- cornelius.maize %>%
#'   filter(env == "E01")
#'
#' mod4 <- lm(yield ~ gen, data = cornelius_maize)
#' herit(object = mod4, geno.term = "gen")
#'
#' # Note it errors because there are no residuals - cannot calculate heritability
#'
#' @import dplyr
#' @import lme4
#'
#' @export
#'
herit <- function(object, ...) {

  UseMethod(generic = "herit", object = object)

} # Close the function
