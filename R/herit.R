#' Calculate heritability from fitted model
#'
#' @description
#' Calculates the heritability (on an entry-mean basis) of a trait given a fitted model object
#'
#' @param object A model object. See \emph{Details} for accepted classes. of class \code{lm}, \code{aov}, \code{glm}, or \code{lmer}.
#' @param geno.term A character vector giving the model term for genotypes.
#' @param env.term A character vector giving the model term for environments. If \code{NULL} (default),
#' heritability is calculated within one environment. Note that \code{env.term = NULL} will
#' ignore multiple environments, if present.
#' @param ge.term A character vector giving the model term for genotype-by-environment interaction.
#'
#' @details
#'
#' Accepted classes for \code{object} are:
#'
#' \itemize{
#'   \item{\code{lm}}
#'   \item{\code{aov}}
#'   \item{\code{glm}}
#'   \item{\code{lmer}}
#'   \item{\code{mmer}}
#' }
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



#' @rdname herit
#' @export
herit.lm <- function(object, geno.term, env.term = NULL, ge.term = NULL) {

  ## Error handling
  # Pull out the terms
  object_terms <- attr(terms(object), "term.labels")

  arg_terms <- c(geno.term, env.term, ge.term)

  # Make sure the terms in the arguments are in the model
  args_not_in_model <- setdiff(arg_terms, object_terms)
  # Error out
  if (length(args_not_in_model) > 0)
    stop("The argument(s) '", paste0(args_not_in_model, collapse = ", "), "' were in the arguments but not in the model object.")

  # Make sure all the object terms are in the arguments
  terms_not_in_args <- setdiff(object_terms, arg_terms)
  # Send a warning
  if (length(terms_not_in_args) > 0)
    warning("The term(s) '", paste0(terms_not_in_args, collapse = ", "), "' were in the model object but not in the arguments.")

  # First pull out the model frame
  object_mf <- model.frame(object)

  # Get anova table
  anova_table <- anova(object)

  # Number of unique genos and envs
  n_geno <- n_distinct(object_mf[,geno.term])

  ## If env.term is NULL, the number of environments is 1 and the number of gen-env
  ## combinations is the same as the number of genotypes
  if (!is.null(env.term)) {

    # Get the estimate of the mean squares
    MS_geno <- anova_table[geno.term, "Mean Sq"]
    MS_ge <- anova_table[ge.term, "Mean Sq"]
    MS_error <- anova_table["Residuals", "Mean Sq"]

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    # Count the number of environments per genotype
    e_j <- apply(rep_table, MARGIN = 1, FUN = function(geno) sum(geno > 0))

    # Harmonic means of plots and environments
    e_h <- harm_mean(e_j)
    p_h <- harm_mean(rep_table)

    # Was V_GE empty? Correctly calculate heritability
    if (length(MS_ge) == 0) {

      # Estimate the genetic variance
      V_G <- (MS_geno - MS_error) / (p_h)
      # Residual variance
      V_R <- MS_error

      V_GE <- NA

      h <- V_G / (V_G + (V_R / p_h))

    } else {
      # Calculate normally

      # Estimate the genetic variance
      V_G <- (MS_geno - MS_ge) / (p_h * e_h)
      # Estimate GxE variance
      V_GE <- (MS_ge - MS_error) / p_h
      # Residual variance
      V_R <- MS_error

      h <- V_G / (V_G + (V_GE / e_h) + (V_R / (p_h * e_h)))

    }

  } else {
    # Estimate heritability within one environment

    # Get the estimate of the mean squares
    MS_geno <- anova_table[geno.term, "Mean Sq"]
    MS_error <- anova_table["Residuals", "Mean Sq"]

    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    p_h <- harm_mean(rep_table)

    # Estimate the genetic variance
    V_G <- (MS_geno - MS_error) / p_h
    # V_GE is NULL
    V_GE <- NA
    # Residual variance
    V_R <- MS_error

    # Heritability on an entry-mean basis
    h <- V_G / (V_G + (V_R / p_h))

  }

  # Assemble data.frame and output
  data.frame(Term = c("heritability", "V_G", "V_GE", "V_R"), Estimate = c(h, V_G, V_GE, V_R))


} # Close the function


#' @rdname herit
#' @export
herit.lmerMod <- function(object, geno.term, env.term = NULL, ge.term = NULL) {

  ## Error handling
  # Pull out the terms
  fixed_terms <- attr(terms(object), "term.labels")
  random_terms <- names(object@cnms)

  object_terms <- c(fixed_terms, random_terms)

  # Make sure the terms in the arguments are in the model
  arg_terms <- c(geno.term, env.term, ge.term)

  # Make sure the terms in the arguments are in the model
  args_not_in_model <- setdiff(arg_terms, object_terms)
  # Error out
  if (length(args_not_in_model) > 0)
    stop("The argument(s) '", paste0(args_not_in_model, collapse = ", "), "' were in the arguments but not in the model object.")

  # Make sure all the object terms are in the arguments
  terms_not_in_args <- setdiff(object_terms, arg_terms)
  # Send a warning
  if (length(terms_not_in_args) > 0)
    warning("The term(s) '", paste0(terms_not_in_args, collapse = ", "), "' were in the model object but not in the arguments.")

  # First pull out the model frame
  object_mf <- model.frame(object)

  # Get anova table
  anova_table <- anova(object)
  # Get the VarCorr table
  vcor <- as.data.frame(VarCorr(object))

  # Number of unique genos and envs
  n_geno <- n_distinct(object_mf[,geno.term])

  ## If env.term is NULL, the number of environments is 1 and the number of gen-env
  ## combinations is the same as the number of genotypes
  if (!is.null(env.term)) {

    # Get the estimates of the variances of the random effects
    V_G <- subset(vcor, grp == geno.term, vcov, drop = TRUE)
    V_GE <- subset(vcor, grp == ge.term, vcov, drop = TRUE)
    V_R <- subset(vcor, grp == "Residual", vcov, drop = TRUE)

    ## Calculate the harmonic mean of the number of environments


    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    # Count the number of environments per genotype
    e_j <- apply(rep_table, MARGIN = 1, FUN = function(geno) sum(geno > 0))

    # Harmonic means of plots and environments
    e_h <- harm_mean(e_j)
    p_h <- harm_mean(rep_table)

    # Was V_GE empty? Correctly calculate heritability
    if (length(V_GE) == 0) {
      V_GE <- NA

      h <- V_G / (V_G + (V_R / p_h))

    } else {
      # Calculate normally

      h <- V_G / (V_G + (V_GE / e_h) + (V_R / (p_h * e_h)))

    }

  } else {
    # Estimate heritability within one environment

    # Get the estimates of the variances of the random effects
    V_G <- subset(vcor, grp == geno.term, vcov, drop = TRUE)
    V_R <- subset(vcor, grp == "Residual", vcov, drop = TRUE)

    # V_GE is NA
    V_GE <- NA

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    p_h <- harm_mean(rep_table)

    # Calculate heritability
    h <- V_G / (V_G + (V_R / p_h))

  }

  # Assemble data.frame and output
  data.frame(Term = c("heritability", "V_G", "V_GE", "V_R"), Estimate = c(h, V_G, V_GE, V_R))


} # Close the function
