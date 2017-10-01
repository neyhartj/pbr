#' Calculate the LSD using a fitted model object
#'
#' @description
#' Calculate the least-significant difference among the means of genotypes from
#' a fitted model object.
#'
#' @param object A model object of class \code{lm}, \code{aov}, \code{glm}, or \code{lmer}.
#' @param geno.term A character vector giving the model term for genotypes.
#' @param env.term A character vector giving the model term for environments. If \code{NULL} (default),
#' heritability is calculated within one environment. Note that \code{env.term = NULL} will
#' ignore multiple environments, if present.
#' @param ge.term A character vector giving the model term for genotype-by-environment interaction.
#' @param alpha The significance level at which to calculate the LSD.
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
#' LSD(object = mod, geno.term = "gen", env.term = "env", ge.term = "gen:env", alpha = 0.05)
#'
#' # Fit a mixed-effects model with lmer
#' mod1 <- lmer(yield ~ (1|gen) + env + (1|gen:env), data = gauch_soy1)
#' # Calculate heritability
#' LSD(object = mod1, geno.term = "gen", env.term = "env", ge.term = "gen:env")
#'
#' @import dplyr
#' @import lme4
#'
#' @export
#'
#'
LSD <- function(object, ...) {

  UseMethod(generic = "LSD", object = object)

} # Close the function





#' @rdname herit
#' @export
LSD.lm <- function(object, geno.term, env.term = NULL, ge.term = NULL, alpha = 0.05) {

  ## Error handling
  # alpha must be between 0 and 1
  if (!(alpha > 0 & alpha < 1))
    stop("The argument 'alpha' must be greater than 0 and less than 1.")


  # Pass the object through the heritability function to estimate the variance
  # components
  herit_out <- herit(object = object, geno.term = geno.term, env.term = env.term, ge.term = ge.term)

  # First pull out the model frame
  object_mf <- model.frame(object)

  # Number of unique genos and envs
  n_geno <- n_distinct(object_mf[,geno.term])

  ## If env.term is NULL, the number of environments is 1 and the number of gen-env
  ## combinations is the same as the number of genotypes
  if (!is.null(env.term)) {

    n_env <- n_distinct(object_mf[,env.term])

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    # Calculate the average number of reps
    mu_rep <- avg_reps(rep_table)

    # Extract the variance estimates
    V_GE <- subset(x = herit_out, Term == "V_GE", Estimate, drop = TRUE)
    V_R <- subset(x = herit_out, Term == "V_R", Estimate, drop = TRUE)


    # Was V_GE empty? Correctly calculate LSD
    if (is.na(V_GE)) {

      # Calculate the LSD at the alpha level
      t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1))

      LSD <- t_alpha * sqrt(2 * (V_R / mu_rep) )

    } else {
      # Calculate normally

      # Calculate the LSD at the alpha level
      t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1) * (n_env - 1))

      LSD <- t_alpha * sqrt(2 * ( (V_R / (mu_rep * n_env)) + (V_GE / n_env) ) )

    }

  } else {

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    # Calculate the average number of reps
    mu_rep <- avg_reps(rep_table)

    # Extract the variance estimates
    V_R <- subset(x = herit_out, Term == "V_R", Estimate, drop = TRUE)

    # Calculate the LSD at the alpha level
    t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1))

    LSD <- t_alpha * sqrt(2 * (V_R / mu_rep) )

  }

  # Return the LSD estimate
  return(LSD)


} # Close the function



#' @rdname LSD
#' @export
LSD.lmerMod <- function(object, geno.term, env.term = NULL, ge.term = NULL, alpha = 0.05) {

  ## Error handling
  # alpha must be between 0 and 1
  if (!(alpha > 0 & alpha < 1))
    stop("The argument 'alpha' must be greater than 0 and less than 1.")


  # Pass the object through the heritability function to estimate the variance
  # components
  herit_out <- herit(object = object, geno.term = geno.term, env.term = env.term, ge.term = ge.term)

  # First pull out the model frame
  object_mf <- model.frame(object)

  # Number of unique genos and envs
  n_geno <- n_distinct(object_mf[,geno.term])

  ## If env.term is NULL, the number of environments is 1 and the number of gen-env
  ## combinations is the same as the number of genotypes
  if (!is.null(env.term)) {

    n_env <- n_distinct(object_mf[,env.term])

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    # Calculate the average number of reps
    mu_rep <- avg_reps(rep_table)

    # Extract the variance estimates
    V_GE <- subset(x = herit_out, Term == "V_GE", Estimate, drop = TRUE)
    V_R <- subset(x = herit_out, Term == "V_R", Estimate, drop = TRUE)


    # Was V_GE empty? Correctly calculate LSD
    if (is.na(V_GE)) {

      # Calculate the LSD at the alpha level
      t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1))

      LSD <- t_alpha * sqrt(2 * (V_R / mu_rep) )

    } else {
      # Calculate normally

      # Calculate the LSD at the alpha level
      t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1) * (n_env - 1))

      LSD <- t_alpha * sqrt(2 * ( (V_R / (mu_rep * n_env)) + (V_GE / n_env) ) )

    }

  } else {

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    # Calculate the average number of reps
    mu_rep <- avg_reps(rep_table)

    # Extract the variance estimates
    V_R <- subset(x = herit_out, Term == "V_R", Estimate, drop = TRUE)

    # Calculate the LSD at the alpha level
    t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1))

    LSD <- t_alpha * sqrt(2 * (V_R / mu_rep) )

  }

  # Return the LSD estimate
  return(LSD)

} # Close the function
