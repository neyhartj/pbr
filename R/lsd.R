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
#' lsd(object = mod, geno.term = "gen", env.term = "env", ge.term = "gen:env", alpha = 0.05)
#'
#' # Fit a mixed-effects model with lmer
#' mod1 <- lmer(yield ~ (1|gen) + env + (1|gen:env), data = gauch_soy1)
#' # Calculate heritability
#' lsd(object = mod1, geno.term = "gen", env.term = "env", ge.term = "gen:env")
#'
#' @import dplyr
#' @import lme4
#'
#' @export
#'
#'
lsd <- function(object, ...) {

  UseMethod(generic = "lsd", object = object)

} # Close the function





#' @rdname lsd
#' @export
lsd.lm <- function(object, geno.term, env.term = NULL, ge.term = NULL, alpha = 0.05) {

  ## Error handling
  # alpha must be between 0 and 1
  if (!(alpha > 0 & alpha < 1))
    stop("The argument 'alpha' must be greater than 0 and less than 1.")

  # Handle the term arguments
  geno.term <- as.character(geno.term)
  env.term <- as.character(env.term)
  ge.term <- as.character(ge.term)


  ## Run an anova and convert to a data.frame
  # Get variance components
  anova_table <- as.data.frame(anova(object))
  # Add the terms
  anova_table <- cbind(term = row.names(anova_table), anova_table)

  # Are the terms in the argument in the anova table
  stopifnot(geno.term %in% anova_table$term)
  stopifnot(env.term %in% anova_table$term)
  stopifnot(ge.term %in% anova_table$term)

  # First pull out the model frame
  object_mf <- model.frame(object)

  # Number of unique genos
  n_geno <- n_distinct(object_mf[,geno.term])

  ## If env.term is NULL, the number of environments is 1 and the number of gen-env
  ## combinations is the same as the number of genotypes
  if (!is.null(env.term)) {

    # Use the model frame to find the harmonic number of environments and reps
    form <- as.formula(paste("~", geno.term, "+", env.term))
    plot_table <- xtabs(form, object_mf)

    n_rep <- harm_mean(plot_table)
    n_env <- harm_mean(apply(X = ifelse(test = plot_table > 1, 1, 0), MARGIN = 1, FUN = sum))

    ## Extract mean squares estimates of sources of variance
    MS_GE <- subset(anova_table, term == ge.term, `Mean Sq`, drop = TRUE)
    V_R <- subset(anova_table, term == "Residuals", `Mean Sq`, drop = TRUE)

    # Was MS_GE empty? Correctly calculate LSD
    if (is.na(MS_GE)) {

      # Calculate the LSD at the alpha level
      t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1))

      LSD <- t_alpha * sqrt(2 * (V_R / n_rep) )

    } else {
      # Calculate normally

      # Calculate the LSD at the alpha level
      t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1) * (n_env - 1))

      LSD <- t_alpha * sqrt((2 * MS_GE) / (n_rep * n_env))

    }

    # Otherwise calculate the LSD without V_GE
  } else {

    # Use the model frame to find the harmonic number of reps and reps
    form <- as.formula(paste("~", geno.term))
    plot_table <- xtabs(form, object_mf)

    n_rep <- harm_mean(plot_table)

    ## Extract mean squares estimates of sources of variance
    V_R <- subset(anova_table, term == "Residuals", `Mean Sq`, drop = TRUE)

    # Calculate the LSD at the alpha level
    t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1))

    LSD <- t_alpha * sqrt(2 * (V_R / n_rep) )

  }

  # Return the LSD estimate
  return(LSD)


} # Close the function



#' @rdname lsd
#' @export
lsd.lmerMod <- function(object, geno.term, env.term = NULL, ge.term = NULL, alpha = 0.05) {

  ## Error handling
  # alpha must be between 0 and 1
  if (!(alpha > 0 & alpha < 1))
    stop("The argument 'alpha' must be greater than 0 and less than 1.")

  # Handle the term arguments
  geno.term <- as.character(geno.term)
  env.term <- as.character(env.term)
  ge.term <- as.character(ge.term)


  ## Get the variance estimates
  vcor <- as.data.frame(VarCorr(object))

  # Are the terms in the argument in the anova table
  stopifnot(ge.term %in% vcor$grp)

  # First pull out the model frame
  object_mf <- model.frame(object)

  # Number of unique genos
  n_geno <- n_distinct(object_mf[,geno.term])


  ## If env.term is NULL, the number of environments is 1 and the number of gen-env
  ## combinations is the same as the number of genotypes
  if (!is.null(env.term)) {

    # Use the model frame to find the harmonic number of environments and reps
    form <- as.formula(paste("~", geno.term, "+", env.term))
    plot_table <- xtabs(form, object_mf)

    n_rep <- harm_mean(plot_table)
    n_env <- harm_mean(apply(X = ifelse(test = plot_table > 1, 1, 0), MARGIN = 1, FUN = sum))

    ## Extract mean squares estimates of sources of variance
    V_GE <- subset(vcor, grp == ge.term, vcov, drop = TRUE)
    V_R <- subset(vcor, grp == "Residual", vcov, drop = TRUE)

    # Was V_GE empty? Correctly calculate LSD
    if (is.na(V_GE)) {

      # Calculate the LSD at the alpha level
      t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1))

      LSD <- t_alpha * sqrt(2 * (V_R / n_rep) )

    } else {
      # Calculate normally

      # Calculate the LSD at the alpha level
      t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1) * (n_env - 1))

      LSD <- t_alpha * sqrt(2 * ( (V_R / (n_rep * n_env)) + (V_GE / n_env) ))

    }

    # Otherwise calculate the LSD without V_GE
  } else {

    # Use the model frame to find the harmonic number of reps and reps
    form <- as.formula(paste("~", geno.term))
    plot_table <- xtabs(form, object_mf)

    n_rep <- harm_mean(plot_table)

    ## Extract the residual variance
    V_R <- subset(vcor, grp == "Residual", vcov, drop = TRUE)

    # Calculate the LSD at the alpha level
    t_alpha <- qt(p = 1 - (alpha / 2), df = (n_geno - 1))

    LSD <- t_alpha * sqrt(2 * (V_R / n_rep) )

  }

  # Return the LSD estimate
  return(LSD)

} # Close the function
