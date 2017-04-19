#'
#' @describeIn herit
#'
herit.lm <- function(object, geno.term, env.term = NULL, ge.term = NULL) {

  ## Error handling
  # Pull out the terms
  object_terms <- attr(terms(object), "term.labels")

  arg_terms <- c(geno.term, env.term, ge.term)

  # Make sure the terms in the arguments are in the model
  args_not_in_model <- setdiff(arg_terms, object_terms)
  # Error out
  if (length(args_not_in_model) > 0)
    stop("The argument(s) '", args_not_in_model, "' were in the arguments but not in the model object.")

  # Make sure all the object terms are in the arguments
  terms_not_in_args <- setdiff(object_terms, arg_terms)
  # Error out
  if (length(terms_not_in_args) > 0)
    stop("The term(s) '", terms_not_in_args, "' were in the model object but not in the arguments.")

  # First pull out the model frame
  object_mf <- model.frame(object)

  # Get anova table
  anova_table <- anova(object)

  # Number of unique genos and envs
  n_geno <- n_distinct(object_mf[,geno.term])

  ## If env.term is NULL, the number of environments is 1 and the number of gen-env
  ## combinations is the same as the number of genotypes
  if (!is.null(env.term)) {

    n_env <- n_distinct(object_mf[,env.term])

    # Get the estimate of the mean squares
    MS_geno <- anova_table[geno.term, "Mean Sq"]
    MS_ge <- anova_table[ge.term, "Mean Sq"]
    MS_error <- anova_table["Residuals", "Mean Sq"]

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    # Calculate the average number of reps
    mu_rep <- avg_reps(rep_table)

    # Was V_GE empty? Correctly calculate heritability
    if (length(MS_ge) == 0) {

      # Estimate the genetic variance
      V_G <- (MS_geno - MS_error) / (mu_rep)
      # Residual variance
      V_R <- MS_error

      V_GE <- NA

      h <- V_G / (V_G + (V_R / mu_rep))

    } else {
      # Calculate normally

      # Estimate the genetic variance
      V_G <- (MS_geno - MS_ge) / (mu_rep * n_env)
      # Estimate GxE variance
      V_GE <- (MS_ge - MS_error) / mu_rep
      # Residual variance
      V_R <- MS_error

      h <- V_G / (V_G + (V_GE / n_env) + (V_R / (mu_rep * n_env)))

    }

  } else {
    # Estimate heritability within one environment

    # Get the estimate of the mean squares
    MS_geno <- anova_table[geno.term, "Mean Sq"]
    MS_error <- anova_table["Residuals", "Mean Sq"]

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term)])

    # Calculate the average number of reps
    mu_rep <- avg_reps(rep_table)

    # Estimate the genetic variance
    V_G <- (MS_geno - MS_error) / mu_rep
    # V_GE is NULL
    V_GE <- NA
    # Residual variance
    V_R <- MS_error

    # Heritability on an entry-mean basis
    h <- V_G / (V_G + (V_R / mu_rep))

  }

  # Assemble data.frame and output
  data.frame(Term = c("heritability", "V_G", "V_GE", "V_R"), Estimate = c(h, V_G, V_GE, V_R))


} # Close the function
#'
#' @describeIn herit
#'
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
  # Error out
  if (length(terms_not_in_args) > 0)
    stop("The term(s) '", paste0(terms_not_in_args, collapse = ", "), "' were in the model object but not in the arguments.")

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

    ## Determine the number of replicates per genotype per environment
    # Table of reps per genotype per environment
    rep_table <- table(object_mf[,c(geno.term, env.term)])

    # Calculate the average number of reps
    mu_rep <- avg_reps(rep_table)
    # Number of environments
    n_env <- n_distinct(object_mf[,env.term])

    # Was V_GE empty? Correctly calculate heritability
    if (length(V_GE) == 0) {
      V_GE <- NA

      h <- V_G / (V_G + (V_R / mu_rep))

    } else {
      # Calculate normally

      h <- V_G / (V_G + (V_GE / n_env) + (V_R / (mu_rep * n_env)))

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
    rep_table <- table(object_mf[,c(geno.term)])

    # Calculate the average number of reps
    mu_rep <- avg_reps(rep_table)

    # Calculate heritability
    h <- V_G / (V_G + (V_R / mu_rep))

  }

  # Assemble data.frame and output
  data.frame(Term = c("heritability", "V_G", "V_GE", "V_R"), Estimate = c(h, V_G, V_GE, V_R))


} # Close the function












