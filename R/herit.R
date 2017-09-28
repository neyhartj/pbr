#' Calculate heritability from fitted model
#'
#' @description
#' Calculates the heritability (on an entry-mean basis) of a trait given a fitted model object
#'
#' @param object A model object. See \emph{Details} for accepted classes.
#' @param exp A quoted expression used to calculate the heritability. For instance,
#' \code{"geno / (geno + (Residual / r))"}.
#' @param ... Other arguments to pass. This is generally a list of other objects
#' that are in \code{exp}, but may not be found in \code{object}. For instance, if
#' the argument \code{exp = "geno / (geno + (Residual / r))"} was passed, you would
#' also pass the argument \code{r = 2}.
#'
#' @details
#'
#' Accepted classes for \code{object} are:
#'
#' \itemize{
#'   \item{\code{lm} - CURRENTLY UNAVAILABLE}
#'   \item{\code{lmerMod}}
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
#' # Set the number of reps and number of environments
#' n_r <- 4
#' n_e <- 36
#'
#' # Fit a linear model using lm
#' mod <- lmer(yield ~ (1|gen) + env + (1|gen:env), data = gauch_soy1)
#' # Calculate heritability
#' herit(object = mod, exp = "gen / (gen + (gen:env / n_e) + (Residual / n_r))",
#'       n_r = n_r, n_e = n_e)
#'
#' @import dplyr
#' @import lme4
#' @importFrom stringr str_replace_all
#'
#' @export
#'
herit <- function(object, exp, ...) {

  UseMethod(generic = "herit", object = object)

} # Close the function



# herit.lm <- function(object, exp, ...) {
#
#   # Remove colons from the expression
#   exp1 <- str_replace_all(string = exp, pattern = ":", replacement = "")
#
#   # Get variance components
#   anova_table <- as.data.frame(anova(object))
#
#   # Capture the other arguments
#   other_args <- list(...)
#
#   # Create a new environment to evaluate expressions and assign values
#   exp_env <- new.env(parent = parent.frame())
#
#   # Extract other terms and create objects
#   for (obj in names(other_args)) assign(x = obj, value = other_args[[obj]], envir = exp_env)
#
#   # Parse the equation text
#   exp_parsed <- parse(text = exp1)
#
#   # Extract the variance components from the model and assign to individual objects
#   for (term in var_comp$grp) {
#     # Remove colons
#     term_name <- str_replace_all(string = term, pattern = ":", replacement = "")
#
#     assign(x = term_name, value = subset(var_comp, grp == term, vcov, drop = TRUE),
#            envir = exp_env)
#
#   }
#
#   # Evaluate the exp expression and return the heritability
#   with(exp_env, eval(exp_parsed))
#
#
#
#
#
#   ## Error handling
#   # Pull out the terms
#   object_terms <- attr(terms(object), "term.labels")
#
#   arg_terms <- c(geno.term, env.term, ge.term)
#
#   # Make sure the terms in the arguments are in the model
#   args_not_in_model <- setdiff(arg_terms, object_terms)
#   # Error out
#   if (length(args_not_in_model) > 0)
#     stop("The argument(s) '", paste0(args_not_in_model, collapse = ", "), "' were in the arguments but not in the model object.")
#
#   # Make sure all the object terms are in the arguments
#   terms_not_in_args <- setdiff(object_terms, arg_terms)
#   # Send a warning
#   if (length(terms_not_in_args) > 0)
#     warning("The term(s) '", paste0(terms_not_in_args, collapse = ", "), "' were in the model object but not in the arguments.")
#
#   # First pull out the model frame
#   object_mf <- model.frame(object)
#
#   # Get anova table
#   anova_table <- anova(object)
#
#   # Number of unique genos and envs
#   n_geno <- n_distinct(object_mf[,geno.term])
#
#   ## If env.term is NULL, the number of environments is 1 and the number of gen-env
#   ## combinations is the same as the number of genotypes
#   if (!is.null(env.term)) {
#
#     # Get the estimate of the mean squares
#     MS_geno <- anova_table[geno.term, "Mean Sq"]
#     MS_ge <- anova_table[ge.term, "Mean Sq"]
#     MS_error <- anova_table["Residuals", "Mean Sq"]
#
#     ## Determine the number of replicates per genotype per environment
#     # Table of reps per genotype per environment
#     rep_table <- table(object_mf[,c(geno.term, env.term)])
#
#     # Count the number of environments per genotype
#     e_j <- apply(rep_table, MARGIN = 1, FUN = function(geno) sum(geno > 0))
#
#     # Harmonic means of plots and environments
#     e_h <- harm_mean(e_j)
#     p_h <- harm_mean(rowSums(rep_table))
#
#     # Was V_GE empty? Correctly calculate heritability
#     if (length(MS_ge) == 0) {
#
#       # Estimate the genetic variance
#       V_G <- (MS_geno - MS_error) / (p_h)
#       # Residual variance
#       V_R <- MS_error
#
#       V_GE <- NA
#
#       h <- V_G / (V_G + (V_R / p_h))
#
#     } else {
#       # Calculate normally
#
#       # Estimate the genetic variance
#       V_G <- (MS_geno - MS_ge) / (p_h * e_h)
#       # Estimate GxE variance
#       V_GE <- (MS_ge - MS_error) / p_h
#       # Residual variance
#       V_R <- MS_error
#
#       h <- V_G / (V_G + (V_GE / e_h) + (V_R / (p_h * e_h)))
#
#     }
#
#   } else {
#     # Estimate heritability within one environment
#
#     # Get the estimate of the mean squares
#     MS_geno <- anova_table[geno.term, "Mean Sq"]
#     MS_error <- anova_table["Residuals", "Mean Sq"]
#
#     # Table of reps per genotype per environment
#     rep_table <- table(object_mf[,c(geno.term, env.term)])
#
#     p_h <- harm_mean(rowSums(rep_table))
#
#     # Estimate the genetic variance
#     V_G <- (MS_geno - MS_error) / p_h
#     # V_GE is NULL
#     V_GE <- NA
#     # Residual variance
#     V_R <- MS_error
#
#     # Heritability on an entry-mean basis
#     h <- V_G / (V_G + (V_R / p_h))
#
#   }
#
#   # Assemble data.frame and output
#   data.frame(Term = c("heritability", "V_G", "V_GE", "V_R"), Estimate = c(h, V_G, V_GE, V_R))
#
#
# } # Close the function


#' @rdname herit
#' @export
herit.lmerMod <- function(object, exp, ...) {

  # Remove colons from the expression
  exp1 <- str_replace_all(string = exp, pattern = ":", replacement = "")

  # Get variance components
  var_comp <- as.data.frame(VarCorr(object))

  # Capture the other arguments
  other_args <- list(...)

  # Create a new environment to evaluate expressions and assign values
  exp_env <- new.env(parent = parent.frame())

  # Extract other terms and create objects
  for (obj in names(other_args)) assign(x = obj, value = other_args[[obj]], envir = exp_env)

  # Parse the equation text
  exp_parsed <- parse(text = exp1)

  # Extract the variance components from the model and assign to individual objects
  for (term in var_comp$grp) {
    # Remove colons
    term_name <- str_replace_all(string = term, pattern = ":", replacement = "")

    assign(x = term_name, value = subset(var_comp, grp == term, vcov, drop = TRUE),
                                    envir = exp_env)

  }

  # Assign the parsed expression to the environment
  assign(x = "exp_parsed", value = exp_parsed, envir = exp_env)

  # Evaluate the exp expression and return the heritability
  eval(expr = exp_parsed, envir = exp_env)

}
