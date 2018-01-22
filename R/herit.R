#' Calculate heritability from fitted model
#'
#' @description
#' Calculates the heritability (on an entry-mean basis) of a trait given a fitted model object
#'
#' @param object A model object. See \emph{Details} for accepted classes.
#' @param exp A quoted expression used to calculate the heritability. For instance,
#' \code{"geno / (geno + (Residual / r))"}. Variables in the expression should be variance
#' components derived from restricted maximum likelihood (REML) or from the mean squares
#' from an ANOVA table.
#' @param ms_exp A named list of expressions used to calculate the variance components.
#' The names of the list should be names of the variance components used in the heritability
#' expression (\code{exp}). Note the expected residual variance is equal to the observed
#' mean squares for the residuals. See \emph{Examples} for examples.
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
#'   \item{\code{lm}}
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
#' lm_mod <- lm(yield ~ gen + env + gen:env + rep %in% env, data = gauch_soy1)
#'
#' # Variance components from a fixed effects model are derived from the ANOVA table.
#' # The function also required expressions to calculate the variance components
#' ms_exp <- list("gen:env" = "(gen:env - Residuals) / n_r",
#'                "gen" = "(gen - gen:env) / (n_r * n_e)")
#'
#' exp = "gen / (gen + (gen:env / n_e) + (Residuals / n_r))"
#'
#' herit(object = lm_mod, exp = exp, ms_exp = ms_exp, n_r = n_r, n_e = n_e)
#'
#' # Fit a linear model using lmer
#' lmer_mod <- lmer(yield ~ (1|gen) + env + (1|gen:env) + (1|rep:env), data = gauch_soy1)
#' # Calculate heritability
#' herit(object = lmer_mod, exp = "gen / (gen + (gen:env / n_e) + (Residual / n_r))",
#'       n_r = n_r, n_e = n_e)
#'
#' @import dplyr
#' @import lme4
#' @importFrom stringr str_replace_all
#' @importFrom broom tidy
#'
#' @export
#'
herit <- function(object, exp, ...) {

  UseMethod(generic = "herit", object = object)

} # Close the function


herit.lm <- function(object, exp, ms_exp, ...) {

  # Remove colons from the expressions
  ms_exp1 <- lapply(ms_exp, str_replace_all, ":", "")
  exp1 <- str_replace_all(string = exp, pattern = ":", replacement = "")

  # Get variance components
  anova_table <- tidy(anova(object))

  # Make sure the terms are in the anova table
  if(!all(names(ms_exp1) %in% anova_table$term))
    stop("The names in 'ms_exp' do not match the terms in the model.")

  # Change the names
  names(ms_exp1) <- str_replace_all(string = names(ms_exp1), pattern = ":", replacement = "")

  # Capture the other arguments
  other_args <- list(...)

  # Parse if necessary
  if (all(length(other_args) == 1, sapply(other_args, class) == "list")) {
    other_args <- do.call("c", other_args)
  }

  # Parse the equation text
  ms_exp_parsed <- lapply(X = ms_exp1, function(x) parse(text = x))
  exp_parsed <- parse(text = exp1)


  # Create a new environment to evaluate expressions and assign values
  exp_env <- new.env(parent = parent.frame())

  # Extract other terms and create objects
  for (obj in names(other_args)) assign(x = obj, value = other_args[[obj]], envir = exp_env)

  # Extract the mean squares and send to the environment
  for (tm in anova_table$term) {
    # Remove colons
    term_name <- str_replace_all(string = tm, pattern = ":", replacement = "")

    assign(x = term_name, value = subset(anova_table, term == tm, meansq, drop = TRUE),
                                      envir = exp_env)

  }

  # Evaluate the mean squares expressions
  ms_eval <- lapply(ms_exp_parsed, eval, envir = exp_env)

  # Assign names
  for (tm in names(ms_eval))  assign(x = tm, value = ms_eval[[tm]], envir = exp_env)

  # Evaluate the exp expression and return the heritability
  eval(exp_parsed, envir = exp_env)

} # Close the function


#' @rdname herit
#' @export
herit.lmerMod <- function(object, exp, ...) {

  # Remove colons from the expression
  exp1 <- str_replace_all(string = exp, pattern = ":", replacement = "")

  # Get variance components
  var_comp <- as.data.frame(VarCorr(object))

  # Capture the other arguments
  other_args <- list(...)

  # Parse if necessary
  if (all(length(other_args) == 1, sapply(other_args, class) == "list")) {
    other_args <- do.call("c", other_args)
  }

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
