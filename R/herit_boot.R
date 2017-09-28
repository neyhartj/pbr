#' Heritability confidence interval by bootstrap
#'
#' @param object A model object. See \emph{Details} for accepted classes.
#' @param exp A quoted expression used to calculate the heritability. For instance,
#' \code{"geno / (geno + (Residual / r))"}.
#' @param boot.reps The number of bootsrap replicates.
#' @param alpha The significance level for the confidence interval.
#' @param ... Other arguments to pass. This is generally a list of other objects
#' that are in \code{exp}, but may not be found in \code{object}. For instance, if
#' the argument \code{exp = "geno / (geno + (Residual / r))"} was passed, you would
#' also pass the argument \code{r = 2}.
#'
#' @importFrom purrr map
#'
#' @export
#'
herit_boot <- function(object, exp, boot.reps = 1000, ...) {

  UseMethod(generic = "herit_boot", object)

}


# herit_boot.lm <- function(object, geno.term, env.term = NULL, ge.term = NULL,
#                           boot.reps = 1000, alpha = 0.95) {
#
#   # Grab the model frame
#   object_mf <- model.frame(object)
#
#   # Grab the model formula
#   object_form <- formula(object)
#
#   # Sample bootstrap replicates
#   boot_df <- bootstrap(data = object_mf, n = boot.reps)
#
#   # Iterate over samples
#   boot_herit <- map(boot_df$strap, ~ lm(object_form, data = .)) %>%
#     map(herit, geno.term = geno.term, env.term = env.term, ge.term = ge.term)
#
#   # Convert to data.frame
#   boot_herit_df <- boot_herit %>%
#     lapply("[[", 2) %>%
#     data.frame()
#
#   # Rename
#   names(boot_herit_df) <- paste("rep", seq(ncol(boot_herit_df)), sep = "")
#
#   # Mutate and gather
#   boot_herit_df1 <- boot_herit_df %>%
#     mutate(Term = boot_herit[[1]]$Term) %>%
#     gather(rep, estimate, -Term)
#
#   # Summarize and return
#   boot_herit_df1 %>%
#     group_by(Term) %>%
#     summarize(mean = mean(estimate), ci_lower = quantile(estimate, 1 - alpha), ci_upper = quantile(estimate, alpha)) %>%
#     as.data.frame()
#
# } # Close function


#' @rdname herit_boot
#' @export
herit_boot.lmerMod <- function(object, exp, boot.reps = 1000, ...) {

  # Calculate heritability
  base_herit <- herit(object = object, exp = exp, ... = ...)

  # Extract other arguments
  other_args <- list(...)

  # Extract the model frame
  mf <- model.frame(object)

  # Simulate responses
  sim_response <- simulate(object = object, nsim = boot.reps)

  # Map over the simulations and refit the model
  sim_fits <- sim_response %>%
    map(refit, object = object)

  # Calculate heritability
  sim_herit <- sim_fits %>%
    map_dbl(herit, exp = exp, ... = ...)

  # Caculate the standard error
  se <- sd(sim_herit)

  # Return the heritability and the standard error
  data.frame(heritability = base_herit, se = se)

} # Close the function


