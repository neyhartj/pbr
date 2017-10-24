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
#' @details
#' This function implements model-based bootstrapping to obtain the standard error
#' and confidence interval around an estimate of the heritability. The function
#' uses the original fitted model to simulate new observations using the same
#' model parameters. The model is then re-fitted using the new observations,
#' and the heritability is re-calculated. This is repeated \emph{n} times.
#'
#' See \code{\link[herit]{pbr}} for accepted classes for \code{herit_boot}.
#'
#' @return
#' A \code{data.frame} with the following values:
#'
#' \describe{
#'   \item{heritability}{The estimate of the heritability using the original data.}
#'   \item{se}{The standard error of the estimate, calculated as the standard deviation
#'   among the bootstrap replicates.}
#'   \item{bias}{The bias of the original heritability estimate, calculated as the difference
#'   between the mean of the bootstrapped estimates and the original estimate.}
#'   \item{ci_lower}{The lower limit of the confidence interval.}
#'   \item{ci_upper}{The upper limit of the confidence interval.}
#' }
#'
#' @importFrom purrr map
#' @importFrom lme4 bootMer
#'
#' @export
#'
herit_boot <- function(object, exp, ms_exp, boot.reps = 1000, alpha = 0.05, ...) {

  UseMethod(generic = "herit_boot", object)

}


herit_boot.lm <- function(object, exp, ms_exp, boot.reps = 1000, alpha = 0.05, ...) {

  # Extract other arguments
  other_args <- list(...)

  # Extract the model frame
  mf <- model.frame(object)

  # Create a function to calculate heritability from the model object
  herit_use <- function(x) herit(object = x, exp = exp, ms_exp = ms_exp, other_args)

  ## Perform bootstrapping
  # First simulate from the model
  sim_response <- simulate(object = object, nsim = boot.reps)

  # Update the model frame with the data
  mf_update <- lapply(X = sim_response, FUN = function(x) {
    mf[,1] <- x
    return(mf)
  })

  # Update and fit the models with the new data
  sim_refit <- lapply(X = mf_update, function(x) update(object, data = x))

  # Calculate heritability on the original data
  base_herit <- herit_use(object)

  # Calculate the heritability for each new model object
  boot_herit <- map_dbl(sim_refit, herit_use)

  # Calculate the standard error
  se <- sd(boot_herit)
  # Calculate the bias
  bias <- mean(boot_herit) - base_herit

  # Calculate the confidence interval
  ci_lower = (2 * base_herit) - quantile(boot_herit, 1 - (alpha / 2))
  ci_upper = (2 * base_herit) - quantile(boot_herit, (alpha / 2))

  # Return the heritability and the standard error
  data.frame(heritability = base_herit, se = se, bias = bias, ci_lower = ci_lower,
             ci_upper = ci_upper, row.names = NULL)

} # Close function


#' @rdname herit_boot
#' @export
herit_boot.lmerMod <- function(object, exp, boot.reps = 1000, alpha = 0.05, ...) {

  # Extract other arguments
  other_args <- list(...)

  # Extract the model frame
  mf <- model.frame(object)

  # Create a function to calculate heritability from the model object
  herit_use <- function(x) herit(object = x, exp = exp, ... = other_args)

  # Perform bootstrapping
  boot_herit <- bootMer(x = object, FUN = herit_use, nsim = boot.reps, type = "parametric")

  # Extract the original fit
  base_herit <- boot_herit$t0
  # Calculate the standard error
  se <- sd(boot_herit$t)
  # Calculate the bias
  bias <- mean(boot_herit$t) - base_herit

  # Calculate the confidence interval
  ci_lower = (2 * base_herit) - quantile(boot_herit$t, 1 - (alpha / 2))
  ci_upper = (2 * base_herit) - quantile(boot_herit$t, (alpha / 2))

  # Return the heritability and the standard error
  data.frame(heritability = base_herit, se = se, bias = bias, ci_lower = ci_lower,
             ci_upper = ci_upper, row.names = NULL)

} # Close the function


