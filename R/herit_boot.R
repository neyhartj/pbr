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
#' #' # Use the gauch.soy dataset
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
#' ms.exp <- list("gen:env" = "(gen:env - Residuals) / n_r",
#'                "gen" = "(gen - gen:env) / (n_r * n_e)")
#'
#' exp = "gen / (gen + (gen:env / n_e) + (Residuals / n_r))"
#'
#' herit(object = lm_mod, exp = exp, ms.exp = ms.exp, n_r = n_r, n_e = n_e)
#'
#' # Fit a linear model using lmer
#' lmer_mod <- lmer(yield ~ (1|gen) + (1|env) + (1|gen:env) + (1|env:rep), data = gauch_soy1)
#' # Calculate heritability
#' herit(object = lmer_mod, exp = "gen / (gen + (gen:env / n_e) + (Residual / n_r))",
#'       n_r = n_r, n_e = n_e)
#'
#'
#' @importFrom lme4 bootMer
#'
#' @export
#'
herit_boot <- function(object, ...) {

  UseMethod(generic = "herit_boot", object = object)

}


herit_boot.lm <- function(object, exp, ms.exp, boot.reps = 1000, alpha = 0.05, ...) {

  # Extract other arguments
  other_args <- list(...)

  # Parse if necessary
  if (all(length(other_args) == 1, sapply(other_args, class) == "list")) {
    other_args <- do.call("c", other_args)
  }

  # Extract the model frame
  mf <- model.frame(object)

  # Create a function to calculate heritability from the model object
  herit_use <- function(x) herit(object = x, exp = exp, ms.exp = ms.exp, other_args)

  ## Perform bootstrapping
  # First simulate from the model
  sim_response <- simulate(object = object, nsim = boot.reps)

  # Get the column position of the response variable in the model.frame
  resp_pos <- attr(terms(mf), "response")

  # Update the model frame with the data
  mf_update <- lapply(X = sim_response, FUN = function(x) {
    mf[,resp_pos] <- x
    return(mf)
  })

  # Update and fit the models with the new data
  sim_refit <- lapply(X = mf_update, function(x) update(object, data = x))

  # Calculate heritability on the original data
  base_herit <- herit_use(object)
  # Calculate the heritability for each new model object
  boot_herit <- lapply(sim_refit, herit_use)


  ## Calculate standard error, bias, and confidence intervals for heritability
  ## and the variance components

  # First pull out the heritability
  h_list <- sapply(X = boot_herit, "[[", "heritability")

  # Calculate the standard error
  h_se <- sd(h_list)
  # Calculate the bias
  h_bias <- mean(h_list) - base_herit$heritability

  # Calculate the confidence interval
  h_ci <- quantile(h_list, probs = c((alpha / 2), 1 - (alpha / 2)))


  ## Now do the same for the variance components
  var_comp_list <- lapply(X = boot_herit, "[[", "var_comp")
  # Condense into a usable matrix
  var_comp_mat <- sapply(X = var_comp_list, FUN = "[[", "variance")
  row.names(var_comp_mat) <- var_comp_list[[1]]$source

  ## Apply a function by rows
  var_comp_se <- apply(X = var_comp_mat, MARGIN = 1, FUN = sd)
  var_comp_bias <- rowMeans(var_comp_mat) - base_herit$var_comp$variance
  var_comp_ci <- apply(X = var_comp_mat, MARGIN = 1, FUN = quantile, probs = c((alpha / 2), 1 - (alpha / 2)))


  ## Package everything together
  h_summary <- data.frame(heritability = base_herit$heritability, std.error = h_se, bias = h_bias, t(h_ci),
                          check.names = FALSE)
  var_comp_summary <- data.frame(base_herit$var_comp, std.error = var_comp_se,
                                 bias = var_comp_bias, t(var_comp_ci), row.names = NULL,
                                 stringsAsFactors = FALSE, check.names = FALSE)

  # Return a list
  list(heritability = h_summary, var_comp = var_comp_summary)

} # Close function


#' @rdname herit_boot
#' @export
herit_boot.lmerMod <- function(object, exp, boot.reps = 1000, alpha = 0.05, ...) {

  # Extract other arguments
  other_args <- list(...)

  # Parse if necessary
  if (all(length(other_args) == 1, sapply(other_args, class) == "list")) {
    other_args <- do.call("c", other_args)
  }

  # Extract the model frame
  mf <- model.frame(object)

  # Create a function to calculate heritability from the model object
  herit_use <- function(object) herit(object = object, exp = exp, other_args)


  ## Perform bootstrapping
  # First simulate from the model
  sim_response <- simulate(object = object, nsim = boot.reps)

  # Get the column position of the response variable in the model.frame
  resp_pos <- attr(terms(mf), "response")

  # Update the model frame with the data
  mf_update <- lapply(X = sim_response, FUN = function(x) {
    mf[,resp_pos] <- x
    return(mf)
  })

  # Update and fit the models with the new data
  sim_refit <- lapply(X = mf_update, function(x) update(object, data = x))

  # Calculate heritability on the original data
  base_herit <- herit_use(object)
  # Calculate the heritability for each new model object
  boot_herit <- lapply(sim_refit, herit_use)


  ## Calculate standard error, bias, and confidence intervals for heritability
  ## and the variance components

  # First pull out the heritability
  h_list <- sapply(X = boot_herit, "[[", "heritability")

  # Calculate the standard error
  h_se <- sd(h_list)
  # Calculate the bias
  h_bias <- mean(h_list) - base_herit$heritability

  # Calculate the confidence interval
  h_ci <- quantile(h_list, probs = c((alpha / 2), 1 - (alpha / 2)))


  ## Now do the same for the variance components
  var_comp_list <- lapply(X = boot_herit, "[[", "var_comp")
  # Condense into a usable matrix
  var_comp_mat <- sapply(X = var_comp_list, FUN = "[[", "variance")
  row.names(var_comp_mat) <- var_comp_list[[1]]$source

  ## Apply a function by rows
  var_comp_se <- apply(X = var_comp_mat, MARGIN = 1, FUN = sd)
  var_comp_bias <- rowMeans(var_comp_mat) - base_herit$var_comp$variance
  var_comp_ci <- apply(X = var_comp_mat, MARGIN = 1, FUN = quantile, probs = c((alpha / 2), 1 - (alpha / 2)))


  ## Package everything together
  h_summary <- data.frame(heritability = base_herit$heritability, std.error = h_se, bias = h_bias, t(h_ci),
                          check.names = FALSE)
  var_comp_summary <- data.frame(base_herit$var_comp, std.error = var_comp_se,
                                 bias = var_comp_bias, t(var_comp_ci), row.names = NULL,
                                 stringsAsFactors = FALSE, check.names = FALSE)

  # Return a list
  list(heritability = h_summary, var_comp = var_comp_summary)

} # Close the function


