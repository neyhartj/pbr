#' Pull main and interaction effects from a model
#'
#' @param object A fitted model object.
#' @param gen.col The name of the column in \code{data} containing the genotype factor.
#' @param env.col The name of the column in \code{data} containing the environment factor.
#'
#' @importFrom effects Effect allEffects
#'
#' @export
#'
pull_effects <- function(object, gen.col = "gen", env.col = "env") {

  UseMethod(generic = "pull_effects", object = object)

}

#' @rdname pull_effects
pull_effects.lm <- function(object, gen.col = "gen", env.col = "env") {

  # Get the coefficients
  fit_coef <- coef(object)
  mf <- model.frame(object)

  ## Calculate effects and means
  # Grand mean
  grand_mean <- fit_coef[1]

  all_terms <- attr(terms(mf), "term.labels")
  ## Designate the terms as main or interacting
  main_terms <- all_terms[!grepl(pattern = ":", x = all_terms)]

  ## Stop if one of the terms is not gen.col:env.col
  gxe_term <- paste(gen.col, env.col, sep = ":")
  if (!gxe_term %in% all_terms) {
    stop("The GxE term is not properly specified in the model. It must take the
         form 'gen.col:env.col'")
  }

  ## Get the main effects
  effs <- lapply(X = main_terms, FUN = Effect, mod = object)
  effs1 <- lapply(X = effs, FUN = "[[", "fit")
  effs2 <- setNames(object = lapply(X = effs1, FUN = `-`, grand_mean), nm = c(gen.col, env.col))

  # Convert to data.frames
  gen_eff <- data.frame(gen = levels(mf[[gen.col]]), gen_main_eff = c(effs2[[gen.col]]),
                        row.names = NULL, stringsAsFactors = FALSE)
  env_eff <- data.frame(env = levels(mf[[env.col]]), env_main_eff = c(effs2[[env.col]]),
                        row.names = NULL, stringsAsFactors = FALSE)

  ## Get the interaction effects (gxe)
  int_eff <- allEffects(object)[[gxe_term]]
  int_eff1 <- cbind(int_eff$x, int_eff = int_eff$fit)

  ## Combine
  all_effects <- merge(x = merge(x = int_eff1, y = gen_eff, by = gen.col), y = env_eff, by = env.col)
  # Recalculate interaction effects
  all_effects$int_eff <- all_effects$int_eff - grand_mean - all_effects$gen_main_eff - all_effects$env_main_eff
  # Add grand mean
  all_effects$grand_mean <- grand_mean

  # Return
  return(all_effects)

}




#' The Sites Regression (SREG, or GGE) model
#'
#' @description Fit a sites-regression model (SREG, or GGE) to analyze genotype-
#' environment interactions in plant breeding experiments.
#'
#' @param formula A formula specifying the model.
#' @param data A data frame in which the variables specified in the formula will be found.
#' If missing, the variables are searched for in the standard way.
#' @param gen.col The name of the column in \code{data} containing the genotype factor.
#' @param env.col The name of the column in \code{data} containing the environment factor.
#' @param n.terms The number of multiplicative terms to include in the model.
#' Defaults to 1 (i.e. SREG_1).
#'
#' @examples
#' # Use the gauch.soy dataset
#' data("gauch.soy")
#'
#' # Filter
#' gauch_soy1 <- droplevels(subset(gauch.soy, env %in% c("A77", "A80", "A83")))
#'
#' sreg_fit <- sreg(formula = yield ~ gen + env + gen:env, data = gauch_soy1)
#'
#' @importFrom tidyr spread_
#'
#' @export
#'
sreg <- function(formula, data, gen.col = "gen", env.col = "env", n.terms = 1) {

  stopifnot(inherits(formula, "formula"))

  # Make sure the env.col is among the names in data
  stopifnot(gen.col %in% names(data))
  stopifnot(env.col %in% names(data))

  # Edit the formula environment
  environment(formula) <- environment()

  ## Create a model.frame
  mf <- model.frame(formula = formula, data = data, drop.unused.levels = TRUE)

  # Get the model terms
  mod_terms <- terms(mf)
  pred_vars <- attr(mod_terms, "predvars")
  all_terms <- attr(mod_terms, "term.labels")

  ## Stop if one of the terms is not gen.col:env.col
  gxe_term <- paste(gen.col, env.col, sep = ":")
  if (!gxe_term %in% all_terms) {
    stop("The GxE term is not properly specified in the model. It must take the
         form 'gen.col:env.col'")
  }


  ## Vector of predictor variables
  ## The first two names (for a single response) are "list" and the response variable
  var_names <- as.character(pred_vars)[-1]
  resp_var_names <- var_names[attr(mod_terms, "response")]
  pred_var_names <- setdiff(var_names, resp_var_names)

  ## Fit the fixed-effect model
  # First set the contrasts
  contr <- as.list(rep("contr.sum", length(pred_var_names)))
  names(contr) <- pred_var_names

  ## Set the contrasts
  for (var in pred_var_names) {
    C_mat <- C(object = mf[[var]], contr = "contr.sum")
    # Set the column names of the contrast matrix
    contrasts(C_mat) <- `colnames<-`(x = contrasts(C_mat), value = head(levels(C_mat), -1))

    mf[[var]] <- C_mat
  }


  fit <- lm(formula = formula, data = mf)

  ## Get the effects
  all_effects <- pull_effects(object = fit, gen.col = gen.col, env.col = env.col)



  ##############
  ##############

  ## Create a two-way table of genotype + gxe effects
  gge_effects <- all_effects
  gge_effects$gge <- gge_effects$gen_main_eff + gge_effects$int_eff

  gge_mat <- spread_(data = gge_effects[,c(gen.col, env.col, "gge")], key_col = env.col, value_col = "gge")
  row.names(gge_mat) <- gge_mat[[1]]
  gge_tab <- as.matrix(gge_mat[,-1])
  dimnames(gge_tab) <- setNames(dimnames(gge_tab), c(gen.col, env.col))


  ### Singular value decomp
  gge_svd <- svd(x = gge_tab)
  # Conver the singular values to a matrix
  gge_svd$d <- matrix(gge_svd$d, nrow = 1)

  ## Split in to a list of length(multi_svd$d) - one for each component
  n_comp <- ncol(gge_svd$d)
  gge_svd_list <- lapply(X = seq(n_comp), FUN = function(n) lapply(X = gge_svd, "[", ,n, drop = FALSE))

  ## Pull out the desired components
  comp_tokeep <- gge_svd_list[seq(n.terms)]
  comp_residual <- gge_svd_list[setdiff(seq(n_comp), seq(n.terms))]

  # If n.term is 1, the eigenvectors are the primary effects
  if (n.terms == 1) {
    gen_prim_eff <- setNames(object = c(comp_tokeep[[1]]$u), nm = paste(rownames(gge_tab), sep = ""))
    env_prim_eff <- setNames(object = c(comp_tokeep[[1]]$v), nm = paste(colnames(gge_tab), sep = ""))

  } else {
    gen_prim_eff <- env_prim_eff <- NULL

  }


  ## For each component, generate PC scores by multplying the eigenvectors by the
  ## square root of the singular value
  pc_score_list_tokeep <- lapply(X = comp_tokeep, FUN = function(comp) {
    gen_score <- comp$u * sqrt(c(comp$d))
    env_score <- comp$v * sqrt(c(comp$d))
    # Multiply them
    as.vector(gen_score) %*% t(as.vector(env_score))
  })

  ## Do the same for the residual components
  pc_score_list_residual <- lapply(X = comp_residual, FUN = function(comp) {
    gen_score <- comp$u * sqrt(c(comp$d))
    env_score <- comp$v * sqrt(c(comp$d))
    # Multiply them
    as.vector(gen_score) %*% t(as.vector(env_score))
  })


  ## Sum the desired PC scores
  theta <- Reduce(f = `+`, x = pc_score_list_tokeep)
  # Add dimnames
  dimnames(theta) <- dimnames(gge_tab)

  ## Sum the residual scores
  rho <- Reduce(f = `+`, x = pc_score_list_residual)
  # Add dimnames
  dimnames(rho) <- dimnames(gge_tab)


  ## Create a data.frame with the theta and rho values
  theta_df <- as.data.frame(as.table(theta), responseName = "theta")
  rho_df <- as.data.frame(as.table(rho), responseName = "rho")
  # Combine
  gge_svd_df <- cbind(theta_df, rho = rho_df$rho)

  ## Combine with all_effects
  all_effects1 <- merge(x = all_effects, y = gge_svd_df, by = c(gen.col, env.col))
  ## Add the response values
  mf_effects <- merge(x = data, y = all_effects1)


  ## The fitted values are the main effects + the svd effects
  fitted_values <- mf_effects$grand_mean + mf_effects$env_main_eff + mf_effects$theta
  # The residuals are the response values - the fitted values
  residual_values <- mf_effects[[resp_var_names]] - fitted_values

  ## Subset the mf_effects df for unique G-E combinations
  mf_effects1 <- unique(mf_effects[,c(gen.col, env.col, "grand_mean", "gen_main_eff", "env_main_eff",
                                      "int_eff", "theta")])

  ## Predicted values
  mf_effects1$ybar <- mf_effects1$grand_mean + mf_effects1$gen_main_eff + mf_effects1$env_main_eff +
    mf_effects1$int_eff
  mf_effects1$fitted <- mf_effects1$grand_mean + mf_effects1$env_main_eff + mf_effects1$theta
  mf_effects1$residual <- mf_effects1$ybar - mf_effects1$fitted

  ## Add the primary genotype and environment effects
  mf_effects2 <- merge(mf_effects1, as.data.frame(env_prim_eff), by.x = env.col, by.y = "row.names")
  mf_effects3 <- merge(mf_effects2, as.data.frame(gen_prim_eff), by.x = gen.col, by.y = "row.names")

  ## Assemble an output list
  sreg <- list(
    fit = fit,
    coefficients = coef(fit),
    residuals = residual_values,
    fitted.values = fitted_values,
    primary.effect = setNames(list(gen_prim_eff, env_prim_eff), c(gen.col, env.col)),
    effects = mf_effects3,
    n.terms = n.terms
  )

  # Add class and return
  structure(sreg, class = "sreg")

}
