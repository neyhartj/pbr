#' Compare two models with a likelihood ratio test
#'
#' @description
#' Performs a likelihood ratio test to compare two models
#'
#' @param model1 A fitted model object.
#' @param model2 See \code{model1}.
#'
#' @export
#'
#'
lr_test <- function(model1, model2) {

  ## Are both models of the same class?
  stopifnot(class(model1) == class(model2))

  # Remove classes
  classes_keep <- c("lm", "aov", "merMod", "lmerMod")
  # Stop if not classes
  model1Classes <- class(model1) %in% classes_keep
  model2Classes <- class(model2) %in% classes_keep

  if (!all(model1Classes, model2Classes))
    stop("Models are not one of classes: ", paste(classes_keep, collapse = ", "))

  model_list <- list(model1 = model1, model2 = model2)

  ## Get the number of terms from each model
  terms_list <- lapply(model_list, terms)
  terms_list1 <- lapply(terms_list, attr, "term.labels")
  n_terms <- sapply(X = terms_list1, length)

  # List of df
  df_list <- sapply(model_list, df.residual)

  # Degrees of freedom
  full_model <- names(which.min(df_list))
  red_model <- names(which.max(df_list))

  ## Get the log-likelihoods and store as list
  ll_list <- sapply(X = model_list, FUN = logLik)

  # Calculate the likelihood ratio
  lr <- -2 * (ll_list[red_model] - ll_list[full_model])
  # Calculate pvalue
  df <- df_list[red_model] - df_list[full_model]
  p_value <- pchisq(q = lr, df = df, lower.tail = FALSE)

  # Export a data.frame
  data.frame(
    full_model = full_model,
    df = df,
    statistic = lr,
    p_value = p_value,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

}
