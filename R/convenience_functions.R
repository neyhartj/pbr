#' Convenience functions
#'
#' @description
#' Convenience functions. These are generally not to be called by the user
#'
avg_reps <- function(x) {

  # x must be of class 'table'
  stopifnot(is.table(x))

  R_1 <- sum(x)
  R_2 <- sum(x^2)

  numerator <- (R_1 - (R_2 / R_1))

  # Unique genotypes
  n <- sum(x >= 1)

  # Solve
  numerator / (n - 1)

} # Close function
