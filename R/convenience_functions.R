#' Convenience functions
#'
#' @description
#' Convenience functions. These are generally not to be called by the user
#'
#' @details
#'
#' The formula for the harmonic mean given in the function \code{harm_mean} is:
#'
#' \eqn{e_h = n / \sum{\frac{1}{e_j}}}
#'
#' By default, the function will remove Infinite values and calculate the sum
#' using non-missing values.
#'
#' @export
#'
harm_mean <- function(x) {
  inv_x <- 1 / x
  inv_x[is.infinite(inv_x)] <- NA
  1 / mean(inv_x, na.rm = TRUE)
}
