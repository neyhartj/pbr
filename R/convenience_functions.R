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
#' @export
harm_mean <- function(x) length(x) / sum(1 / x)
