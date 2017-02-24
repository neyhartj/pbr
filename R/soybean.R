#' Yield data from a New York soybean trial
#'
#' @description
#' Total yield, in kilograms per hectare, of 7 soybean varieties grown in 10 environments between
#' 1977 and 1988.
#'
#' @format
#' A data frame with 280 observations and the following 4 variables:
#' \describe{
#'    \item{environment}{The growing environemnt, abbreviated at LYY, where L is location and Y is year. A factor with 10 levels}
#'    \item{variety}{Soyebean variety names. A factor with 7 levels.}
#'    \item{replicate}{The replicate number of the variety observation in the environment. A factor with 4 levels}
#'    \item{yield}{Yield in kilograms per hectare.}
#' }
#'
#' @source
#' Gauch, H.G. 1992. Statistical analysis of regional yield trials: AMMI analysis of factorial
#' designs. Elsevier, Amsterdam, the Netherlands.
#'
"soybean"
