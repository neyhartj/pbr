#' Environmental distance matrix
#'
#' @description
#' Calculates the distance between environments based on the observed phenotypes
#' of genotypes shared between environments.
#'
#' @param x A \code{matrix} of phenotypic values of genotypes in each environment.
#' The matrix should be of dimensions \code{i} x \code{j} where \code{i} is the
#' number of genotypes and \code{j} is the number of environments. The matrix
#' column names should be environment names.
#'
#' @details
#' The distance between environments (\eqn{D_jj}) cannot be less than \code{0}
#' and cannot exceed \code{4}. A distance close to \code{0} indicates that the relative
#' performance of genotypes is similar (i.e. the correlation between genotype
#' performance across environments is \eqn{r_jj = 1}) while a distance
#' close to \code{2} indicates that the correlation between genotype performance
#' is close to \eqn{r_jj = 0}. Finally, as \eqn{D_jj} approaches 4, crossover
#' interactions are abundant and the correlation of genotype performance is close
#' to \eqn{r_jj = -1}.
#'
#' @return
#' A \code{dist} matrix of distances between environments. This \code{dist} object
#' can be used by clustering functions such as \code{\link[stats]{hclust}}.
#'
#' @examples
#'
#' \dontrun{
#' data("barley")
#'
#' # Convert to useable matrix
#' barley_mat <- barley %>%
#'   mutate(env = paste(site, year, sep = "_")) %>%
#'   select(variety, env, yield) %>%
#'   spread(env, yield) %>%
#'   select(-1) %>%
#'   as.matrix()
#'
#' # Calculate distance
#' barley_dist <- dist_env(barley_mat)
#'
#' # Environment dendrogram
#' barley_clust <- hclust(barley_dist, method = "average")
#' plot(barley_clust)
#' }
#'
#' @import dplyr
#'
#' @export
#'
#'
dist_env <- function(x) {

  # Verify that x is a matrix
  if (!is.matrix(x)) stop("x must be a matrix.")

  # Verify columnnames exist
  if (is.null(colnames(x))) stop("x must have column names (i.e. environment names).")

  # Number of environments
  n.env <- ncol(x)

  # Find all pairwise combinations
  env.pairwise <- expand.grid(colnames(x), colnames(x))

  # Apply a function over these pairs
  D_jj <- apply(X = env.pairwise, MARGIN = 1, FUN = function(env.pairs) {

    # Subset the data
    x[,as.character(as.matrix(env.pairs))] %>%

      # In each environment, calculate the deviation from each observation to the mean of
      # all observation and divide by the sd of all observations
      # Then find the mean squared differences. This is the distance
      apply(MARGIN = 2, FUN = function(p) (p - mean(p, na.rm = T)) / sd(p, na.rm = T)) %>%
      # Subset the complete cases
      na.omit() %>%
      # Square the difference
      apply(MARGIN = 1, FUN = function(pairs) diff(pairs)^2) %>%
      mean()

  })

  # Create matrix
  env_mat <- matrix(data = D_jj, nrow = n.env, ncol = n.env, dimnames = list(colnames(x), colnames(x)))

  # Output as distance matrix
  as.dist(env_mat)

}
