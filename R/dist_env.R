#' Environmental distance matrix
#'
#' @description
#' Calculates the distance between environments based on the observed phenotypes
#' of genotypes shared between environments.
#'
#' @param x A \code{data.frame} of phenotypic observations on genotypes in environments. The data.frame
#' should be in tidy format, where every phenotypic observation has a separate row.
#' @param gen.col The column name that contains the names of genotypes.
#' @param env.col The column name that contains the names of environments.
#' @param pheno.col The column name that contains the phenotypic observations.
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
#' data("yang.barley")
#'
#' # Calculate distance
#' yang_barley_dist <- dist_env(x = yang.barley, env.col = "site")
#'
#' # Environment dendrogram
#' yang_barley_clust <- hclust(yang_barley_dist, method = "average")
#' plot(yang_barley_clust)
#' }
#'
#' @import dplyr
#'
#' @export
#'
#'
dist_env <- function(x, gen.col = "gen", env.col = "env", pheno.col = "yield") {

  # Verify that x is a data.frame
  stopifnot(is.data.frame(x))

  # Convert to data.frame
  x <- as.data.frame(x)

  # Verify columnnames exist
  if (!all(c(gen.col, env.col, pheno.col) %in% colnames(x)))
    stop("The values of gen.col or env.col or pheno.col are not columns in the x data.frame.")

  # Pull out environment names
  env <- x[,env.col] %>%
    unique() %>%
    as.character()

  # Number of environments
  n_env <- length(env)

  # Iterate over pairs
  D_ij <- combn(x = env, m = 2, FUN = function(env_pairs) {

    # Subset the data for the environments and common genotypes
    ## First set lazy evaluations of the filter criteria for environments
    env_filter_criteria <- lazyeval::interp(~ col %in% env_pairs, col = as.name(env.col))
    # Now for genotypes
    geno_filter_criteria <- lazyeval::interp(~ col %in% x_geno, col = as.name(gen.col))
    ## Do the same for the summary
    summarize_exp <- lazyeval::interp(~ mean(val), val = as.name(pheno.col))

    x_sub <- x %>%
      filter_(env_filter_criteria) %>%
      # Calculate the mean of each genotype in that environment
      group_by_(.dots = c(env.col, gen.col)) %>%
      summarize_(value = summarize_exp)

    # Find the genotypes that are common (i.e. 2 obs)
    x_geno <- x_sub %>%
      # Filter genotypes with less than 2 obs
      group_by_(.dots = gen.col) %>%
      filter(n() == 2) %>%
      distinct(gen.col) %>%
      pull()

    # In each environment, calculate the deviation from each observation to the mean of
    # all observation and divide by the sd of all observations (`scale` function)
    # Then find the mean squared differences. This is the distance
    x_sub %>%
      group_by_(.dots = env.col) %>%
      mutate(t_ij = scale(value)) %>%
      # Filter for genotypes that are found in both environments
      filter_(geno_filter_criteria) %>%
      group_by_(.dots = gen.col) %>%
      summarize(diff_sq = diff(t_ij)^2) %>%
      pull(diff_sq) %>%
      # Since the difference between t_ij is only for the common genotypes, the mean
      # at this point is division by the number of common genotypes
      mean()

  })

  # Empty matrix
  D_mat <- matrix(0, nrow = n_env, ncol = n_env, dimnames = list(env, env))

  # Add data
  D_mat[lower.tri(D_mat)] <- D_ij

  # Output as distance matrix
  as.dist(D_mat)

} # Close the fuction
