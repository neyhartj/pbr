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
#'
#' @importFrom utils combn
#' @importFrom stats as.dist
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
  env_names <- unique(as.character(x[[env.col]]))
  # Number of environments
  n_env <- length(env_names)

  ## Make a table of the observations in the dataset
  x_table <- table(x[c(env.col, gen.col)])
  # Are any greater than 1?
  if (any(x_table > 1)) stop ("Only one observation of each genotype in each environment should be included in the dataset.")

  # Order on environment, then genotype
  x1 <- x[order(x[[env.col]], x[[gen.col]]), , drop = FALSE]

  # Iterate over pairs
  D_ij <- combn(x = env_names, m = 2, FUN = function(env_pairs) {

    ## Subset the data for this environment pair
    env_index <- x1[[env.col]] %in% env_pairs
    x_sub <- x1[env_index, , drop = FALSE]

    # Find the common genotypes
    geno_count <- table(x_sub[[gen.col]])
    x_geno <- names(geno_count[geno_count == 2])

    # Subset the data again
    x_sub1 <- x_sub[x_sub[[gen.col]] %in% x_geno, , drop = FALSE]

    # Scale the phenotypic values by the mean and sd
    pheno_scale <- tapply(X = x_sub1[[pheno.col]], INDEX = x_sub1[[env.col]],
                          function(x) as.numeric(scale(x)), simplify = FALSE)

    # Find the squared differences between genotypes
    pheno_scale_diff <- tapply(X = unlist(pheno_scale, use.names = FALSE),
                               INDEX = x_sub1[[gen.col]],
                               function(x) diff(x)^2, simplify = FALSE)

    # Calculate the mean among these squared differences and return
    mean(unlist(pheno_scale_diff))

  })

  # Empty matrix
  D_mat <- matrix(0, nrow = n_env, ncol = n_env, dimnames = list(env_names, env_names))

  # Add data
  D_mat[lower.tri(D_mat)] <- D_ij

  # Output as distance matrix
  as.dist(D_mat)

} # Close the fuction
