#' Perform an LSD test among means
#'
#' @description
#' Performs an least-significant difference test among group means, for example
#' among genotype means.
#'
#' @param df A \code{data.frame} with two columns. The first column includes the names
#' of elements in the group, and the second column is the mean of that element
#' in the group.
#' @param LSD The estimated least-significant difference (i.e. determined by the
#' \code{\link{pbr}[LSD]} function.
#' @param alpha The significance level for hypothesis testing.
#' @param p.adj The p-value adjustment method for multiple comparison testing.
#'
#'
#' @export
#'
LSD_test <- function(df, LSD, alpha = 0.05,
                     p.adj = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {

  ## Error handling
  df <- as.data.frame(df)

  # The df must have two columns
  if (ncol(df) > 2)
    stop("The 'df' input must have two columns: the group name and the mean of each group.")

  # The first column must be a character and the second a numeric
  df[,1] <- as.character(df[,1])
  df[,2] <- as.numeric(df[,2])

  # Alpha must be between 0 and 1 (not inclusive)
  if (!(alpha > 0 & alpha < 1))
    stop("The argument 'alpha' must be greater than 0 and less than 1.")

  # Form pairwise comparisons of means
  groups <- df[,1]

  group_pairs <- combn(x = groups, m = 2)

  group_means <- apply(X = group_pairs, MARGIN = 2, FUN = function(grp)
    diff(subset(x = df, df[,1] %in% as.character(grp), 2, drop = TRUE)))

  # Bind
  group_comparisons <- data.frame(data.frame(t(group_pairs), stringsAsFactors = FALSE),
                                  group_means,
                                  row.names = apply(X = group_pairs, MARGIN = 2, FUN = paste, collapse = "-")) %>%
    structure(names = c(paste("group", seq(2), sep = ""), "diff"))

  ## Use the LSD to determine whether the differences are significant
  group_comparisons$sig <- abs(group_comparisons$diff) >= LSD

  # Return the test
  return(group_comparisons)

} # Close the function

