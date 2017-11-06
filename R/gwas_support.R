#' Calculate marker scores in association study
#'
#' @param M_test Marker genotype matrix
#' @param model The model name.
#' @param snp_info A data.frame of marker names, chromosomes, and postitions. This
#' is used for models "G" and "QG".
#' @param P3D Logical - are variance components known or should they be estimated per marker?
#' @param Hinv List of inverse of the phenotype vcov matrix (from the mixed-model output).
#' Must be of length 1 or the number of chromosomes.
#' @param y Vector of response values
#' @param X List of incidence matrices of fixed effects (must be full-rank, i.e. from svd).
#' Must of length 1 or the number of chromosomes.
#' @param Z0 Incidence matrix of random main effects (i.e. genotypes)
#' @param K0 List of covariance matrices of random main effects. Must be of length
#' 1 or the number of chromosomes.
#' @param X_fixed Incidence matrix of environment fixed effects (but not including intercept)
#' @param Z_rand Incidence matrix of genotype random effects (but not including intercept)
#'
#' @details
#' Calculate the p-value for each marker in an association study. The main effect
#' of each marker is tested using a Wald test.
#'
#' @importFrom EMMREML emmreml
#' @importFrom purrr map
#' @import dplyr
#'
#'
score_calc <- function(M_test, model, snp_info, P3D, Hinv, y, X, Z0, K0, X_fixed, Z_rand) {

  # Create a list of markers by chromosome or just a list
  # of markers
  mar_list <- snp_info %>%
    split(.[,2]) %>%
    map(1)

  # Subset from the X and K matrix those chromosomes in the mar_list, if possible
  if (model %in% c("G", "QG")) {
    K_use <- K0[names(mar_list)]
    Hinv_use <- Hinv[names(mar_list)]
    X_use <- ifelse(model == "QG", X[names(mar_list)], X)

  } else {
    X_use <- X
    K_use <- K0
    Hinv_use <- Hinv

  }

  Z_use <- Z0

  # number of observations
  n <- length(y)

  # Map over the list of markers and the X or K matrix
  scores <- pmap(list(mar_list, X_use, K_use, Hinv_use), .f = function(mar, x, k, h) {

    # Subset the marker matrix for those markers
    m <- M_test[,mar,drop = FALSE]

    # Apply function over the column of the M matrix (i.e. markers)
    apply(X = m, MARGIN = 2, FUN = function(snp) {

      ## First calculate the main effect of the marker
      # Model matrix of SNP main effect
      X_snp_main <- Z_rand %*% snp
      X_use1 <- cbind(x, X_snp_main)

      ## Re-estimate variance components?
      if (!P3D) {
        # Fit the model and get the inverse of the phenotypic vcov matrix
        fit <- emmreml(y = y, X = X_use1, Z = Z0, K = k)
        H2inv <- H_inv(Vu = fit$Vu, Ve = fit$Ve, n = n, Z = Z0, K = k)

      } else {
        H2inv <- h

      }

      # W matrix and inverse
      W <- crossprod(X_use1, H2inv %*% X_use1)
      Winv <- try(solve(W), silent = TRUE)

      # If the inverse is successful, calculate the p-value for that SNP using
      # an F test
      if (class(Winv) != "try-error") {
        # Number of fixed terms
        p <- ncol(X_use1)
        v2 <- n - p
        # Location of marker fixed terms
        p_mar <- seq(ncol(x) + 1, p)

        # Vector of fixed effect coefficients
        beta <- Winv %*% crossprod(X_use1, H2inv %*% y)
        # Residuals
        resid <- y - X_use1 %*% beta
        # Estimate of sum of squared residuals
        s2 <- as.double(crossprod(resid, H2inv %*% resid))/v2
        # VCOV of fixed effect coefficients
        CovBeta <- s2 * Winv

        beta0 <- beta[p_mar,, drop = FALSE]
        # The Wald-test statistic is the square of the coefficient divided by the variance of that coefficient
        w0 <- beta0^2 / diag(CovBeta)[p_mar]
        # Under the NULL of beta = 0, the w statistic is chi-square distributed with 1 df
        p_main <- pchisq(q = w0, df = 1, lower.tail = FALSE)


        colnames(beta0) <- "beta0_hat"

        # Output data.frame
        out_df <- data.frame(term = "main_effect", df = 1, W_statistic = w0, p_value = p_main,
                             row.names = NULL, stringsAsFactors = FALSE)

        # Output list
        return(list(test = out_df, beta0_hat = beta0, beta1_hat = NA))

      } else {

        out_df <- data.frame(term = "main_effect", df = 1, W_statistic = NA, p_value = NA,
                             row.names = NULL, stringsAsFactors = FALSE)

        return(list(test = out_df, beta0_hat = NA))

      } }) }) # Close the apply and pmap functions


  # Combine the data.frames and return
  map_df(scores, ~mutate(map_df(., "test"), marker = names(.))) %>%
    mutate(estimate = unlist(map(scores, ~map_dbl(., "beta0_hat")))) %>%
    select(marker, term, df, estimate, names(.))

} # CLose the function


#' Plot the output of an association analysis
#'
#' @param x An object of class \code{gwas} or a list of such objects.
#' @param fdr.level The false discovery rate cutoff for declaring adjusted p-values
#' significant.
#' @param type The type of plot to create. Can be \code{"manhattan"} or
#' \code{"qq"}.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
plot_gwas <- function(x, fdr.level = 0.05, type = c("manhattan", "qq")) {

  ## Is x a list of lists or a list?
  x_class <- class(x)
  # If gwas, convert to list, or keep the same
  if (x_class == "gwas") {
    x1 <- list(x)
  } else{
    x1 <- x
  }

  # Match the plot type
  type <- match.arg(type)

  # Iterate over the list of gwas outputs
  plot_data <- map_df(x1, ~mutate(.$scores, model = .$metadata$model))

  # Extract the column names for chromosome and position
  chrom_name <- colnames(plot_data)[2]
  pos_name <- colnames(plot_data)[3]

  # Adjust p-values using the qvalue function
  plot_data_adj <- plot_data %>%
    mutate_at(vars(chrom_name), as.factor) %>%
    group_by(model, trait) %>%
    mutate(p_value_adj = p.adjust(p_value, "fdr"),
           neg_log_p_adj = -log10(p_value_adj),
           neg_log_fdr = -log10(fdr.level)) %>%
    ungroup()

  # Extract p_values and estimate the expected p_value
  p_values <- plot_data_adj %>%
    select(model, trait, p_value) %>%
    group_by(model, trait) %>%
    arrange(p_value) %>%
    mutate(exp_p_value = ppoints(n = n())) %>%
    mutate_at(vars(contains("p")), ~-log10(.))


  # Plot
  if (type == "manhattan") {
    g <- plot_data_adj %>%
      ggplot(aes(x = eval(as.name(pos_name)), y = neg_log_p_adj)) +
      geom_point(aes(col = eval(as.name(chrom_name)))) +
      geom_hline(aes(yintercept = neg_log_fdr), lty = 2, lwd = 1) +
      facet_grid(trait + model ~ chrom, scales = "free", switch = "x") +
      ylab("-log10(q)") +
      xlab("Position") +
      scale_color_discrete(guide = FALSE) +
      theme_bw()

  } else {

    g <- p_values %>%
      ggplot(aes(x = exp_p_value, y = p_value, col = model)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1) +
      facet_wrap(~ trait, ncol = 2) +
      ylab("Observed -log10(p)") +
      xlab("Expected -log10(p)") +
      theme_bw()

  }

  # Return the graph
  invisible(print(g))

} # CLose the function



#' Calculate the inverse of the H matrix
#'
#' @param Vu Estimated variance of random effects
#' @param Ve Estimated variance of residuals
#' @param n Number of observations
#' @param Z Random effects incidence matrix
#' @param K Covariance matrix for random effects
#'
H_inv <- function(Vu, Ve, n, Z, K) {

  H <- (Z %*% K %*% t(Z)) + ((Ve / Vu) * diag(n))
  solve(H)

} # Close the function

