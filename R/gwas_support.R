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
#' @param Z1 Incidence matrix of gxe random interaction effects.
#' @param K1 List of covariance matrices of random interaction effects. Must be of length
#' 1 or the number of chromosomes.
#' @param X_fixed Incidence matrix of environment fixed effects (but not including intercept)
#' @param Z_rand Incidence matrix of genotype random effects (but not including intercept)
#'
#' @details
#' Calculate the p-value for each marker in an association study. The main effect
#' of each marker and QTLxE interaction is tested using a Wald test.
#'
#' @importFrom sommer mmer
#' @importFrom EMMREML emmreml
#' @importFrom purrr map
#' @import dplyr
#'
#'
score_calc <- function(M_test, model, snp_info, P3D, Hinv, y, X, Z0, K0, Z1 = NULL,
                       K1 = NULL, X_fixed, Z_rand) {

  # Create a list of markers by chromosome or just a list
  # of markers
  mar_list <- snp_info %>%
    split(.[,2]) %>%
    map(1)

  # Subset from the X and K matrix those chromosomes in the mar_list, if possible
  if (model %in% c("G", "QG")) {
    X_use <- X[names(mar_list)]
    K_use <- K0[names(mar_list)]
    Hinv_use <- Hinv[names(mar_list)]

  } else {
    X_use <- X
    K_use <- K0
    Hinv_use <- Hinv

  }

  Z_use <- Z0

  # Map over the list of markers and the X or K matrix
  scores <- pmap(list(mar_list, X_use, K_use, Hinv_use), .f = function(mar, x, k, h) {

    # Subset the marker matrix for those markers
    m <- M_test[,mar]

    # Apply function over the column of the M matrix (i.e. markers)
    apply(X = m, MARGIN = 2, FUN = function(snp) {

      # number of observations
      n <- length(y)

      ## First calculate the main effect of the marker
      # Model matrix of SNP main effect
      X_snp_main <- Z_rand %*% snp
      X_use1 <- cbind(x, X_snp_main)

      ## Re-estimate variance components?
      if (!P3D) {
        # Fit the model
        # The stream should split between simple and Q - and other models
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

