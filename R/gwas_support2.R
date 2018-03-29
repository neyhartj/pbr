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
#' @param Z_rand Incidence matrix of random effects (but not including intercept)
#' @param Z1 Incidence matrix of gxe random effects
#' @param K1 List of covariance matrices of gxe random effects.  Must be of length
#' 1 or the number of chromosomes.
#' @param X_fixed Incidence matrix of environment fixed effects (but not including intercept)

#'
#' @details
#' Calculate the p-value for each marker in an association study. The main effect
#' of each marker and the qxe effect is tested using a Wald test.
#'
#' @importFrom EMMREML emmreml
#' @importFrom purrr map map_df
#' @import dplyr
#' @import stringr
#'
#'
score_calc <- function(M_test, model, snp_info, P3D, Hinv, y, X, Z0, K0, Z_rand,
                       X_fixed = NULL) {

  # Create a list of markers by chromosome or just a list
  # of markers
  mar_list <- lapply(X = split(snp_info, snp_info$chrom), FUN = "[[", "marker")

  # If H_inv is NULL, set to NA
  if (is.null(Hinv)) Hinv <- NA

  # Then split by the presence of G
  if (str_detect(model, "G")) {
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

  ## Test markers
  # Map over the list of markers and the X or K matrix
  scores <- list(mar_list, X_use, K_use, Hinv_use) %>%
    pmap(function(mar, x, k, h) {

    # Subset the marker matrix for those markers
    m <- M_test[,mar,drop = FALSE]
    # Designate the inverse phenotype vcov matrix
    H2inv <- h

    ## Create an incidence matrix of SNP genotypes per phenotypic observation
    # Model matrix of SNP main effect
    X_snp_main <- Z_rand %*% m

    # Apply function over the column of the M matrix (i.e. markers)
    apply(X = X_snp_main, MARGIN = 2, FUN = function(snp_main) {

      ## First calculate the main effect of the marker
      X_use1 <- cbind(x, snp_main)

      # Number of fixed terms
      p <- ncol(X_use1)

      ## Re-estimate variance components?
      if (!P3D) {
        # Fit the model
        fit <- emmreml(y = y, X = X_use1, Z = Z0, K = k, varbetahat = TRUE)

        # Get the beta coefficient for the marker effect
        beta <- c(fit$betahat)
        # Get the variance of the coefficients
        varBeta <- fit$varbetahat

      } else {

        # W matrix and inverse
        W <- crossprod(X_use1, H2inv %*% X_use1)
        Winv <- try(solve(W), silent = TRUE)

        # If the inverse is successful, calculate the p-value for that SNP using
        # an F test
        if (class(Winv) != "try-error") {
          # Vector of fixed effect coefficients
          beta <- c(Winv %*% crossprod(X_use1, H2inv %*% y))
          # Predictions (ignore random effects)
          y_hat <- X_use1 %*% beta
          # Residuals
          resid <- y - y_hat

          # Number of fixed terms
          v2 <- n - p
          # Estimate of sum of squared residuals
          s2 <- as.numeric(crossprod(resid, H2inv %*% resid))/v2
          # VCOV of fixed effect coefficients
          CovBeta <- s2 * Winv
          varBeta <- diag(CovBeta)

        } else {
          # Just return NA
          beta <- varBeta <- NA

        }

      }

      # Test the marker effect
      # Location of marker effects
      mar_p <- seq(ncol(x) + 1, p)
      # Location of marker main effect terms
      p_mar <- mar_p[1]
      p_mar1 <- setdiff(mar_p, p_mar)

      # Main effect test
      beta0 <- beta[p_mar]
      se0 <- sqrt(varBeta[p_mar])
      # The Wald-test statistic is the square of the coefficient divided by the variance of that coefficient
      w0 <- (beta0^2) / (se0^2)
      # Under the NULL of beta = 0, the w statistic is chi-square distributed with 1 df
      p_main <- pchisq(q = w0, df = 1, lower.tail = FALSE)

      # Output data.frame
      data.frame(beta = beta0, se = se0, df = 1, statistic = w0, pvalue = p_main,
                 row.names = NULL, stringsAsFactors = FALSE)

    })

    }) # Close the apply and pmap functions

  # List of markers that were tested (and information)
  scores_df <- bind_cols(select(snp_info, marker:pos), map_df(scores, ~bind_rows(.)))

} # Close the function


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

  # Number of chromosomes
  n_chrom <- n_distinct(plot_data$chrom)
  # Color palette for chromosomes
  chrom_color <- rep(c("grey", "black"), length.out = n_chrom) %>%
    set_names(unique(plot_data$chrom))

  # Adjust p-values using the qvalue function
  plot_data_adj <- plot_data %>%
    mutate_at(vars(chrom_name), as.factor) %>%
    group_by(model, trait, term) %>%
    mutate(p_value_adj = p.adjust(p_value, "fdr"),
           neg_log_p_adj = -log10(p_value_adj),
           neg_log_fdr = -log10(fdr.level)) %>%
    ungroup() %>%
    # assign chromosome color
    mutate(chrom_color = str_replace_all(chrom, chrom_color)) %>%
    # Remove NAs
    filter(!is.na(p_value))

  # Extract p_values and estimate the expected p_value
  p_values <- plot_data_adj %>%
    select(model, trait, term, p_value) %>%
    group_by(model, term, trait) %>%
    arrange(p_value) %>%
    mutate(exp_p_value = ppoints(n = n())) %>%
    mutate_at(vars(contains("p")), ~-log10(.))


  # Plot
  if (type == "manhattan") {
    g <- plot_data_adj %>%
      ggplot(aes(x = eval(as.name(pos_name)), y = neg_log_p_adj)) +
      geom_point(col = plot_data_adj$chrom_color) +
      geom_hline(aes(yintercept = neg_log_fdr), lty = 2, lwd = 1) +
      facet_grid(trait + model + term ~ chrom, scales = "free", switch = "x") +
      ylab("-log10(p)") +
      xlab("Position") +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.spacing.x = unit(0, units = "in")
      )

  } else {

    g <- p_values %>%
      ggplot(aes(x = exp_p_value, y = p_value, col = model)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1) +
      facet_wrap(term ~ trait, ncol = 2) +
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
#' @param weights Weights for the random effect(s). If the number of random
#' effects is 1, weight = 1.
#' @param Z List of random effects incidence matrices
#' @param K List of covariance matrices for random effects
#'
#' @importFrom purrr pmap
#'
H_inv <- function(Vu, Ve, weights, Zlist, Klist) {

  # Number of observations
  n <- nrow(Zlist[[1]])
  # R matrix
  R <- diag(n)
  # Lambda
  lambda <- Ve / Vu

  # Create an empty list
  ZKZlist <- vector("list",  length(Zlist))

  # Iterate
  for (i in seq_along(ZKZlist)) {
    ZKZlist[[i]] <- (Zlist[[i]] %*% Klist[[i]] %*% t(Zlist[[i]])) * weights[i]
  }

  # Sum them

  ZKZt <- purrr::reduce(ZKZlist, `+`)

  # Add the residual covariance
  H <- ZKZt + (lambda * R)
  as.matrix(solve(H))

} # Close the function
