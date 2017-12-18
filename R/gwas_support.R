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
#' @importFrom sommer mmer
#' @importFrom purrr map map_df
#' @import dplyr
#' @import stringr
#'
#'
score_calc <- function(M_test, model, snp_info, P3D, Hinv, test_qxe = FALSE,
                       y, X, Z0, K0, Z_rand, Z1 = NULL, K1 = NULL, X_fixed = NULL) {

  # Create a list of markers by chromosome or just a list
  # of markers
  mar_list <- snp_info %>%
    split(.$chrom) %>%
    map("marker")

  # Subset from the X and K matrix those chromosomes in the mar_list, if possible
  # First split stream by the number of variance components
  if (model %in% c("KKE", "GGE", "QKKE", "QGGE")) {

    # Split by presence of G
    if (str_detect(model, "G")) {
      K_use <- K0[names(mar_list)]
      Hinv_use <- Hinv[names(mar_list)]
      X_use <- ifelse(model == "QGGE", X[names(mar_list)], X)
      K1_use <- K1[names(mar_list)]

    } else {
      X_use <- X
      K_use <- K0
      Hinv_use <- Hinv
      K1_use <- K1

    }

    Z1_use <- Z1

  } else {
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
  }

  Z_use <- Z0

  # number of observations
  n <- length(y)

  ## Test markers
  # Split the stream by the number of variance components
  if (!model %in% c("KKE", "GGE", "QKKE", "QGGE")) {

    # Map over the list of markers and the X or K matrix
    scores <- pmap(list(mar_list, X_use, K_use, Hinv_use), .f = function(mar, x, k, h) {

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
        # Should QxE be tested?
        if (!test_qxe) {
          X_use1 <- cbind(x, snp_main)

        } else {
          # Model matrix of SNP x Environment
          X_snp_qxe <- c(snp_main) * X_fixed
          # Include main effect of SNP, then remove the last QTLxE term to keep model full rank
          X_use1 <- cbind(x, snp_main, X_snp_qxe[,-ncol(X_snp_qxe)])

        }

        ## Re-estimate variance components?
        if (!P3D) {
          # Fit the model and get the inverse of the phenotypic vcov matrix
          fit <- emmreml(y = y, X = X_use1, Z = Z0, K = k)
          H2inv <- H_inv(Vu = fit$Vu, Ve = fit$Ve, n = n, Z = Z0, K = k)

        }

        # W matrix and inverse
        W <- crossprod(X_use1, H2inv %*% X_use1)
        Winv <- try(solve(W), silent = TRUE)

        # If the inverse is successful, calculate the p-value for that SNP using
        # an F test
        if (class(Winv) != "try-error") {
          # Vector of fixed effect coefficients
          beta <- Winv %*% crossprod(X_use1, H2inv %*% y)
          # Predictions (ignore random effects)
          y_hat <- X_use1 %*% beta
          # Residuals
          resid <- y - y_hat

          # Number of fixed terms
          p <- ncol(X_use1)
          v2 <- n - p
          # Estimate of sum of squared residuals
          s2 <- as.numeric(crossprod(resid, H2inv %*% resid))/v2
          # VCOV of fixed effect coefficients
          CovBeta <- s2 * Winv

          # Location of marker effects
          mar_p <- seq(ncol(x) + 1, p)
          # Location of marker main effect terms
          p_mar <- mar_p[1]
          p_mar1 <- setdiff(mar_p, p_mar)

          # Main effect test
          beta0 <- beta[p_mar,, drop = FALSE]
          # The Wald-test statistic is the square of the coefficient divided by the variance of that coefficient
          w0 <- beta0^2 / diag(CovBeta)[p_mar]
          # Under the NULL of beta = 0, the w statistic is chi-square distributed with 1 df
          p_main <- pchisq(q = w0, df = 1, lower.tail = FALSE)

          ## Interaction test
          if (length(p_mar1) != 0) {
            # L matrix
            L <- diag(length(p_mar1))
            # L <- matrix(1, nrow = 1, ncol = length(p_mar_qtlxe))
            # Marker x E betas
            beta1 <- L %*% beta[p_mar1,, drop = FALSE]
            mat <- qr.solve(L %*% CovBeta[p_mar1, p_mar1] %*% t(L))

            # Wald statistic
            w1 <- t(beta1) %*% mat %*% beta1
            # Under the NULL of beta_1 = beta_2 = ... = beta_j = 0, the w statistic is chi-square distributed with j - 1 df
            p_qxe <- pchisq(q = as.numeric(w1), df = length(p_mar1), lower.tail = FALSE)

          } else {
            w1 <- p_qxe <- beta1 <- NA

          }
        } else {
          w1 <- p_qxe <- p_main <- w0 <- beta1 <- beta0 <- NA
          p_mar1 <- numeric()

        }

        # Output a list
        list(score = data.frame(term = c("main_effect", "qxe"), df = c(1, length(p_mar1)),
                                W_statistic = c(w0, w1), p_value = c(p_main, p_qxe)),
             estimates = list(beta0 = beta0, beta1 = beta1))

      }) }) # Close the apply and pmap functions

  } else {

    # Map over the list of markers and the X or K matrix
    scores <- pmap(list(mar_list, X_use, K_use, Hinv_use, K1_use), .f = function(mar, x, k, h, k1) {

      # Subset the marker matrix for those markers
      m <- M_test[,mar,drop = FALSE]
      # Designate the inverse phenotype vcov matrix
      H2inv <- h

      ## Create an incidence matrix of SNP genotypes per phenotypic observation
      # Model matrix of SNP main effect
      X_snp_main <- Z_rand %*% m

      # Apply function over the column of the M matrix (i.e. markers)
      apply(X = X_snp_main, MARGIN = 2, FUN = function(snp_main) {

        # Should QxE be tested?
        if (!test_qxe) {
          X_use1 <- cbind(x, snp_main)

        } else {
          # Model matrix of SNP x Environment
          X_snp_qxe <- c(snp_main) * X_fixed
          # Include main effect of SNP, then remove the last QTLxE term to keep model full rank
          X_use1 <- cbind(x, snp_main, X_snp_qxe[,-ncol(X_snp_qxe)])

        }

        ## Re-estimate variance components?
        if (!P3D) {
          fit <- mmer(Y = y, X = X_use1, Z = list(gen = list(Z = Z_use, K = k),
                                                  gxe = list(Z = Z1_use, K = k1)),
                      silent = T)

          # Calculate the weights (it is simply the ratio of the variance components divided
          # by the sum of the variance components)
          weights <- unlist(fit$var.comp[1:2]) %>% {. / sum(.)}

          # Calculate the sum of the variance components
          Vu <- sum(unlist(fit$var.comp[1:2]))

          H2inv <- H_inv(Vu = Vu, Ve = fit$var.comp$units[1], weights = weights,
                        Zlist = list(Z_use, Z1_use), Klist = list(k, k1))
        }

        # W matrix and inverse
        W <- crossprod(X_use1, H2inv %*% X_use1)
        Winv <- try(solve(W), silent = TRUE)

        # If the inverse is successful, calculate the p-value for that SNP using
        # an F test
        if (class(Winv) != "try-error") {
          # Vector of fixed effect coefficients
          beta <- Winv %*% crossprod(X_use1, H2inv %*% y)
          # Predictions (ignore random effects)
          y_hat <- X_use1 %*% beta
          # Residuals
          resid <- y - y_hat

          # Number of fixed terms
          p <- ncol(X_use1)
          v2 <- n - p
          # Estimate of sum of squared residuals
          s2 <- as.double(crossprod(resid, H2inv %*% resid))/v2
          # VCOV of fixed effect coefficients
          CovBeta <- s2 * Winv

          # Location of marker effects
          mar_p <- seq(ncol(x) + 1, p)
          # Location of marker main effect terms
          p_mar <- mar_p[1]
          p_mar1 <- setdiff(mar_p, p_mar)

          # Main effect test
          beta0 <- beta[p_mar,, drop = FALSE]
          # The Wald-test statistic is the square of the coefficient divided by the variance of that coefficient
          w0 <- beta0^2 / diag(CovBeta)[p_mar]
          # Under the NULL of beta = 0, the w statistic is chi-square distributed with 1 df
          p_main <- pchisq(q = w0, df = 1, lower.tail = FALSE)

          ## Interaction test
          if (length(p_mar1) != 0) {
            # L matrix
            L <- diag(length(p_mar1))
            # L <- matrix(1, nrow = 1, ncol = length(p_mar_qtlxe))
            # Marker x E betas
            beta1 <- L %*% beta[p_mar1,, drop = FALSE]
            row.names(beta1) <- colnames(X_use1)[p_mar1]
            mat <- qr.solve(L %*% CovBeta[p_mar1, p_mar1] %*% t(L))

            # Wald statistic
            w1 <- t(beta1) %*% mat %*% beta1
            # Under the NULL of beta_1 = beta_2 = ... = beta_j = 0, the w statistic is chi-square distributed with j - 1 df
            p_qxe <- pchisq(q = as.numeric(w1), df = length(p_mar1), lower.tail = FALSE)

          } else {
            w1 <- p_qxe <- NA

          }
        } else {
          w1 <- p_qxe <- p_main <- w0 <- beta1 <- beta0 <- NA
          p_mar1 <- numeric()
        }

        # Output a list
        list(score = data.frame(term = c("main_effect", "qxe"), df = c(1, length(p_mar1)),
                                W_statistic = c(w0, w1), p_value = c(p_main, p_qxe)),
             estimates = list(beta0 = beta0, beta1 = beta1))

        })

      }) # Close the apply and pmap functions

  }

  # Transpose the list
  scores_transp <- map(scores, transpose)

  # Combine the data.frames and return
  scores_df <- scores_transp %>%
    map("score") %>%
    map_df(~mutate(bind_rows(.), marker = rep(names(.), each = 2))) %>%
    as_data_frame()

  estimate_df <- scores_transp %>%
    map_df(~as_data_frame(.$estimates) %>%
             mutate(estimate = c("beta0", "beta1")) %>%
             gather(marker, value, -estimate))

  bind_cols(scores_df, estimate_df) %>%
    select(marker, term, estimate = value, df, W_statistic, p_value)

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
