#' Multi-Environment Genomewide Association Analysis
#'
#' @description
#' Performs a genomewide association analysis for the main effect of markers and
#' their interaction with the environment. Most arguments are taken from the
#' \code{\link[GWAS]{rrBLUP}} function.
#'
#' @param pheno A data.frame of phenotypic data.
#' @param geno A data.frame of marker names, positions, and genotypes.
#' @param fixed A character vector of fixed effects.
#' @param K The covariance matrix of random effects (i.e. genotypes).
#' @param n.PC The number of principal components from singular value decomposition
#' of the \code{K} matrix to use to correct for population structure.
#' @param P3D Logical. Population parameters previous determined.
#' @param n.core The number of cores to use when calculating marker scores.
#' @param impute.method The method by which to impute the marker genotypes.
#'
#' @return
#' An object of class \code{gwas}, with marker scores for each trait.
#'
#' @examples
#'
#' data("tr_cap_genos_hmp")
#' data("tr_cap_phenos_met")
#'
#' # Filter the genotypes
#' geno <- subset(tr_cap_genos_hmp,
#'                select = which(names(tr_cap_genos_hmp) %in%
#'                           c("marker", "chrom", "pos", unique(tr_cap_phenos_met$line))))
#'
#' pheno <- tr_cap_phenos_met
#'
#' \dontrun{}
#' # Perform association
#' gwas_out <- gwas_e(pheno = pheno, geno = geno, fixed = "trial")
#' gwas1 <- GWAS()
#'
#' @importFrom rrBLUP A.mat
#' @importFrom EMMREML emmreml
#' @import purrr
#' @import dplyr
#'
#' @export
#'
gwas_e <- function(pheno, geno, fixed = NULL, K = NULL, n.PC = 0, P3D = TRUE,
                   n.core = 1, impute.method = c("mean", "EM", "pass")) {

  ## ERROR
  pheno <- as.data.frame(pheno)
  geno <- as.data.frame(geno)

  # Match the imputation method
  impute.method <- match.arg(impute.method)

  # Number of phenotype columns
  n_pheno <- ncol(pheno) - 1

  # Make sure the fixed effect columns are in the phenotype df
  if (!all(fixed %in% colnames(pheno)))
    stop("The column name in 'fixed' is not in the 'pheno' data.frame.")

  # Number of fixed effects
  n_fixed <- length(fixed)
  # If n_fixed > 1, reject
  if (n_fixed > 1)
    stop("The number of fixed effects should only be 1 (environment).")

  # Name of the random effects column
  rand_name <- colnames(pheno)[1]

  # Convert the random effects (lines) and fixed effects (environment) to a factor
  pheno[,rand_name] <- as.factor(pheno[,rand_name, drop = TRUE])
  pheno[,fixed] <- as.factor(pheno[,fixed, drop = TRUE])

  # Number of traits
  n_trait <- n_pheno - n_fixed
  # Find the names of the traits
  trait_names <- setdiff(x = tail(colnames(pheno), n = -1), y = fixed)

  # Make sure all phenotyped lines are genotyped
  geno_names <- levels(pheno[,rand_name])

  if (!all(geno_names %in% colnames(geno)[-1:-3])) {
    stop("Some genotypes with phenoypic information were not genotyped.")
  }

  ## Internal functions
  ## Make a matrix full rank
  make_full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }

  ## M genotype matrix
  M <- t(subset(geno, select = -1:-3))
  dimnames(M) <- list(row.names(M), geno[, 1, drop = TRUE])

  # Impute markers
  geno_impute <- A.mat(X = M, min.MAF = 0, max.missing = 1, impute.method = impute.method,
                       shrink = FALSE, return.imputed = TRUE)

  M <- geno_impute$imputed


  # SNP info (map, name, etc)
  snp_info <- subset(geno, select = 1:3)
  # Number of markers
  m <- ncol(M)
  # Marker names
  mar_names <- colnames(M)

  ## K covariance matrix
  if (is.null(K)) {
    K <- geno_impute$A
  }
  # PCs for population structure
  if (n.PC > 0) {
    eig_vec <- eigen(K)$vectors
  }

  # Create a list to store output of multiple traits
  trait_scores <- vector("list", n_trait) %>%
    setNames(trait_names)


  # Iterate over phenotypes
  for (i in seq_along(trait_names)) {
    # Notification
    print(paste("GWAS for trait:", trait_names[i]))

    # Matrix for random effects
    ## First formula for the random effects
    rand_form <- as.formula(paste(trait_names[i], paste0("~ -1 +", rand_name)))
    mf <- model.frame(rand_form, pheno, drop.unused.levels = TRUE, na.action = "na.omit")
    Z0 <- model.matrix(rand_form, mf)

    # Re-order the K matrix
    K1 <- K[levels(mf$line), levels(mf$line)]

    # Re-order the marker matrix
    M1 <- M[levels(mf$line),]


    # ## Random effect of GxE
    # rand_form <- as.formula(paste(trait_names[i], paste0("~ -1 +", paste(rand_name, fixed, sep = ":"))))
    # mf <- model.frame(rand_form, pheno, drop.unused.levels = TRUE, na.action = "na.omit")
    # Z1 <- model.matrix(rand_form, mf)
    # # K Matrix of GxE
    # K_Z1 <- diag(ncol(Z1))


    # Formula for the response and fixed
    mf_form <- as.formula(paste(trait_names[i], paste0("~ ", fixed, "- 1", collapse = "+")))

    # Create the model frame
    # Drop NAs and unused factor levels
    mf <- model.frame(mf_form, pheno, drop.unused.levels = TRUE, na.action = "na.omit")

    ## Create model matrices for the fixed effects
    # First vector for the grand mean
    X_mu <- model.matrix(~ 1, mf)

    # If n_fixed is greater than 0, only create a mean vector
    if (n_fixed > 0) {
      X_fixed <- model.matrix(mf_form, mf)

    } else {
      # For completeness, include a matrix of 0s
      X_fixed <- matrix(data = 0, nrow = nrow(X_mu), ncol = 1)

    }

    # If population structure should be corrected via PC, add those vectors to the X matrix
    if (n.PC > 1) {
      X_pc <- Z0 %*% eig_vec[,seq(n.PC), drop = FALSE]

    } else {
      X_pc <- NULL

    }

    # Combine with the mean vector
    X <- make_full(cbind(X_mu, X_fixed, X_pc))

    # vector of response
    y <- model.response(mf)
    # Number of obs
    n <- length(y)

    # Solve
    if (P3D) {

      # Fit the mixed model
      fit <- emmreml(y = y, X = X, Z = Z0, K = K1)
      Hinv <- H_inv(Vu = fit$Vu, Ve = fit$Ve, n = n, Z = Z0, K = K1)

      # fit <- emmremlMultiKernel(y = y, X = X, Zlist = list(Z0, Z1), Klist = list(K1, K_Z1))

#       fit <- mixed.solve(y, X = X, Z = Z, K = K1, return.Hinv = TRUE)
#       Hinv <- fit$Hinv
      print("Variance components estimated. Testing markers.")

    } else {
      Hinv <- NULL

    }


    ## Calculate the scores per marker
    ## Parallelize if possible
    if (n.core > 1) {
      mar_split <- split(mar_names, cut(seq_along(mar_names), breaks = n.core, labels = FALSE))

      # Use mclapply
      scores <- mclapply(X = mar_split, FUN = function(markers) {
        score_calc(M_test = M1[,markers, drop = FALSE], P3D = P3D, Hinv = Hinv,
                   X = X, y = y, X_fixed = X_fixed, Z = Z0, K = K1, qtlxe = TRUE)
      }, mc.cores = n.core)

      # Collapse the list
      scores <- transpose(scores) %>%
        map(bind_rows)

    } else {
      scores <- score_calc(M_test = M1, P3D = P3D, Hinv = Hinv, X = X, y = y,
                           X_fixed = X_fixed, Z = Z0, K = K1, qtlxe = TRUE)

    }

    # Add the marker scores to the matrix
    trait_scores[[trait_names[i]]] <- scores

  } # Close the trait loop

  # Transpose the trait list
  trait_scores <- trait_scores %>%
    transpose() %>%
    map(bind_rows) %>%
    map(~mutate(., trait = rep(trait_names, each = nrow(.) / n_trait))) %>%
    map(~ select(., trait, marker, names(.)))

  # Create a data.frame
  class(trait_scores) <- c("gwas", class(trait_scores))

  return(trait_scores)

} # Close the function






#' Calculate marker scores in association study
#'
#' @param M_test Marker genotype matrix
#' @param P3D Logical - are variance components known or should they be estimated per marker?
#' @param Hinv The inverse of the H matrix (from the mixed-model output)
#' @param y Vector of response values
#' @param X Incidence matrix of fixed effects (must be full-rank)
#' @param X_fixed Incidence matrix of environment fixed effects (but not including intercept)
#' @param Z Incidence matrix of random effects (i.e. genotypes)
#' @param K The covariance matrix of random effects
#' @param qtlxe Logical - should QTLxE effects be tested?
#'
#' @details
#' Calculate the p-value for each marker in an association study. The main effect
#' of each marker is tested using an F-test, while QTLxE interaction is tested using
#' a Wald test.
#'
#' @importFrom EMMREML emmreml
#' @import dplyr
#'
#'
score_calc <- function(M_test, P3D, Hinv, y, X, X_fixed, Z, K, qtlxe = TRUE) {

  # Apply function over the column of the M matrix (i.e. markers)
  scores <- apply(X = M_test, MARGIN = 2, FUN = function(snp) {

    # number of observations
    n <- length(y)

    # Test only the main effect or fit main effect + GxE?
    if (!qtlxe) {

      ## First calculate the main effect of the marker
      # Model matrix of SNP main effect
      X_snp_main <- Z %*% snp
      X1 <- cbind(X, X_snp_main)

      ## First calculate the main effect of the marker
      # Model matrix of SNP main effect
      X_snp_main <- Z %*% snp
      X1 <- cbind(X, X_snp_main)

      ## Re-estimate variance components?
      if (!P3D) {
        # Fit the mixed model
        fit <- emmreml(y = y, X = X1, Z = Z, K = K)
        H2inv <- H_inv(Vu = fit$Vu, Ve = fit$Ve, n = n, Z = Z, K = K)

      } else {
        H2inv <- Hinv

      }

      # W matrix and inverse
      W <- crossprod(X1, H2inv %*% X1)
      Winv <- try(solve(W), silent = TRUE)

      # If the inverse is successful, calculate the p-value for that SNP using
      # an F test
      if (class(Winv) != "try-error") {
        # Number of fixed terms
        p <- ncol(X1)
        v2 <- n - p
        # Location of marker fixed terms
        p_mar <- seq(ncol(X) + 1, p)

        # Vector of fixed effect coefficients
        beta <- Winv %*% crossprod(X1, H2inv %*% y)
        # Residuals
        resid <- y - X1 %*% beta
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

        return(list(test = out_df, beta0_hat = NA, beta1_hat = NA))

      }



    } else {

      ## Now calculate marker effects for each environment
      # Model matrix of SNP main effect
      X_snp_main <- Z %*% snp
      # Model matrix of SNP x Environment
      X_snp_qtle <- c(X_snp_main) * X_fixed

      # Include main effect of SNP, then remove the last QTLxE term to keep model full rank
      X1 <- cbind(X, X_snp_main, X_snp_qtle[,-ncol(X_snp_qtle)])

      ## Re-estimate variance components?
      if (!P3D) {
        # Fit the mixed model
        fit <- emmreml(y = y, X = X1, Z = Z, K = K)
        H2inv <- H_inv(Vu = fit$Vu, Ve = fit$Ve, n = n, Z = Z, K = K)

      } else {
        H2inv <- Hinv

      }

      # W matrix and inverse
      W <- crossprod(X1, H2inv %*% X1)
      Winv <- try(solve(W), silent = TRUE)

      # If the inverse is successful, calculate the p-value for that SNP using
      # a Wald test
      if (class(Winv) != "try-error") {
        # Number of fixed terms
        p <- ncol(X1)
        v2 <- n - p
        # Location of marker main effect term
        p_mar <- ncol(X) + 1
        # Location of marker x env terms
        p_mar_qtlxe <- seq(ncol(X) + 2, p)

        # Vector of fixed effects
        beta <- Winv %*% crossprod(X1, H2inv %*% y)
        # Residuals
        resid <- y - X1 %*% beta
        # Estimate of sum of squared residuals
        s2 <- as.double(crossprod(resid, H2inv %*% resid))/v2
        CovBeta <- s2 * Winv

        ## Test main effect
        beta0 <- beta[p_mar,, drop = FALSE]
        w0 <- beta0^2 / diag(CovBeta)[p_mar]
        # Under the NULL of beta = 0, the w statistic is chi-square distributed with 1 df
        p_main <- pchisq(q = as.numeric(w0), df = 1, lower.tail = FALSE)


        ## Test QTLxE effects
        # L matrix
        L <- diag(length(p_mar_qtlxe))
        # L <- matrix(1, nrow = 1, ncol = length(p_mar_qtlxe))
        # Marker x E betas
        beta1 <- L %*% beta[p_mar_qtlxe]
        mat <- qr.solve(L %*% CovBeta[p_mar_qtlxe, p_mar_qtlxe] %*% t(L))

        # Wald statistic
        w1 <- t(beta1) %*% mat %*% beta1
        # Under the NULL of beta_1 = beta_2 = ... = beta_j = 0, the w statistic is chi-square distributed with j - 1 df
        p_qtlxe <- pchisq(q = as.numeric(w1), df = length(p_mar_qtlxe), lower.tail = FALSE)



        # Estimate of dropped environment
        beta1 <- rbind(beta[p_mar_qtlxe,,drop = FALSE] + c(beta0), beta0)
        # Add names to beta matrices
        dimnames(beta1) <- list(colnames(X_fixed), "beta1_hat")
        colnames(beta0) <- "beta0_hat"

        # Output data.frame
        out_df <- data.frame(term = c("main_effect", "qtl_x_env"),
                             df = c(1, length(p_mar_qtlxe)), W_statistic = c(w0, w1),
                             p_value = c(p_main, p_qtlxe), row.names = NULL, stringsAsFactors = FALSE)

        out_list <- list(test = out_df, beta0_hat = beta0, beta1_hat = beta1)

      } else {

        out_df <- data.frame(term = c("main_effect", "qtl_x_env"),
                             df = c(1, ncol(X1) - ncol(X) - 1), W_statistic = NA,
                             p_value = NA, row.names = NULL, stringsAsFactors = FALSE)

        return(list(test = out_df, beta0_hat = NA, beta1_hat = NA))

      }

    } })

  # Combine the data.frames
  test_df <- bind_rows(lapply(scores, "[[", "test"))
  test_df1 <- data.frame(marker = rep(names(scores), each = ifelse(qtlxe, 2, 1)),
                         test_df, stringsAsFactors = FALSE)

  # Extract main effects
  beta0_hat_df <- sapply(scores, "[[", "beta0_hat") %>%
    data.frame(marker = names(.), beta0_hat = ., row.names = NULL, stringsAsFactors = FALSE)

  if (qtlxe) {

    # Extract QTLxE effects
    beta1_hat_df <- lapply(scores, "[[", "beta1_hat") %>%
      do.call("cbind", .) %>% t() %>%
      data.frame(marker = names(scores), ., row.names = NULL, stringsAsFactors = FALSE) %>%
      setNames(c("marker", colnames(X_fixed)))

  } else {
    beta1_hat_df <- NULL

  }

  # Output list
  return(list(sig_test = test_df1, beta0_hat = beta0_hat_df, beta1_hat = beta1_hat_df))

} # Close the function



#' Plot the output from GWAS
#'
#' @param x An object of class \code{gwas}.
#' @param fdr.level The false discovery rate level.
#'
#' @importFrom tidyr gather
#' @importFrom dplyr mutate group_by ungroup
#' @import ggplot2
#'
#' @export
#'
plot.gwas <- function(x, fdr.level = 0.05) {

  ## Tidy up and adjust the p-values
  gwas_tidy <- gather(x, trait, p_value, -1:-3) %>%
    group_by(trait) %>%
    mutate(p_adj = p.adjust(p = p_value, method = "fdr"),
           neg_log_p_adj = -log10(p_adj),
           neg_log_fdr = -log10(fdr.level)) %>%
    ungroup()

  # Rename the first 3 columns
  colnames(gwas_tidy)[1:3] <- c("marker", "chrom", "pos")

  # Plot
  g <- gwas_tidy %>%
    ggplot(aes(x = pos, y = neg_log_p_adj, col = chrom)) +
    geom_point() +
    geom_hline(aes(yintercept = neg_log_fdr)) +
    facet_grid(trait ~ chrom, switch = "x") +
    ylab("-log10(p)") +
    xlab("Position") +
    scale_color_discrete(guide = FALSE)

  print(g)

} # Close the function


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


