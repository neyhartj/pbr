#' Genomewide Association Analysis for Multiple Enviroments
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
#'
#' @importFrom rrBLUP A.mat
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
    eig.vec <- eigen(K)$vectors
  }

  # Create a matrix to store the marker scores for each trait
  marker_scores <- matrix(data = 0, nrow = m, ncol = n_trait * 2,
                          dimnames = list(mar_names, NULL))
  colnames(marker_scores) <- paste(rep(trait_names, each = 2), c("main", "qtlxe"), sep = "_")


  # Iterate over phenotypes
  for (i in seq_along(trait_names)) {
    # Notification
    print(paste("GWAS for trait:", trait_names[i]))

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

    # Combine with the mean vector
    X <- make_full(cbind(X_mu, X_fixed))

    # vector of response
    y <- model.response(mf)
    # Number of obs
    n <- length(y)

    # Matrix for random effects
    ## First formula for the random effects
    rand_form <- as.formula(paste(trait_names[i], paste0("~ -1 +", rand_name)))
    mf <- model.frame(rand_form, pheno, drop.unused.levels = TRUE, na.action = "na.omit")
    Z <- model.matrix(rand_form, mf)

    # Solve
    if (P3D) {
      Hinv <- mixed.solve(y, X = X, Z = Z, K = K, return.Hinv = TRUE)$Hinv
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
        score_calc(M_test = M[,markers, drop = FALSE], P3D = P3D, Hinv = Hinv,
                   X = X, y = y, X_fixed = X_fixed, Z = Z, K = K, qtlxe = TRUE)
      }, mc.cores = n.core)

      scores <- do.call("cbind", scores)

    } else {
      scores <- score_calc(M_test = M, P3D = P3D, Hinv = Hinv, X = X, y = y,
                           X_fixed = X_fixed, Z = Z, K = K, qtlxe = TRUE)

    }

    # Add the marker scores to the matrix
    j <- (i * 2 - 1)
    marker_scores[,c(j, j + 1)] <- scores

  } # Close the trait loop

  # Create a data.frame
  out <- data.frame(snp_info, marker_scores)
  class(out) <- c("gwas", class(out))

  return(out)

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
#' @importFrom rrBLUP mixed.solve
#'
#'
score_calc <- function(M_test, P3D, Hinv, y, X, X_fixed, Z, K, qtlxe = TRUE) {

  # Apply function over the column of the M matrix (i.e. markers)
  apply(X = M_test, MARGIN = 2, FUN = function(snp) {

    # Number of observations
    n <- length(y)

    ## First calculate the main effect of the marker
    # Model matrix of SNP main effect
    X_snp_main <- Z %*% snp
    X1 <- cbind(X, X_snp_main)

    ## Re-estimate variance components?
    if (!P3D) {
      H2inv <- mixed.solve(y = y, X = X1, Z = Z, K = K, return.Hinv = TRUE)$Hinv

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
      # Location of marker fixed terms
      p_mar <- p_mar <- seq(ncol(X) + 1, p)
      # Other
      v1 <- 1
      v2 <- n - p

      # Vector of fixed effects
      beta <- Winv %*% crossprod(X1, H2inv %*% y)
      # Residuals
      resid <- y - X1 %*% beta
      # Estimate of sum of squared residuals
      s2 <- as.double(crossprod(resid, H2inv %*% resid))/v2
      CovBeta <- s2 * Winv

      # F-statistic for the marker main effect
      Fstat <- beta[p_mar]^2 / diag(CovBeta)[p_mar]
      # Quantiles for beta distribution
      q <- v2/(v2 + v1 * Fstat)
      main_effect_pvalue <- unname(pbeta(q, v2/2, v1/2))

    } else {
      main_effect_pvalue <- NA

    }

    if (qtlxe) {

      ## Now calculate marker effects for each environment
      # Model matrix of SNP x Environment
      X_snp_qtle <- c(X_snp_main) * X_fixed
      X1 <- cbind(X, X_snp_qtle)

      ## Re-estimate variance components?
      if (!P3D) {
        H2inv <- mixed.solve(y = y, X = X1, Z = Z, K = K, return.Hinv = TRUE)$Hinv

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
        # Location of marker fixed terms
        p_mar <- p_mar <- seq(ncol(X) + 1, p)

        # Vector of fixed effects
        beta <- Winv %*% crossprod(X1, H2inv %*% y)
        # Residuals
        resid <- y - X1 %*% beta
        # Estimate of sum of squared residuals
        s2 <- as.double(crossprod(resid, H2inv %*% resid))/v2
        CovBeta <- s2 * Winv

        # L matrix
        L <- diag(length(p_mar))
        # Marker x E betas
        mar_beta <- L %*% beta[p_mar]
        mat <- qr.solve(L %*% CovBeta[p_mar, p_mar] %*% t(L))

        # Chi-square statistic
        chisqustat <- t(mar_beta) %*% mat %*% mar_beta
        # P- value
        qtlxe_pvalue <- as.numeric(pchisq(q = chisqustat, df = length(p_mar), lower.tail = FALSE))

      } else {
        qtlxe_pvalue <- NA

      }

      # Return the main effect and qtlxe p-values
      return(c(main_effect = main_effect_pvalue, qtlxe = qtlxe_pvalue))

    } else {
      return(c(main_effect = main_effect_pvalue))
    }
  })

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
