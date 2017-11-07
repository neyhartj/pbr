#' Genomewide Assocation Analysis Through Linear Mixed-Models
#'
#' @description
#' Performs a genomewide association analysis for the main effect of QTL. Most
#' arguments are taken from the \code{\link[GWAS]{rrBLUP}} function.
#'
#' @param pheno A data.frame of phenotypic data.
#' @param geno A data.frame of marker names, positions, and genotypes.
#' @param fixed A formula definining the fixed effects in the \code{pheno} input.
#' @param model The type of association model to run. Can be one of: \code{"simple"},
#' \code{"K"}, \code{"Q"}, \code{"QK"}, \code{"G"}, \code{"QG"}, ... . See \emph{Details}
#' for more information on these models.
#' @param n.PC The number of principal components from singular value decomposition
#' of the \code{K} matrix to use to correct for population structure. If \code{0},
#' no correction for population structure is performed.
#' @param P3D Logical. Population parameters previous determined.
#' @param n.core The number of cores to use when calculating marker scores.
#' @param impute.method The method by which to impute the marker genotypes.
#'
#' @return
#' An object of class \code{gwas}, with marker scores for each trait.
#'
#' @details
#' A brief description of model types follows:
#'
#' \describe{
#'   \item{simple}{Tests for marker effect without correcting for background polygenic
#'   effect or for population structure. \cr
#'   Model: \eqn{y = X\beta + S\alpha + e}}
#'   \item{K}{Tests for marker effect while correcting for background polygenic
#'   effect, but not population structure. \cr
#'   Model: \eqn{y = X\beta + Zu + S\alpha + e}, where \eqn{u ~ N(0, K\sigma^2_G)}}
#'   \item{Q}{Tests for marker effect while correcting for population structure,
#'   but not background polygenic effect. \cr
#'   Model: \eqn{y = X\beta + Qv + S\alpha + e}}
#'   \item{QK}{Tests for marker effect while correcting for background polygenic
#'   effect and for population structure. \cr
#'   Model: \eqn{y = X\beta + Zu + Qv + S\alpha + e}, where \eqn{u ~ N(0, K\sigma^2_G)}}
#'   \item{G}{Tests for marker effect while correcting for the background polygenic
#'   effect of chromosomes other than the one on which the marker resides. Does
#'   not correct for population structure. \cr
#'   Model: \eqn{y = X\beta + Zu + S\alpha + e}, \cr where \eqn{u ~ N(0, K\sigma^2_G)} and K
#'   is the relationship matrix obtained using markers on all chromosomes excluding
#'   the one on which the putative QTL resides.}
#'   \item{QG}{Tests for marker effect while correcting for the background polygenic
#'   effect of chromosomes other than the one on which the marker resides. Does
#'   correct for population structure. \cr
#'   Model: \eqn{y = X\beta + Zu + Qv + S\alpha + e}, \cr where \eqn{u ~ N(0, K\sigma^2_G)} and K
#'   is the relationship matrix obtained using markers on all chromosomes excluding
#'   the one on which the putative QTL resides.}
#' }
#'
#' @return
#' An object of class \code{gwas}, with marker scores for each trait.
#'
#' @examples
#' data("tr_cap_genos_hmp")
#' data("tr_cap_phenos_met")
#'
#' # Filter the genotypes
#' geno <- tr_cap_genos_hmp %>%
#'   select(which(names(.) %in% c("marker", "chrom", "pos", unique(tr_cap_phenos_met$line))))
#'
#' pheno <- tr_cap_phenos_met
#'
#' \dontrun{}
#' # Perform association using the simple model versus the Q+K model
#' gwas_simple_out <- gwas(pheno = pheno, geno = geno, fixed = ~ trial, model = "simple")
#' gwas_QK_out <- gwas(pheno = pheno, geno = geno, fixed = ~ trial, model = "QK", n.PC = 2)
#'
#' # Compare with rrBLUP GWAS function for the Q+K model
#'
#'
#' # Use a subset of the data to fit all models
#' # Use the first 25 markers of the first two chromosomes
#' geno <- geno %>%
#'   filter(chrom %in% unique(chrom)[1:2]) %>%
#'   group_by(chrom)
#'
#' # Use phenotypes from the first two trials
#' pheno <- pheno %>%
#'   filter(trial %in% unique(trial)[1:2])
#'
#' models <- c("simple", "K", "Q", "QK", "G", "QG")
#'
#' gwas_out <- map(models, ~gwas(pheno = pheno, geno = geno, fixed = ~trial, model = ., n.PC = 2))
#'
#'
#' @importFrom rrBLUP A.mat
#' @importFrom EMMREML emmreml
#' @import purrr
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import parallel
#'
#' @export
#'
gwas <- function(pheno, geno, fixed = NULL, model = c("simple", "K", "Q", "QK", "G", "QG"),
                 n.PC = 0, P3D = TRUE, n.core = 1, impute.method = c("mean", "EM", "pass")) {

  ## ERROR
  pheno <- as.data.frame(pheno)
  geno <- as.data.frame(geno)

  # Rename the first three columns of the 'geno' input
  colnames(geno)[1:3] <- c("marker", "chrom", "pos")

  # n.PC cannot be less than 0
  stopifnot(n.PC >= 0)

  # Make sure 'fixed' is a formula
  stopifnot(class(fixed) == "formula")
  # Extract the terms in the fixed formula
  fixed_terms <- attr(terms(fixed), "term.labels")

  # Make sure the fixed effect columns are in the phenotype df
  stopifnot(all(fixed_terms %in% colnames(pheno)))

  # Match the imputation method
  model <- match.arg(model)
  impute.method <- match.arg(impute.method)

  # Number of phenotype columns
  n_pheno <- ncol(pheno) - 1

  # Name of the random effects column
  rand_name <- colnames(pheno)[1]

  # Convert the random effects (lines) and fixed effects (environment) to a factor
  pheno[,rand_name] <- as.factor(pheno[,rand_name, drop = TRUE])
  pheno[,fixed_terms] <- as.factor(pheno[,fixed_terms, drop = TRUE])

  # Number of fixed terms
  n_fixed <- length(fixed_terms)
  # Number of traits
  n_trait <- n_pheno - n_fixed
  # Find the names of the traits
  trait_names <- setdiff(x = tail(colnames(pheno), n = -1), y = fixed_terms)

  # Make sure all phenotyped lines are genotyped
  geno_names <- levels(pheno[,rand_name])
  stopifnot(all(geno_names %in% colnames(geno)))


  ## Internal functions
  ## Make a matrix full rank
  make_full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }

  ## M genotype matrix - subset for the levels of genotypes in the pheno df
  M <- t(subset(geno, select = which(names(geno) %in% geno_names)))
  dimnames(M) <- list(row.names(M), geno[, 1, drop = TRUE])

  # Impute markers
  geno_impute <- A.mat(X = M, min.MAF = 0, max.missing = 1, impute.method = impute.method,
                       shrink = FALSE, return.imputed = TRUE)

  # Extract the imputed markers
  M <- geno_impute$imputed

  ## SNP info (map, name, etc)
  snp_info <- subset(geno, select = 1:3)
  # Number of markers
  m <- ncol(M)
  # Marker names
  mar_names <- colnames(M)

  # Re-order the marker matrix
  M1 <- M[geno_names,]


  if (model %in% c("K", "Q", "QK")) {

    ## Create covariance matrices - only if the specified model is given
    # Extract the whole matrix from the imputed markers
    K_all <- geno_impute$A

    # Reorder the K matrices
    K0 <- list(K_all[geno_names, geno_names])


  } else if (model %in% c("G", "QG")) {

    # If only one chromosome is present, error
    stopifnot(n_distinct(snp_info$chrom) > 1)

    ## Create relationship matrices for markers on each chromsome, if called
    # First get the markers names per chromosome, and then collect the marker names
    # for all chromosomes except the ith chromosome
    mar_names_chr <- snp_info %>%
      split(.$chrom) %>%
      map(~.$marker) %>%
      map(~setdiff(mar_names, .))
    # Extract the genotypes for those markers and calculate relationship matrices
    K_chr <- mar_names_chr %>%
      map(~M[,.]) %>%
      map(~A.mat(., min.MAF = 0, max.missing = 1, shrink = FALSE))

    K0_chr <- K_chr %>%
      map(~.[geno_names, geno_names])

  }

  ## PCs for population structure
  # PCs for population structure
  if (n.PC > 0) {
    # If n.PC > 0, but the model is not Q, QK, or QK, do not use PCs
    if (!str_detect(model, "Q")) {
      warning("The specified model is not one of 'Q', 'QK', or 'QG'. No PCs will be used as covariates.")

    } else if (model == "QG") {
      # If the model is the QG model, model population structure per chromosome
      eig_vec_list <- map(K_chr, ~eigen(.)$vector)

    } else {
      eig_vec_list <- list(eigen(K_all)$vectors)

    }
  } else {
    # If n.PC is equal to 0, but the model contains 'Q', error out
    if (model %in% c("Q", "QK", "QG")) stop("'n.PC' == 0, but the model is one of 'Q', 'QK', or 'QG'.")

  }

  # Create a list to store output of multiple traits
  trait_scores <- vector("list", n_trait) %>%
    setNames(trait_names)

  # Iterate over phenotypes
  for (i in seq_along(trait_names)) {
    # Notification
    cat(paste("GWAS for trait: ", trait_names[i], ", using model: ", model, "\n", sep = ""))

    # Matrix for random effects
    ## First formula for the random effects
    rand_form <- as.formula(paste(trait_names[i], paste0("~ -1 +", rand_name)))
    mf <- model.frame(rand_form, pheno, drop.unused.levels = FALSE, na.action = "na.omit")
    Z0 <- model.matrix(rand_form, mf)

    ## Fixed effect model matrices
    # Formula for the response and fixed
    fixed_form <- as.formula(paste(trait_names[i], paste0("~ ", as.character(fixed)[-1], "- 1", collapse = "+")))
    mf <- model.frame(fixed_form, pheno, drop.unused.levels = TRUE, na.action = "na.omit")
    # Vector of the grand mean
    X_mu <- model.matrix(~ 1, mf)

    # If n_fixed is not greater than 0, only create a mean vector
    if (n_fixed > 0) {
      X_fixed <- model.matrix(fixed_form, mf)

    } else {
      # For completeness, include a matrix of 0s
      X_fixed <- matrix(data = 0, nrow = nrow(X_mu), ncol = 1)

    }

    # If population structure should be corrected via PC, add those vectors to the X matrix
    if (n.PC > 1) {
      # If a G model, create n_chrom Q matrices
      if (model == "QG") {
        Q_chr <- map(eig_vec_list, ~ Z0 %*% .[,seq(n.PC), drop = FALSE])

      } else if (model == "Q") {
        Q <- Z0 %*% eig_vec_list[[1]][,seq(n.PC), drop = FALSE]

      } else {
        Q <- NULL

      }
    } else {
      Q <- NULL

    }

    # Combine the Q matrix to make the fixed effect matrix
    if (model == "QG") {
      X_model <- map(Q_chr, ~make_full(cbind(X_mu, X_fixed, .)))

    } else {
      # Combine with the mean vector
      X_model <- list(make_full(cbind(X_mu, X_fixed, Q)))

    }

    ## vector of response
    y <- model.response(mf)
    # Number of obs
    n <- length(y)

    # Solve
    if (P3D) {

      # Fit the model
      # The stream should split between simple and Q - and other models
      if (model %in% c("simple", "Q")) {
        Z_model <- diag(length(y))
        K_model <- list(diag(ncol(Z_model)))

      } else if (model %in% c("K", "QK")) {
        Z_model <- Z0
        K_model <- K0

      } else {
        Z_model <- Z0
        K_model <- K0_chr

      }

      fit <- pmap(list(X_model, K_model), ~emmreml(y = y, X = .x, Z = Z_model, K = .y))
      Hinv <- pmap(list(fit, K_model),
                   ~H_inv(Vu = .x$Vu, Ve = .x$Ve, n = n, Z = Z_model, K = .y)) %>%
        set_names(names(K_model))

      cat("Variance components estimated. Testing markers.\n")

    } else {
      Hinv <- NULL

    }

    ## Calculate the scores per marker
    ## Parallelize if possible
    if (n.core > 1) {

      # How many cores are there to use?
      avail_cores <- detectCores()

      # Split the marker data
      mar_split <- snp_info %>%
        mutate(core = sort(rep(seq(n.core), length.out = nrow(.)))) %>%
        split(.$core)

      stopifnot(n.core <= avail_cores)

      # Detect the operating system
      if (str_detect(string = sessionInfo()$running, "Windows")) {

        # Create the cluster
        cluster <- makeCluster(n.core)
        # Extract functions to send to the clusters
        env_to_export <- list(as.list(loadNamespace("pbr")), as.list(loadNamespace("purrr")),
                      as.list(loadNamespace("EMMREML")),
                      list(M1 = M1, P3D = P3D, Hinv = Hinv, X_model = X_model,
                      Z_model = Z_model, y = y, K_model = K_model, model = model,
                      Z0 = Z0))
        env_to_export <- unlist(env_to_export, recursive = FALSE)

        # Export the functions and data
        clusterExport(cl = cluster,
                      varlist = c("M1", "P3D", "Hinv", "X_model", "y", "Z_model",
                               "K_model", "model", "Z0", "H_inv", "score_calc",
                               "pmap", "map"), envir = as.environment(env_to_export))

              # Run the code
        scores <- clusterApply(cl = cluster, x = mar_split, fun = function(mar_df) {
          score_calc(M_test = M1[,mar_df$marker, drop = FALSE], snp_info = mar_df,
                     P3D = P3D, Hinv = Hinv, X = X_model, y = y, Z0 = Z_model,
                     K0 = K_model, model = model, Z_rand = Z0) })

        # Shut down the cluster
        stopCluster(cluster)

      } else {
        # Otherwise use regular parallel
        scores <- mclapply(X = mar_split, FUN = function(mar_df) {
          score_calc(M_test = M1[,mar_df$marker, drop = FALSE], snp_info = mar_df,
                     P3D = P3D, Hinv = Hinv, X = X_model, y = y, Z0 = Z_model,
                     K0 = K_model, model = model, Z_rand = Z0)

        }, mc.cores = n.core)

      }

      # Collapse the list
      scores <- do.call("rbind", scores)

    } else {
      scores <- score_calc(M_test = M1, snp_info = snp_info, P3D = P3D, Hinv = Hinv,
                           X = X_model, y = y, Z0 = Z_model, K0 = K_model, model = model,
                           Z_rand = Z0)

    }

    # Add the marker scores to the matrix
    trait_scores[[trait_names[i]]] <- scores

  } # Close the trait loop

  # Edit the trait scores
  scores_df <- pmap(list(trait_scores, trait_names), ~mutate(.x, trait = .y)) %>%
    bind_rows() %>%
    select(trait, marker, names(.)) %>%
    left_join(snp_info, ., by = "marker") %>%
    arrange(trait)

  # Create a list with metadata and output as a gwas file
  out <- list(
    scores = scores_df,
    metadata = list(fixed = fixed, model = model, P3D = P3D, n.PC = n.PC,
                    impute.method = impute.method)
  )

  return(structure(out, class = "gwas"))

} # Close the function










