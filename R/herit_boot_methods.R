#'
#' @describeIn herit_boot
#'
herit_boot.lm <- function(object, geno.term, env.term = NULL, ge.term = NULL,
                          boot.reps = 1000, alpha = 0.95) {

  # Grab the model frame
  object_mf <- model.frame(object)

  # Grab the model formula
  object_form <- formula(object)

  # Sample bootstrap replicates
  boot_df <- bootstrap(data = object_mf, n = boot.reps)

  # Iterate over samples
  boot_herit <- map(boot_df$strap, ~ lm(object_form, data = .)) %>%
    map(herit, geno.term = geno.term, env.term = env.term, ge.term = ge.term)

  # Convert to data.frame
  boot_herit_df <- boot_herit %>%
    lapply("[[", 2) %>%
    data.frame()

  # Rename
  names(boot_herit_df) <- paste("rep", seq(ncol(boot_herit_df)), sep = "")

  # Mutate and gather
  boot_herit_df1 <- boot_herit_df %>%
    mutate(Term = boot_herit[[1]]$Term) %>%
    gather(rep, estimate, -Term)

  # Summarize and return
  boot_herit_df1 %>%
    group_by(Term) %>%
    summarize(mean = mean(estimate), ci_lower = quantile(estimate, 1 - alpha), ci_upper = quantile(estimate, alpha)) %>%
    as.data.frame()

} # Close function
#'
#' @describeIn herit_boot
#'
herit_boot.lmerMod <- function(object, geno.term, env.term = NULL, ge.term = NULL,
                               boot.reps = 1000, alpha = 0.95) {

  # Grab the model frame
  object_mf <- model.frame(object)

  # Grab the model formula
  object_form <- formula(object)

  # Sample bootstrap replicates
  boot_df <- bootstrap(data = object_mf, n = boot.reps)

  # Iterate over samples
  boot_herit <- map(boot_df$strap, ~ lmer(object_form, data = .)) %>%
    map(herit, geno.term = geno.term, env.term = env.term, ge.term = ge.term)

  # Convert to data.frame
  boot_herit_df <- boot_herit %>%
    lapply("[[", 2) %>%
    data.frame()

  # Rename
  names(boot_herit_df) <- paste("rep", seq(ncol(boot_herit_df)), sep = "")

  # Mutate and gather
  boot_herit_df1 <- boot_herit_df %>%
    mutate(Term = boot_herit[[1]]$Term) %>%
    gather(rep, estimate, -Term)

  # Summarize and return
  boot_herit_df1 %>%
    group_by(Term) %>%
    summarize(mean = mean(estimate), ci_lower = quantile(estimate, 1 - alpha), ci_upper = quantile(estimate, alpha)) %>%
    as.data.frame()

} # Close function


