#' Mixed Model Finlay-Wilkinson Regression
#'
#' @description
#' Fits Finlay-Wilkinson using mixed models
#'
#' @param y The vector of the response variables
#' @param GEN Vector of genotype names
#' @param ENV Vector of environmental names
#'
#' @examples
#' # Load data
#' data("gauch.soy")
#'
#' fw_out <- fw_mixed(y = gauch.soy$yield, GEN = gauch.soy$gen,
#'                    ENV = gauch.soy$env)
#'
#' @import lme4
#' @importFrom purrr map map_dbl
#'
#' @export
#'
fw_mixed <- function(formula, data, subset, weights) {

  # Convert genotypes and environments to factors
  GEN <- factor(GEN)
  ENV <- factor(ENV)

  # Indicators for GEN and ENV
  G_ID <- as.numeric(GEN)
  E_ID <- as.numeric(ENV)

  # Levels for GEN and ENV
  GEN_levels <- levels(GEN)
  ENV_levels <- levels(ENV)

  # Number of GEN and ENV
  n_GEN <- length(GEN_levels)
  n_ENV <- length(ENV_levels)

  # Contrasts
  fGENc <- GEN
  attr(fGENc, "contrasts") <- contr.sum(n_GEN)
  fENVc <- ENV
  attr(fENVc, "contrasts") <- contr.sum(n_ENV)

  # Model frame
  mf <- model.frame(y ~ fENVc + fGENc)
  # Base model
  lmer0 <- lmer(formula = y ~ (1|fGENc) + (1|fENVc), data = mf)
  # Random effects
  blups <- ranef(lmer0)

  # BLUPs of genotype and environmental effects
  h <- as.matrix(blups$fENVc) %>%
    as.numeric() %>%
    setNames(ENV_levels)
  h <- c(h, -sum(h, na.rm = T))

  # Add the environmental effect to the model frame
  mf2 <- cbind(mf, h = h[E_ID])

  # Control for lmer
  lmer_control <- lmerControl(check.nlev.gtr.1 = "ignore", )

  # Split the data.frame and apply the model
  suppressWarnings(fits2 <- mf2 %>%
    split(.$fGENc) %>%
    map(~ lmer(y ~ (0 + h | fGENc), data = ., control = lmer_control)))

  # Extract the coefficients (genotype mean)
  g <- fits2 %>%
    map_dbl(~ fixef(.))
  b <- fits2 %>%
    map_dbl(~ as.numeric(ranef(.)[[1]]))
  var_e <- fits2 %>%
    map_dbl(~ sigma(.)^2)
  df <- fits2 %>%
    map_dbl(~ df.residual(.))

  var_e_weighted <- sum(var_e * df) / sum(df)

  # Create list
  outlist <- list(y = y, GEN = GEN, ENV = ENV, GEN_levels = GEN_levels,
                  ENV_levels = ENV_levels, mu = 0, g = as.data.frame(g),
                  b = as.data.frame(b), var_e = as.data.frame(var_e),
                  var_e_weighted = var_e_weighted)

  return(outlist)

} # Close function
