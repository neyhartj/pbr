## Heritability simulation
##
## Let's say we have 3 environments in which an augmented design is planted. There are 100 test genotypes and 10 check genotypes, each ## replicated 3 times. We need to figure out how to estimate heritability properly

# Libraries
library(tidyverse)
library(stringr)
library(broom)
library(regress)
library(lme4)
library(pbr)

# How many of each thing?
n_env <- 3
n_test <- 100
n_check <- 10
n_check_reps <- 3

h2_j <- c(0.7, 0.8, 0.5)
varR <- (1 / h2_j) - 1
names(varR) <- str_c("E", seq(n_env))

# Simulate genotypic values
test_geno <- rnorm(n = n_test, mean = 0, sd = 1)
check_geno <- rnorm(n = n_check, mean = 0, sd = 1)

# Combine
geno <- data_frame(line = c(str_c("T", seq(n_test)), rep(str_c("C", seq(n_check)), each = n_check_reps)),
                   rep = c(rep("R1", n_test), rep(str_c("R", seq(n_check_reps)), n_check)),
                   geno = c(test_geno, rep(check_geno, each = n_check_reps)))

# Add random deviations
p_obs <- bind_cols(geno, map_df(sqrt(varR), ~rnorm(n = nrow(geno), mean = 0, sd = .))) %>%
  gather(env, res, -line, -rep, -geno) %>%
  mutate(p_obs = geno + res)

# Adjust line names to denote checks vs test
p_obs1 <- p_obs %>%
  mutate(test = ifelse(str_detect(line, "^T"), line, "00check"),
         check = ifelse(str_detect(line, "^C"), line, "00test"),
         te = interaction(test, env),
         ce = interaction(check, env),
         re = interaction(rep, env)) %>%
  mutate_if(is.character, as.factor)



# The null model will be fit across environments
full_mod <- lmer(formula = p_obs ~ check + env + (1|test) + (1|te) + (1|ce) + (1|re),
                    data = p_obs1)

var_comp_full <- as_data_frame(VarCorr(full_mod)) %>%
  select(grp, vcov) %>%
  spread(grp, vcov) %>%
  unlist()

H_full <- var_comp_full["test"] / (var_comp_full["test"] + (var_comp_full["te"] / 3) + (var_comp_full["Residual"]))

## Two stage
##
## Stage one - estimate BLUEs

stage_one <- p_obs1 %>%
  group_by(env) %>%
  mutate(r = pbr::harm_mean(table(line))) %>%
  do({

    fit <- lm(p_obs ~ test + check, data = .)

    coef_test <- coef(fit) %>%
      data_frame(term = names(.), coef = .) %>%
      mutate(coef = coef + coef[1]) %>%
      filter(term != "(Intercept)", !str_detect(term, "^check"))

    V_j <- vcov(fit)[coef_test$term, coef_test$term]

    # Change the line term names
    coef_test <- coef_test %>%
      mutate(term = str_replace_all(term, "line", ""))

    line_terms <- pull(coef_test, term) %>%
      str_replace_all(pattern = "line", "")

    # Replace the names in the vcov matrix
    dimnames(V_j) <- list(line_terms, line_terms)

    # Extract the residual variance
    V_R <- sigma(fit) ^ 2

    # Return the coefficients and the variance of the adjusted means
    data_frame(fit = list(fit), BLUE = list(coef_test), V_j = list(V_j), V_R = V_R, r = .$r[[1]]) })


# For each environment, fit a new model with the line BLUEs and calculate heritability based on the measured residual variance
stage_one_herit <- stage_one %>%
  ungroup() %>%
  mutate(herit1 = var(BLUE[[1]]$coef) / (var(BLUE[[1]]$coef) + (V_R / r)),
         herit2 = var(BLUE[[1]]$coef) / (var(BLUE[[1]]$coef) + (V_R / 1)),
         herit1_bias = abs(herit1 - h2_j),
         herit2_bias = abs(herit2 - h2_j))

# Now find the heritability across environments
stage_one_data <- stage_one_herit %>%
  select(env, BLUE) %>%
  unnest() %>%
  mutate_at(vars(env:term), as.factor) %>%
  mutate(te = interaction(term, env))

R <- stage_one_herit %>%
  pull(V_j) %>%
  map(diag) %>%
  map(diag) %>%
  reduce(., Matrix::bdiag) %>%
  as.matrix()

stage_mod <- regress(formula = coef ~ env,
                     Vformula = ~term + te + R,
                     data = stage_one_data, identity = FALSE, pos = rep(T, 3))

var_comp <- stage_mod$sigma

H <- var_comp["term"] / (var_comp["term"] + (var_comp["te"] / n_env) + (var_comp["R"] / 1))


### Try to calculate heritability as defined by Piepho and Mohring 2007

oats <- john.alpha

## Fit the model - BLUEs
fit1 <- lmer(yield ~ 1 + rep + gen +  (1|rep:block), oats)

fit1_V <- diag(vcov(fit1))[-(1:3)]
vblue <- mean(fit1_V)


## BLUPs
fit2 <- lmer(yield ~ 1 + rep + (1|gen) +  (1|rep:block), oats)

blups <- ranef(fit2, condVar = T)$gen
vblup <- 2 * mean(attr(blups, "postVar"))

varG <- c(VarCorr(fit2)[["gen"]])


## Herit
# v1
(h2 <- varG / (varG + (vblue / 2)))

# V2
(h2.2 <- 1 - vblup / (2 * varG))





