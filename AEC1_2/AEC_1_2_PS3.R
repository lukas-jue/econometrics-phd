# Advanced Econometrics I - Part 2
# Problem Set 3 - Due 23rd January 2020

library(mvtnorm)
library(bbmle)
library(dplyr)
library(mdscore)

set.seed(20200118)

n <- 500
means <- c(0, 0)
cov_mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
multnorm <- rmvnorm(n, mean = means, sigma = cov_mat)
colnames(multnorm) <- c("x", "z")
x <- multnorm[,1]
z <- multnorm[,2]
epsilon <- rnorm(n)

plot(x, z)

y = 1 + 1.2*x + 1.5*z + 2*epsilon

# 1. Benchmark OLS

model_ols <- lm(y ~ x + z)
summary(model_ols)

ols_beta_hat <- model_ols$coefficients
ols_beta_se <- sqrt(diag(vcov(model_ols)))

# 2. Compute MLE Estimates


LL <- function(beta0, beta1, beta2, mu, sigma) {
  e = y - x * beta1 - z * beta2 - beta0
  #
  e = suppressWarnings(dnorm(e, mu, sigma, log = TRUE))
  #
  -sum(e)
}

model_mle <- mle2(LL, start = list(beta0 = 1, beta1 = 1, beta2 = 1, mu = 0, sigma = 1))
summary(model_mle)
mle_beta_hat <- model_mle@coef[1:3]
mle_beta_se <- sqrt(diag(vcov(model_mle)))

beta_hat_comp <- rbind(ols_beta_hat, mle_beta_hat)
se_comp <- rbind(ols_beta_se, mle_beta_se[1:3])
rownames(se_comp) <- c("ols_beta_se", "mle_beta_se")
sigma_hat <- se_comp^2
rownames(sigma_hat) <- c("ols_beta_sigma", "mle_beta_sigma")

# check error variances
epsilon_hat_ols <- y - x * ols_beta_hat["x"] - z * ols_beta_hat["z"] - ols_beta_hat["(Intercept)"]
epsilon_hat_mle <- y - x * mle_beta_hat["beta1"] - z * mle_beta_hat["beta2"] - mle_beta_hat["beta0"]

(epsilon_hat_ols %*% epsilon_hat_ols) / n 
(epsilon_hat_mle %*% epsilon_hat_mle) / n # mle has slightly higher variance -> due to numerical optimization based on guesses

var(epsilon_ols)
var(epsilon_mle) # here, error variances are identical: What's the difference to (e %*% e)/n?

# 3. Concentrated Log-Likelihood

sigma_hat <- (epsilon_hat_mle %*% epsilon_hat_mle) / n

LL_conc <- function(beta0, beta1, beta2){ # wrong! Need to fix!
  e = y - x * beta1 - z * beta2 - beta0
  -n/2 * (log(2*pi + log(e %*% e + 1)))
  }
optim(LL_conc)

LL_conc(beta0 = 1, beta1 = 1, beta2 = 1)


# 4. Restrict parameter
model_mle_rest <- mle2(LL, start = list(beta0 = 1, beta1 = 1, beta2 = 1, mu = 0, sigma = 1), fixed = list(beta2 = 1.3))
summary(model_mle_rest)
summary(model_mle)


# 5. Evaluate Null Hypothesis

## 5.1 Likelihood Ratio Test
  # Test ln(theta_hat) = ln(theta_hat_restricted)
  # If null is true, 
  # LR = -2(ln(L_r) - ln(L_u)) follows Chi-squ dist X(q)

LL_u <- model_mle@details$value
LL_r <- model_mle_rest@details$value

LR <- -2 * (LL_r - LL_u) # Shouldn't LL_u be smaller than LL_r?
anova(model_mle, model_mle_rest) # same result
# critical value at alpha = 0.05: Chi(1) = 3.841
# |LR| > 3.841 => H0 is not rejected


## 5.2 Wald Test
  # only based on unrestricted model
  # If restrictions are valid, c(theta_hat) - q is approx. 0
  # Compute W according to equation 37
var_beta2 <- model_mle@vcov["beta2", "beta2"]
W <- (model_mle@coef["beta2"] - 1.3)^2 / var_beta2
# W = 18.2, same critical value as above, so null is not rejected

## 5.3 Lagrange Multiplier Test
  # only based on restricted model

## 5.4 F Test
  # Note: In this case (one restriction) equal to Wald test
W <- (model_mle@coef["beta2"] - 1.3)^2 / var_beta2
